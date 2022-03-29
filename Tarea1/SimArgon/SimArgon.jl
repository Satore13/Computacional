include("../Verlet.jl")
Frame, Particle = Verlet.BFrame, Verlet.BParticle
#= F(r) = 4ε/r * [12 (σ/ r)^12 - 6(σ / r)^6], tomamos σ = ε = 1

=#
using Serialization, CairoMakie
global L = 10.0 #Tamaño del receptáculo

function norm(v::Vector{N})::AbstractFloat where N <: Number
    sum = zero(N)
    for vi in v
        sum += vi^2
    end
    return sqrt(sum)
end

function update(f::Frame)
    #Creamos el array que contendrá las partículas nuevas
    new_particles = Particle[]

    #Función cuyo fin es wrappear en una coordenada a [0,L]
    function wrap_pos(x)
        return mod(x, L)
    end

    #Pasamos a través de todas las partículas y vamos guardando las actualizadas en el array
    for p in f.particles
        
        new_x = [wrap_pos(xi) for xi in p.x]
        push!(new_particles, Particle(new_x, p.v, p.mass))
    end
    return Frame(f.t, new_particles)
end

function get_shortest_distance(x::Vector{Float64}, other_x::Vector{Float64})::Vector{Float64}
    #Vectores de traslación a las 8 casillas circundantes
    traslaciones = [[0.0, 0.0], [0.0, L], [0.0, -L], [L, 0.0], [-L, 0.0], [L, L], [L, -L], [-L, L], [-L, -L]]

    #El vector más largo posible -por limitacione de periodicidad - es el que va del origen al origen enfrentado a él
    shortest_vector = [L, L]
    shortest_norm = norm(shortest_vector)

    #Comprobamos de todas las posiciones trasladadas cuál es la más próxima
    for trasl in traslaciones
        t_other_x = other_x + trasl
        current_v = x - t_other_x
        if norm(current_v) < shortest_norm
            shortest_norm = norm(current_v)
            shortest_vector = current_v
        end 
    end

    return shortest_vector
end

function get_shortest_distance(p::Particle, other_p::Particle)::Vector{Float64}
    get_shortest_distance(p.x, other_p.x)
end

function random_particle()::Particle
    mass = 1.0
    x = rand(Float64, 2) .* L
    v = [0.0, 0.0]
    return Particle(x, v, mass)
end

function acceleration(p, f)::Vector{Float64}
    mass = 1 # acc = force
    total_force = [0.0, 0.0]

    for other_p in f.particles
        if other_p === p
            continue
        end
        r = get_shortest_distance(p, other_p)
        norm_r = norm(r)
        if norm_r >= 3.0
            continue
        end
        total_force += 4*(12/norm_r^14 - 6/norm_r^8) .* r
    end
    return total_force
end

#Parametros de la simulación
const FPS = 30

function run()
    @show step = 0.002
    @show population_size = 10
    @show simulation_length = 50
    number_of_real_frames = simulation_length / step

    #Generamos las partículas en su estado inicial
    initial_particles = Particle[]
    for _ in 1:population_size
        push!(initial_particles, random_particle())
    end

    #Creamos el primer Frame con t=0.0
    current_frame = Frame(0.0, initial_particles)

    #Guardamos el primer Frame
    saved_frames = Frame[]
    push!(saved_frames, current_frame)

    #Ejecutamos las simulaciones
    for _ in 1:number_of_real_frames
        current_frame = Verlet.stepFrame(current_frame, step, acceleration, update)
        push!(saved_frames, current_frame)
    end

    #Y por último guardamos la simulación
    serialize("Tarea1/SimArgon/sim.out", (step, saved_frames))
end

function render()
    #Primero cargamos la simulación
    println("Cargando archivo...")
    (step, data) = deserialize("Tarea1/SimArgon/sim.out")
    println("Archivo cargado")
    @show FPS
    @show step
    #Creamos la figura y el axis donde dibujar
    fig = Figure()
    ax_main = Axis(fig[1,1])
    limits!(ax_main, (0, L), (0, L))

    #Contabilizador del frame actual
    current_frame_i = Observable{Int64}(1)

    #= Tenemos que crear un array de puntos que dependan de current_frame_i y nos den la coordenada en cada instante de las partículas =#
    points = @lift([Point2f(p.x[1], p.x[2]) for p in data[$current_frame_i].particles])
    scatter!(ax_main, points)


    number_of_frames = size(data)[1]
    println("Comenzado el renderizado...")
    last_percentage = 0.0

    #Renderizamos un frame de cada 1/(step * FPS) ya que esto nos da el ratio entre frames a renderizar y frame reales
    #Grabamos la simulación
    record(fig, "Tarea1/SimArgon/sim.mp4", range(1, number_of_frames, step = round(Int64, FPS^-1 / step)); framerate = FPS) do row_number
        current_frame_i[] = row_number

        current_percentage = (row_number/number_of_frames*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
    end

end