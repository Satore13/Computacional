include("../Verlet.jl")
Frame, Particle = Verlet.BFrame, Verlet.BParticle
using Serialization, CairoMakie
const FPS = 30

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
        if x == NaN
            println("Están desapareciendo partículas")
        end
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
    #Vectores de traslación a las 8 casillas circundantes incluyendo la identidad
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
    θ = rand(Float64)*2π
    v = 0.0 .*[cos(θ), sin(θ)]
    return Particle(x, v, mass)
end

function acceleration(p, f)::Vector{Float64}
    mass = 1 # acc = force
    total_force = [0.0, 0.0]

    for other_p in f.particles
        if other_p === p
            continue
        end
        total_force += force_between_p(p, other_p)
    end
    if total_force == NaN
        println("ERROR: Fuerza infinita")
    end
    return total_force

end

function force_between_p(p::Particle, other_p::Particle)::Vector{Float64}
    r = get_shortest_distance(p, other_p)
    norm_r = norm(r)

    
    if norm_r >= 3
        return [0.0, 0.0]
    end

    if norm_r <= 0.8
        
        norm_r = 0.8
        r ./= norm(r)
        r .*= 0.8
    end

    return 4*(12/norm_r^14 - 6/norm_r^8) .* r
end

function plot_force()
    ax = Axis(Figure()[1, 1])
    x = collect(range(0, L/2, length = 500))
    y = [force_between_p(Particle( [xi, 0.0], [0.0, 0.0], 1.0), Particle([0.0, 0.0], [0.0, 0.0], 1.0))[1] for xi in x]
    
    lines!(ax, x, y)
    current_figure()
end

function plot_potential()
    ax = Axis(Figure()[1, 1])
    x = collect(range(0, L/2, length = 500))
    y = [potential_energy(Particle( [xi, 0.0], [0.0, 0.0], 1.0), Particle([0.0, 0.0], [0.0, 0.0], 1.0))[1] for xi in x]
    
    lines!(ax, x, y)
    current_figure()
end

function potential_energy(p::Particle, other_p::Particle)::Float64
    r = get_shortest_distance(p, other_p)
    norm_r = norm(r)

    if other_p === p return 0 end #Una partícula no puede tener energía potencial consigo misma!
    
    if norm_r >= 3
        return 0.0
    end

    if norm_r <= 0.8
        #Potential is equal to potential accumulated till 0.8 plus lineal grow for constant force in 0-0.8 regime
        return  4*(1/0.8^12 - 1/0.8^6)  +  (4*(12/0.8^13 - 6/0.8^7) * (0.8 - norm_r))
    end

    return 4*(1/norm_r^12 - 1/norm_r^6)

end
function potential_energy(f::Frame)::Float64
    return 1/2 * sum([potential_energy(p, op) for p in f.particles, op in f.particles])
end

function kinetic_energy(p::Particle)
    return 0.5 * norm(p.v)^2
end
function kinetic_energy(f::Frame)
    return sum([kinetic_energy(p) for p in f.particles])
end

function plot_energies(filename::String)
    plot_energies(deserialize("Tarea1/SimArgon/" * filename)[1])
end
function plot_energies(data::Vector{Frame})
    fig = Figure()
    ax = Axis(fig[1, 1])

    t = [f.t for f in data]

    pE = [potential_energy(f) for f in data]
    kE = [kinetic_energy(f) for f in data]

    tE = pE .+ kE

    lines!(ax, t, pE, label = "Energía Potencial")
    lines!(ax, t, kE, label = "Energía Cinética")
    lines!(ax, t, tE, label = "Energía Total")
    Legend(fig[:, 2], ax)
    return fig
end
total_energy(f::Frame) = kinetic_energy(f) + potential_energy(f)

function run(;output::String, initial_particles::Vector{Particle}, side_length::Real, duration::Float64)
    step = 0.002
    number_of_real_frames = duration/ step
    time_between_saved_frames = FPS^-1
    global L = Float64(side_length)
    last_percentage = 0.0

    #Creamos el primer Frame con t=0.0
    current_frame = Frame(0.0, initial_particles)

    #Guardamos el primer Frame
    saved_frames = Frame[]
    push!(saved_frames, current_frame)

    #Ejecutamos las simulaciones
    for i in 1:number_of_real_frames
        current_percentage = i/number_of_real_frames *100 
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(Int64,current_percentage), "%\r")
        end
        current_frame = Verlet.stepFrame(current_frame, step, acceleration, update)
        if current_frame.t - last(saved_frames).t > time_between_saved_frames
            push!(saved_frames, current_frame)
        end 
    end

    #Y por último guardamos la simulación
    output_dir = "Tarea1/SimArgon/" * string(output)*  ".out"
    serialize(output_dir, (saved_frames, side_length))
end

function render(; input::String, output::String)
    #Primero cargamos la simulación
    println("Cargando archivo...")
    data = deserialize("Tarea1/SimArgon/" * input)
    side_length = data[2]
    data = data[1]

    println("Archivo cargado")
    @show FPS
    #Creamos la figura y el axis donde dibujar
    fig = Figure()
    ax_main = Axis(fig[1,1], aspect = 1)
    ax_main
    limits!(ax_main, (0, side_length), (0, side_length))

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
    record(fig, "Tarea1/SimArgon/" * output, range(1, number_of_frames); framerate = round(Int64, FPS)) do row_number
        current_frame_i[] = row_number

        current_percentage = (row_number/number_of_frames*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
    end

end

#= Generadores de las distribuciones de partículas que nos interesan=#

function distribucion_uniforme(side_length::Integer)

    #Generamos las partículas en su estado inicial
    initial_particles = Particle[]
    for (i, j) in Iterators.product(1:side_length,1:side_length)
        dist_betw = 1.0
        x = dist_betw * i
        y = dist_betw * j
        push!(initial_particles, Particle([x, y], [0.01, 0.01], 1.0))
    end
    return initial_particles
end