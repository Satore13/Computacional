include("Verlet.jl")
Particle = Verlet.BParticle
Frame = Verlet.BFrame
using Serialization
using CairoMakie, DataStructures
function forceOnPendulum(p::Particle, f::Frame)::Vector{Float64}
    return [-p.mass*9.8*sin(p.x[1])]
end


function wrapAngle(frame)
    new_particles = Particle[]

    for p in frame.particles
        new_angle = p.x[1]
        if new_angle > π
            new_angle = new_angle - 2π
        elseif new_angle < -π
            new_angle = new_angle + 2π
        end    
        push!(new_particles, Particle([new_angle], p.v, p.a, p.mass, p.tags))
    end
    

    return Frame(frame.t, new_particles)
end

function run(θ0, v0)
    FPS = 25
    sec_b_f = FPS^-1
    f = Frame(0.0, [Particle([θ0], [v0], 1.0)])
    frames = Frame[]
    last_frame_t = 0
    for i in 1:1000
        f = Verlet.stepFrame(f, 0.01, forceOnPendulum, wrapAngle)
        if f.t - last_frame_t > sec_b_f
            push!(frames, f)  
            last_frame_t = f.t
        end
    end
    println(".")
    serialize("Tarea1/PenduloReal/"*string(round(θ0, digits = 2)) *"_"* string(v0)*".out", frames)
end 

function runAll()
    run.([0.0, π/4, π/2, π*3/4, 2.967, π, π, π], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 5])
end

function animate(name::String)
    filenames = "Tarea1/PenduloReal/" .* ["0.0_0.0", "0.79_0.0", "1.57_0.0", "2.36_0.0", "2.97_0.0", "3.14_0.0", "3.14_0.1", "3.14_5.0"] .* ".out"
    #Leer los datos y guardarlos en un DataFrame
    println("Cargando archivo...")
    set_of_data = deserialize.(filenames)
    println("Archivo cargado")
    @show typeof(set_of_data)
    current_frame_n = Observable{Int64}(1)

    figure = Figure(resolution = (1600, 400))

    titles = "θ₀ = " .* ["0", "45", "90", "135", "170", "180", "180", "180"] .*"°"

    axises = [Axis(figure[1, i], aspect = 1, title = titles[i]) for i in 1:8]
    axises_ph = [Axis(figure[2, i], aspect = 1) for i in 1:8]
    x = []
    y = []
    p = []
    θ = []
    for frames in set_of_data
        push!(x, [sin(f.particles[1].x[1]) for f in frames])
        push!(y, [-cos(f.particles[1].x[1]) for f in frames])
        push!(p,  [f.particles[1].v[1] for f in frames])
        push!(θ, [f.particles[1].x[1] for f in frames])
    end 

    for ax in axises
        xlims!(ax, (-2.5, 2.5))
        ylims!(ax, (-2.5, 2.5))
    end
    for ax in axises_ph
        xlims!(ax, (-10, 10))
        ylims!(ax, (-10, 10))
    end

    series_points_in_ph = []

    for (index, (ax, axph)) in enumerate(zip(axises, axises_ph))
        xt = @lift(x[index][$current_frame_n])
        yt = @lift(y[index][$current_frame_n])
        pt = @lift(p[index][$current_frame_n])
        θt = @lift(θ[index][$current_frame_n])

        point = @lift(Point2f[($xt, $yt)])
        push!(series_points_in_ph, Observable(Point2f[]))
        leading_point = @lift(Point2f[($θt, $pt)])
        scatter!(ax, point)
        scatter!(axph, series_points_in_ph[index], size = 1)
        scatter!(axph, leading_point, color = :red)
    end

    
    last_percentage = 0.0
    number_of_rows = length(set_of_data[1])[1]
    record(figure, "Tarea1/PenduloReal/"*name*".mp4", range(1, number_of_rows); framerate = 25) do row_number
        current_frame_n[] = row_number

        for (index, ph) in enumerate(series_points_in_ph)
            push!(ph[], Point2f(θ[index][row_number], p[index][row_number]))
        end
        current_percentage = (row_number/number_of_rows*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
    end
end