include("Verlet.jl")
Particle = Verlet.Particle
Frame = Verlet.Frame
@enum Planetas begin
    Tierra=1
    Mercurio
    Venus
    Marte
    Jupiter
    Saturno
    Urano
    Neptuno
end
import Base.string
string(p::Planetas) = string(Int(p))
using CSV
function norm(x::Vector{Float64})::Float64
    normsq = 0.0
    for xi in x
        normsq += xi^2
    end
    return √(normsq)
end
function gravitationalAcc(p::Particle, f::Frame)
    m_sol = 1
    return (1 / norm(p.x)^3) .* (-p.x)
end
@time begin
    function run()
        FPS = 24
        
        f = Frame(0.0, [Particle([1.0, 0.0], [0.0, 1.0], 3.002513826043238e-6)])
        simulation_data = Verlet.dataFrameRowFromFrame(f)

        time_between_frames = FPS^-1
        time_of_last_frame = 0.0

        length_of_sim = 5
        step = 0.001
        number_of_frames = length_of_sim/step

        show(number_of_frames)
        for i in 1:number_of_frames
            f = Verlet.stepFrame(f, step, gravitationalAcc)
            if f.t - time_of_last_frame > time_between_frames
                Verlet.addFrameRowToDataFrame!(simulation_data, f)
                time_of_last_frame = f.t
            end
        end
        CSV.write("datos_sistema_solar.out", simulation_data)
    end
    run()
end


