include("Verlet.jl")
Particle = Verlet.Particle
Frame = Verlet.Frame
using CSV
using DataFrames

function forceOnPendulum(p::Particle, f::Frame)::Vector{Float64}
    return [-p.mass*9.8*sin(p.x[1])]
end

@time begin
    N = 1
    FPS = 30
    seconds_between_frames = FPS^-1
    global f = Frame(0.0, [Particle([π*8/9], [0.0], 1.0)])
    simulation_data = Verlet.dataFrameRowFromFrame(f)
    global seconds_of_last_frame = 0.0
    for i in 1:10000
        global f = Verlet.stepFrame(f, 0.001, forceOnPendulum)
        if f.t - seconds_of_last_frame > seconds_between_frames
            Verlet.addFrameRowToDataFrame!(simulation_data, f)
            global seconds_of_last_frame = f.t
        end
    end
    println(simulation_data)
    CSV.write("datos_simulacion.out",simulation_data)


end 