include("PenduloDoble.jl")
using GLMakie
#En este archivo irán exclusivamente las funciones de ploteo y/o ajustes
function animar_simulacion(sim::Simulation, filename::String = "out.mp4")
    filename = "Tarea4/videosPD/" * filename



    r1 = Observable{Point2f}(r1r2_de_polares(sim.video[1], sim.parameters)[1])
    r2 = Observable{Point2f}(r1r2_de_polares(sim.video[1], sim.parameters)[2])
    cuerdas = @lift [Point2f(0.0, 0.0), $r1, $r2]

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1)
    scatter!(ax, r1, color = :black)
    scatter!(ax, r2, color = :black)
    lines!(ax, cuerdas, color = :black)

    side = sim.parameters[:l1] + sim.parameters[:l2]
    xlims!(ax, (-side * 1.2, side * 1.2))
    ylims!(ax, (-side * 1.2, side * 1.2))

    @show sim.fps
    @show sim.video[end].time
    record(fig, filename, eachindex(sim.video), framerate = Int64(sim.fps)) do i
        current_frame = sim.video[i]
        @show current_frame.time
        r1[], r2[] = r1r2_de_polares(current_frame, sim.parameters)
    end
end