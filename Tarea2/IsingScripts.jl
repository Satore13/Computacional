include("ModeloIsing.jl")
using GLMakie
using Serialization
function plotear_EyMvspMC(datos::Vector{Red}, guardadas_cada::Integer = 1)
    N = size(data[1])
    fig = Figure()
    ax_E = Axis(fig[1, 1], xlabel = "pMC")
    ax_M = Axis(fig[2, 1], xlabel = "pMC")
    x = collect(1:size(datos, 1))
    x = (x .- 1 ) .*guardadas_cada
    yE = [energia(red)/N^2 for red in datos]
    yM = [abs(magnetizacion(red))/N^2 for red in datos]

    lE = lines!(ax_E, x, yE, color = :red)
    lM = lines!(ax_M, x, yM, color = :black)

    ylims!(ax_M, (0 , 1))
    ylims!(ax_E, (-2, 0))

    Legend(fig[:, 2], [lE, lM], ["Energía", "Magnetización"])
    Label(fig[0, :], "Red $N"* "x"* "$N", textsize = 25)
    return fig
end



function plotear_configuracion(datos::Vector{Red})
    N = size(datos[1])
    datos_s = [getfield.(red.Nudos, :val) for red in datos]
    x = collect(1:N)
    y = collect(1:N)


    fig = Figure()
    ax = Axis(fig[1, 1])
    index = Observable{Int64}(1)

    z = @lift [datos_s[$index][i, j] for i in x, j in y]
    @lift heatmap!(ax, $z, colormap = Reverse(:greys))

    record(fig, "Tarea2/sim.mp4", 1:size(datos_s, 1),framerate = 10) do i
        index[] = i
    end

end