include("ModeloIsing.jl")
using GLMakie
using Serialization
function plotear_EyMvspMC(datos::Vector{Red}, guardadas_cada::Integer = 1)
    N = size(data[1])
    fig = Figure()
    ax_E = Axis(fig[1, 1], xlabel = "pMC")
    ax_M = Axis(fig[2, 1], xlabel = "pMC")
    x = collect(1:size(datos, 1))
    x[2:end] = (x[2:end] .- 1)*guardadas_cada
    yE = [energia(red)/N^2 for red in datos]
    yM = [abs(magnetizacion(red))/N^2 for red in datos]

    lE = lines!(ax_E, x, yE, color = :red)
    lM = lines!(ax_M, x, yM, color = :black)

    hlines!(ax_E, -2.0)
    hlines!(ax_M, 1.0)
    Legend(fig[:, 2], [lE, lM], ["Energía", "Magnetización"])
    Label(fig[0, :], "Red $N"* "x"* "$N", textsize = 25)
    return fig
end



function plotear_configuracion(datos::Vector{Red})
    N = size(datos[1])
    datos_s = [getfield.(red.Nudos, :val) for red in datos]

    x = collect(1:N)
    y = collect(1:N)

    z = [datos_s[1][i, j] for i in x, j in y]
    heatmap(z, colormap = Reverse(:greys))
end