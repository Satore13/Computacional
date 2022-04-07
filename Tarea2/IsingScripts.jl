include("ModeloIsing.jl")
using CairoMakie

function plotear_EyMvspMC(datos::Vector{Red})
    fig = Figure()
    ax_E = Axis(fig[1, 1], xlabel = "pMC")
    ax_M = Axis(fig[2, 1], xlabel = "pMC")
    x = collect(1:size(datos, 1))
    yE = [energia(red) for red in datos]
    yM = [magnetizacion(red) for red in datos]

    lE = lines!(ax_E, x, yE, color = :red)
    lM = lines!(ax_M, x, yM, color = :black)
    Legend(fig[:, 2], [lE, lM], ["Energía", "Magnetización"])
    return fig
end

function plotear_configuracion(datos::Vector{Red})
    N = size(datos[1])
    datos_s = [getfield.(red.Nudos, :val) for red in datos]

    x = collect(1:N)
    y = collect(1:N)

    z = [datos_s[1][i, j] for i in x, j in y]

    heatmap(z)

end