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
