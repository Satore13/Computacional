include("PenduloDoble.jl")
using CairoMakie
using Printf
#En este archivo irán exclusivamente las funciones de ploteo y/o ajustes
function animar_simulacion(sim::Simulation ,filename::String = "out.mp4"; length = nothing)
    filename = "Tarea4/videosPD/" * filename

    video = isnothing(length) ? sim.video : sim.video[1:length]

    r1 = Observable{Point2f}(r1r2_de_polares(video[1], sim.parameters)[1])
    r2 = Observable{Point2f}(r1r2_de_polares(video[1], sim.parameters)[2])
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
    @show video[end].time
    record(fig, filename, eachindex(video), framerate = Int64(sim.fps)) do i
        current_frame = video[i]
        print(" Tiempo del frame actual: $(round(current_frame.time, digits = 2))\r")
        r1[], r2[] = r1r2_de_polares(current_frame, sim.parameters)
    end
    println()
end


function plt_error_rel_vs_tiempo(sims::Vector{Simulation})
    fig = Figure()
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = "Tiempo[s]", ylabel = L"Error relativo $\xi_{rel}$")

    for sim in sims
        x = getfield.(sim.video, :time)
        y = energia.(sim.video, Ref(sim.parameters))
        y = abs.(y[1] .- y) ./ 30
        y = y[2:end]
        x = x[2:end]
        x = x[y .!= 0]
        y = y[y .!= 0]
        
        label = L"θ_1(0) = θ_2(0) = %$(Int64(round(rad2deg(sim.video[1].y[1]), digits = 0)))"
        #label = L"$Δt = 10^{%$(Int64(log10(sim.h)))}$ s"
        
        
        lines!(ax, x, y, label = label)
    end
    axislegend(ax, position = :lt)
    #Obtenemos todos los tiempos y los guardamos en un arreglo
    fig
end
function plt_error_rel_vs_tiempo_cum(sims::Vector{Simulation})
    fig = Figure()
    ax = Axis(fig[1, 1],xscale = log10, yscale = log10, xlabel = "Tiempo[s]", ylabel = L"Error relativo cumulativo $\xi_{c}$")

    for sim in sims
        x = getfield.(sim.video, :time)
        y = energia.(sim.video, Ref(sim.parameters))
        y = abs.(y[1] .- y) ./ 30
        for i in eachindex(y[1:(end-1)])
            y[i + 1] = y[i] + y[i + 1]
        end
        y = y[2:end]
        x = x[2:end]
        label = L"θ_1(0) = θ_2(0) = %$(Int64(round(rad2deg(sim.video[1].y[1]), digits = 0)))"
        #label = L"$Δt = 10^{%$(Int64(log10(sim.h)))}$ s"
        lines!(ax, x, y, label = label)
    end
    axislegend(ax, position = :lt)
    #Obtenemos todos los tiempos y los guardamos en un arreglo
    fig
end
function plt_error_rel_vs_tiempo(sim::Simulation)
    plt_error_rel_vs_tiempo([sim])
end

function plt_poincare(batch::Vector{Simulation})
    fig = Figure()
    ax = Axis(fig[1, 1])

    for sim in batch
        dots = Point2f[]

        signo_anterior = Int64(sign(sim.video[1].y[1]))

        for v in getfield.(sim.video, :y)

            if signo_anterior != Int64(sign(v[1]))
                push!(dots, Point2f(v[2], v[4]))
            end
            signo_anterior = Int64(sign(v[1]))
        end
        label = L"θ_1 = θ_2 = %$(rad2deg(sim.video[1].y[1]))"
        scatter!(ax, dots, markersize = 3, label = label)
    end
    axislegend(ax)
    fig
end