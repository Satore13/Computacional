include("PenduloDoble.jl")
using CairoMakie

#Aumentamos el tamaño de fuente
set_theme!(Theme(fontsize = 22))

using LaTeXStrings
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
    ax = Axis(fig[1, 1], yscale = log10, xlabel = LaTeXString("Tiempo[s]"), ylabel = L"Error relativo $ξ_{rel}$")

    for sim in sims
        x = getfield.(sim.video, :time)
        y = energia.(sim.video, Ref(sim.parameters))
        y = abs.(y[1] .- y) ./ 30
        y = y[2:end]
        x = x[2:end]
        x = x[y .!= 0]
        y = y[y .!= 0]
        
        label = L"θ_1(0) = θ_2(0) = %$(Int64(round(rad2deg(sim.video[1].y[1]), digits = 0)))°"
        #label = L"$Δt = 10^{%$(Int64(log10(sim.h)))}$ s"
        
        
        lines!(ax, x, y, label = label)
    end
    axislegend(ax, position = :rb)
    fig
end
function plt_error_rel_vs_tiempo_cum(sims::Vector{Simulation})
    fig = Figure()
    ax = Axis(fig[1, 1],xscale = log10, yscale = log10, xlabel = LaTeXString("Tiempo[s]"), ylabel = L"Error relativo cumulativo $ξ_{c}$")

    for sim in sims
        x = getfield.(sim.video, :time)
        y = energia.(sim.video, Ref(sim.parameters))
        y = abs.(y[1] .- y) ./ 30
        for i in eachindex(y[1:(end-1)])
            y[i + 1] = y[i] + y[i + 1]
        end
        y = y[2:end]
        x = x[2:end]
        label = L"θ_1(0) = θ_2(0) = %$(Int64(round(rad2deg(sim.video[1].y[1]), digits = 0)))°"
        #label = L"$Δt = 10^{%$(Int64(log10(sim.h)))}$ s"
        lines!(ax, x, y, label = label)
    end
    axislegend(ax, position = :rb)
    #Obtenemos todos los tiempos y los guardamos en un arreglo
    fig
end
function plt_error_rel_vs_tiempo(sim::Simulation)
    plt_error_rel_vs_tiempo([sim])
end

function plt_poincare(batch::Vector{Simulation}, duration::Union{Integer,Nothing} = nothing)
    fig = Figure()
    ax = Axis(fig[1, 1], 
            xticks = MultiplesTicks(4, π, "π"), 
            yticks = MultiplesTicks(4, π, "π"),
            ylabel = L"$p_2$[kg m²/s]",
            xlabel = L"$θ_2$[rad]",
            aspect = 1.0)

    for sim in batch
        if isnothing(duration)
            svideo = sim.video
        else
            svideo = sim.video[1:duration]
        end
        dots = Point2f[]

        signo_anterior = Int64(sign(svideo[1].y[1]))

        for v in getfield.(svideo, :y)

            if signo_anterior != Int64(sign(v[1]))
                push!(dots, Point2f(v[2], v[4]))
            end
            signo_anterior = Int64(sign(v[1]))
        end
        label = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°,  θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°"
        scatter!(ax, dots, markersize = 4, label = label)
    end
    Legend(fig[1, 2], ax, fontsize = 20)
    fig
end
plt_poincare(batch::Simulation) = plt_poincare([batch])

function plt_proyeccion_posiciones(batch::Vector{Simulation}, duration::Union{Integer,Nothing} = nothing)
    fig = Figure()
    ax = Axis(fig[1, 1], 
            xticks = MultiplesTicks(4, π, "π"), 
            yticks = MultiplesTicks(4, π, "π"),
            ylabel = L"$θ_2$[rad]",
            xlabel = L"$θ_1$[rad]",
            aspect = 1.0)

    for sim in batch
        if isnothing(duration)
            svideo = sim.video
        else
            svideo = sim.video[1:duration]
        end
        dots = Point2f[]

        for frame in svideo
            push!(dots, Point2f(frame.y[1], frame.y[2]))
        end
        newline = "\n"
        
        label = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))° %$(newline)  θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°"
        lines!(ax, dots, markersize = 4, label = label)
    end
    Legend(fig[1, 2], ax, fontsize = 20)
    fig
end

function plot_energia()
    function energia_f(θ1, θ2)
        s = crear_simulacion(θ1 = Float64(θ1*π), θ2 = Float64(θ2*π))
        return energia(s.video[1], s.parameters)
    end

    θ1 = θ2 = Float32.(range(-1, 1, step = 0.01))
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"θ_1", ylabel = L"θ_2", xticks = (-1:0.5:1, [L"-π", L"-\frac{π}{2}", L"0", L"\frac{π}{2}", L"π"]), yticks = (-1:0.5:1, [L"-π", L"-\frac{π}{2}", L"0", L"\frac{π}{2}", L"π"]))
    hm = heatmap!(ax, θ1, θ2, energia_f)
    Colorbar(fig[1, 2], hm, label = LaTeXString("Energía [J]"), ticks = 0:20:60)
    fig
end

function plt_angulos_vs_tiempo(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)

    fig = Figure()
    ax = Axis(fig[1, 1], 
                xlabel = L"$t$[s]",
                ylabel = L"$θ_i$",
                yticks = MultiplesTicks(4, π, "π"),
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    θ1 = Point2f[]
    θ2 = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end
    for frame in svideo
        push!(θ1, Point2f(frame.time, frame.y[1]))
        push!(θ2, Point2f(frame.time, frame.y[2]))

    end
    lines!(ax, θ1, label = L"θ_1")
    lines!(ax, θ2, label = L"θ_2")

    axislegend(ax)

    fig
end

function plt_rel_angulos_vs_tiempo(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)

    fig = Figure()
    ax = Axis(fig[1, 1], 
                yscale = log10,
                xlabel = L"$t$[s]",
                ylabel = L"|\frac{θ_1}{θ_2}|",
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    ratio = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end
    for frame in svideo
        push!(ratio, Point2f(frame.time, abs(frame.y[1] / frame.y[2])))
    end
    popfirst!(ratio)
    lines!(ax, ratio)
    fig
end


function plt_energias_vs_tiempo(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)

    fig = Figure()
    ax = Axis(fig[1, 1], 
                xlabel = L"$t$[s]",
                ylabel = L"$E_i$",
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    θ1 = Point2f[]
    θ2 = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end
    for frame in svideo
        push!(θ1, Point2f(frame.time, energia_1(frame, sim.parameters)))
        push!(θ2, Point2f(frame.time, energia_2(frame, sim.parameters)))

    end
    lines!(ax, θ1, label = L"E_1")
    lines!(ax, θ2, label = L"E_2")

    axislegend(ax)

    fig
end

function plt_rel_energias_vs_tiempo(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)

    fig = Figure()
    ax = Axis(fig[1, 1], 
                xlabel = L"$t$[s]",
                ylabel = L"$\frac{E_1}{E_2}$",
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    ratio = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end
    for frame in svideo
        push!(ratio, Point2f(frame.time, energia_1(frame, sim.parameters) / energia_2(frame, sim.parameters)))

    end
    lines!(ax, ratio)

    fig
end


function plt_peaks_angulos_vs_tiempo(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)

    fig = Figure()
    ax = Axis(fig[1, 1], 
                xlabel = L"$t$[s]",
                ylabel = L"$θ_i$",
                yticks = MultiplesTicks(4, π, "π"),
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    θ1a = Point2f[]
    θ2a = Point2f[]
    θ1na = Point2f[]
    θ2na = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end

    previous_θ1 = 0.0
    previous_θ2 = 0.0

    θ1_growing = true
    θ2_growing = true

    for frame in svideo
        θ1 = frame.y[1]
        θ2 = frame.y[2]
        if θ1_growing && previous_θ1 > θ1 
            θ1_growing = false
            push!(θ1a, Point2f(frame.time, θ1))
        elseif !θ1_growing && previous_θ1 < θ1
            θ1_growing = true
            push!(θ1a, Point2f(frame.time, θ1))
        end
        if θ2_growing && previous_θ2 > θ2
            θ2_growing = false
            push!(θ2a, Point2f(frame.time, θ2))
        elseif !θ2_growing && previous_θ2 < θ2
            θ2_growing = true
            push!(θ2a, Point2f(frame.time, θ2))
        end
        previous_θ1 = θ1
        previous_θ2 = θ2
        push!(θ1na, Point2f(frame.time, frame.y[1]))
        push!(θ2na, Point2f(frame.time, frame.y[2]))
    end
    #lines!(ax, θ1na, label = L"θ_1")
    #lines!(ax, θ2na, label = L"θ_2")
    scatter!(ax, θ1a, label = L"$θ_1$ picos")
    scatter!(ax, θ2a, label = L"$θ_2$ picos")
    

    axislegend(ax)

    fig
end

function plt_ph_ang1(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)
    fig = Figure()
    ax = Axis(fig[1, 1], 
                xlabel = L"$θ_1$[rad]",
                ylabel = L"$\frac{dθ_1}{dt}$[kg m² /s]",
                xticks = MultiplesTicks(4, π, "π"),
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    pos = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end
    for frame in svideo
        push!(pos, Point2f(frame.y[1], angles_velocity(frame, sim.parameters)[1]))
    end
    lines!(ax, pos)
    fig
end

function plt_ph_ang2(sim::Simulation, rang::Union{UnitRange, Nothing} = nothing)
    fig = Figure()
    ax = Axis(fig[1, 1], 
                xlabel = L"$θ_2$[rad]",
                ylabel = L"$\frac{dθ_2}{dt}$[rad/s]",
                xticks = MultiplesTicks(4, π, "π"),
                yticks = MultiplesTicks(4, π, "π"),
                title = L"θ_1(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[1])))°, θ_2(0) = %$(round(digits = 2,rad2deg(sim.video[1].y[2])))°")

    pos = Point2f[]
    if isnothing(rang)
        svideo = sim.video
    else 
        svideo = sim.video[rang]
    end
    for frame in svideo
        push!(pos, Point2f(frame.y[2], angles_velocity(frame, sim.parameters)[2]))
    end
    lines!(ax, pos)
    fig
end
