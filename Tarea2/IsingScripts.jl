include("ModeloIsing.jl")
using GLMakie
using Serialization
using LsqFit
using Measurements
function plotear_EyMvspMC(datos::Vector{Vector{Red}}, colors::Vector{Symbol}; gc::Integer = 1, filename::Union{Nothing, String} = nothing)
    if size(datos, 1) ≠ size(colors, 1)
        error("El array de datos y de colores ha de ser del mismo tamaño")
    end

    fig = Figure()
    ax_E = Axis(fig[1, 1], xlabel = "pMC", title = "Energía")
    ax_M = Axis(fig[2, 1], xlabel = "pMC", title = "Magnetización")
    plts_energia = []
    plts_magnetizacion = []
    labels = String[]
    for (d, c) in zip(datos, colors)
        N = size(d[1])
        T = round(digits = 3, d[1].Temperatura)
        x = collect(1:size(d, 1))
        x = (x .- 1 ) .*gc
        yE = [energia_norm(red) for red in d]
        yM = [abs(magnetizacion_norm(red)) for red in d]
        push!(plts_energia, lines!(ax_E, x, yE, color = c))
        push!(plts_magnetizacion, lines!(ax_M, x, yM, color = c))
        push!(labels, "T = $T, N = $N")
    end
    ylims!(ax_M, (0 , 1))
    ylims!(ax_E, (-2, 0))

    Legend(fig[:, 2], [[plts_energia[i], plts_magnetizacion[i]] for i in eachindex(plts_energia)], labels, merge = true, unique = true)
    return fig
end

function plotear_configuracion(red::Red)
    N = size(red)
    datos_s = getfield.(red.Nudos, :val)
    x = collect(1:N)
    y = collect(1:N)


    fig = Figure()
    ax = Axis(fig[1, 1])
    z = [datos_s[i, j] for i in x, j in y]
    heatmap!(ax, z, colormap = Reverse(:greys))

    return fig
end

function plotear_configuracion_anim(datos::Vector{Red}, filename::String = "Tarea2/sim.mp4")
    N = size(datos[1])
    datos_s = [getfield.(red.Nudos, :val) for red in datos]
    x = collect(1:N)
    y = collect(1:N)


    fig = Figure()
    ax = Axis(fig[1, 1])
    index = Observable{Int64}(1)

    z = @lift [datos_s[$index][i, j] for i in x, j in y]
    @lift heatmap!(ax, $z, colormap = Reverse(:greys))
    @lift print(string($(index)) * "\r")
    record(fig, filename, 1:size(datos_s, 1),framerate = 10) do i
        index[] = i
    end
end

function average_E(data::Vector{Red})
    return sum(energia_norm.(data)) / size(data, 1)
end

function average_M(data::Vector{Red})
    return sum(magnetizacion_norm.(data)) / size(data, 1)
end

function plotear_EvsT!(ax::Axis, label, data::Vector{Tuple{Float64, Float64, Float64}})
    x = getindex.(data, 1)
    yE = getindex.(data, 2)
    scatter!(ax, x,  yE, label = label)
end

function plotear_MvsT!(ax::Axis, label, data::Vector{Tuple{Float64, Float64, Float64}})
    x = getindex.(data, 1)
    yM = getindex.(data, 3)
    scatter!(ax, x,  yM, label = label, color = :black, markersize = 5)
end

function procesar_EyMvsT(n)
    T = collect(range(0.1, 5, step = 0.1))
    output = Tuple{Float64, Float64, Float64}[]
    gc = 10
    for t in T
        t = round(digits = 5, t)
        
        filename = "Tarea2/output/n$n"*"T$t"*"gc$gc"*".out"
        println("Cargando archivo: $filename")
        data = deserialize(filename)[(end - 200):end]
        push!( output, (t, average_E(data), abs(average_M(data))) )
    end
    return output
end

function plotear_EvsTvsn()
    N = [16, 32, 64, 128]
    fig = Figure()
    ax = Axis(fig[1, 1])

    ylims!(ax, (-2.2, -0.2))
    xlims!(ax, (0, 5))

    for n in N
        plotear_EvsT!(ax, "N = $n" , deserialize("Tarea2/output/$(n)E_MvsT.out"))
    end
    axislegend(ax, position = :rb)
    return fig
end
function plotear_MvsTvsn()
    N = [128]
    fig = Figure()
    ax = Axis(fig[1, 1])
    ylims!(ax, (-0.2, 1.2))
    xlims!(ax, (1, 4))
    for n in N
        plotear_MvsT!(ax, "N = $n" , deserialize("Tarea2/output/n$(n)EyMvsT.out"))
    end
    axislegend(ax, position = :rt)
    return fig
end


function plotear_correlacion(datos::Vector{Vector{NTuple{2, Float64}}}, labels::Vector{String})
    fig = Figure()
    ax = Axis(fig[1, 1])
    for (d, l) in zip(datos, labels)
        scatter!(ax, d, label = l)
        lines!(ax, d, label = l)
        ylims!(ax, (0,1.1))
    end
    Legend(fig[:, 2], ax, merge = true, unique = true)
    fig
end

function plotear_correlacion(datos::Vector{NTuple{2, Float64}}, label::String)
    plotear_correlacion([datos], [label])
end

function plotear_correlaciones(N::Int64)
    labels = String[]
    datos = Vector{NTuple{2, Float64}}[]
    temperaturas = [1.1:0.2:2.7 ; Dict(16 => 2.313, 32 => 2.290, 64 => 2.280, 128 => 2.269)[N]]
    N = string(N)
    for t in temperaturas
        t = string(round(t, digits = 3))
        filename = "Tarea2/output/f_corr/f_corr_n$(N)T$(t).out"
        push!(labels, "T = $(t)")
        push!(datos, deserialize(filename))
    end
    plotear_correlacion(datos, labels)
end

function ajuste_funcion_correlacion(datos::Vector{NTuple{2, Float64}})
    @. model(x, p) = exp(-x / p[1]) + p[2]

    x = getindex.(datos, 1)
    y = getindex.(datos, 2)
    p0 = [1.0, 0.0]

    fit = curve_fit(model, x, y, p0)
    return (coef(fit)[1], stderror(fit)[1])
end

function plotear_l_corr_vs_temp()
    fig = Figure()
    ax = Axis(fig[1, 1])
    for (n, Tc, c, l) in zip([32, 64, 128], [2.290, 2.280, 2.269], [:blue, :orange, :red], ["32x32", "64x64", "128x128"])
        ξvsT = NTuple{3, Float64}[]
        for t in sort([1.3:0.2:3.5; Tc])
            filename = "Tarea2/output/f_corr/f_corr_n$(n)T$(t).out"
            datos = deserialize(filename)
            output = ajuste_funcion_correlacion(datos)
            push!(ξvsT, (t, output...))
        end
        x = getindex.(ξvsT, 1)
        y = getindex.(ξvsT, 2)
        dy = getindex.(ξvsT, 3)
        errorbars!(ax, x, y, dy)
        #scatter!(ax, x, y)
        lines!(ax, x, y, color = c, label = l)
    end
    axislegend(ax)
    return fig
end

function ajuste_MvsT!(ax, filename::String)
    data = deserialize(filename)
    x = Float64[]
    y = Float64[]
    for d in data
        if d[3] > 0.7
            push!(x, d[1])
            push!(y, d[3])
        end
    end
    @.  model(x, p) =  (p[1] - p[2]/(sinh(2/x)^4))^p[3]
    
    p0 = [1.0, 0.01, 1.0]
    fit = curve_fit(model, x, y, p0)
    #fit = LsqFit.LsqFitResult(fit.)
    p = coef(fit) .± stderror(fit)
    @show p
    @show Tc = 2/(asinh(  (p[2]/p[1]) ^ (1/4)) )   

    plotear_MvsT!(ax, "", data)

    x = collect(range(x[1], Measurements.value(Tc - 0.001), length = 100))
    y = model(x, coef(fit))
    lines!(ax, x, y, color = :red)
end

function ajuste_MvsT(filename::String)
    ajuste_MvsT!(Axis(Figure()[1, 1]), filename)
    current_figure()
end

function plotear_ajustes()
    filenames = "Tarea2/output/2_5to3_5_segundo/n" .* ["16", "32", "64", "128"] .* "EyMvsT.out"
    titles = ["16x16", "32x32", "64x64", "128x128"]
    fig = Figure()

    axes = [Axis(fig[i, j]) for i in 1:2, j in 1:2]

    @show filenames
    @show axes
    
    for (ax, file, title) in zip(axes,  filenames, titles)
        ajuste_MvsT!(ax, file)
        ax.title = title
    end
    return fig
end

function purgar_datos(filename::String)
    data = deserialize(filename)
    @assert typeof(data) == Vector{Tuple{Float64, Float64, Float64}}
    new_data = empty(data)
    for d in data
        if (d[1] < 2.26) && (d[3] < 0.65)
            continue
        end
        push!(new_data, d)
    end
    serialize(filename, new_data)
end

function plotear_calores_especificos()
    fig = Figure()
    ax = Axis(fig[1, 1])
    Ns = [16, 32, 64, 128]
    xlims!(ax, (2.0, 3.0))
    ylims!(ax, (0.0, 5.0))
    for (n, c, tc) in zip(Ns, [:blue, :red, :green, :black], [2.313, 2.290, 2.280, 2.269])
        filename = "Tarea2/output/calor_especifico/ct_n$(n)2to3.out"
        data = deserialize(filename)
        x = getindex.(data, 1)
        y = getindex.(data, 2)
        scatter!(ax, x, y, color = c, label = "$(n)x$(n)")
        lines!(ax, x, y, color = c, label = "$(n)x$(n)")
        vlines!(ax, tc, color = c, label = "$(n)x$(n)")
    end
    axislegend(ax, merge = true)
    return fig
end

function plotear_tc_vs_L()
    T_c = [2.313, 2.290, 2.280, 2.269] 
    wt = [0.006, 0.002, 0.003, 0.002] .^ -2
    L = [16, 32, 64, 128] 
    
    @. model(x, p) = 1/p[1] - p[2] * x 
 
    p0 = [2.0, 1.0]

    x = L .^ -1
    y = T_c .^ -1

    @show x
    @show y
    @show wt
    fit = curve_fit(model, x, y,wt,   p0)

    @show coef(fit) .± stderror(fit)

    ax = Axis(Figure()[1, 1], ylabel = L"T_c^{-1}", xlabel = L"N^{-1}")

    scatter!(ax, x, y)
    
    x = collect(range(x[begin], x[end], length = 10))
    y = model(x, coef(fit))
    @show x 
    @show y
    lines!(ax, x, y)
    current_figure()
end