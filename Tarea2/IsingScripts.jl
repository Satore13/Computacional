include("ModeloIsing.jl")
using GLMakie
using Serialization
using LsqFit
using Measurements
function plotear_EyMvspMC(datos::Vector{Red}, guardadas_cada::Integer = 1, filename::Union{Nothing, String} = nothing)
    N = size(datos[1])
    T = round(digits = 5, datos[1].Temperatura)
    fig = Figure()
    ax_E = Axis(fig[1, 1], xlabel = "pMC")
    ax_M = Axis(fig[2, 1], xlabel = "pMC")
    x = collect(1:size(datos, 1))
    x = (x .- 1 ) .*guardadas_cada
    yE = [energia_norm(red) for red in datos]
    yM = [abs(magnetizacion_norm(red)) for red in datos]

    lE = lines!(ax_E, x, yE, color = :red)
    lM = lines!(ax_M, x, yM, color = :black)

    ylims!(ax_M, (0 , 1))
    ylims!(ax_E, (-2, 0))

    Legend(fig[:, 2], [lE, lM], ["Energía", "Magnetización"])
    Label(fig[0, :], "Red $(N)x$(N)\t T = $(T)", textsize = 25)
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

function plotear_configuracion_anim(datos::Vector{Red})
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
    scatter!(ax, x,  yM, label = label, color = :black)
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


function plotear_correlacion(red::Red)
    max_val = ceil(Int64, size(red)/2)

    r = collect(1:max_val)
    G = [global_correlation(red, rᵢ) for rᵢ in r]

    scatter(r, G)

end

function ajuste_MvsT!(ax, filename::String)
    data = deserialize(filename)
    x = Float64[]
    y = Float64[]
    for d in data
        if d[1] < 2.26 && d[3] > 0.8
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
    filenames = "Tarea2/output/2_5to3_5_segundo/n" .* ["16", "32", "64"#=, "128"=#] .* "EyMvsT.out"
    titles = ["16x16", "32x32", "64x64"]#, "128x128"]
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
        if (d[1] < 2.26) && (d[3] < 0.6)
            continue
        end
        push!(new_data, d)
    end
    serialize(filename, new_data)
end