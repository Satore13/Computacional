include("ModeloIsing.jl")
using GLMakie
using Serialization
using SciPy
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
    scatter!(ax, x,  yM, label = label)
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


function procesar_EyMvsTzoom(n)
    T = collect(range(1.9, 3.0, length = 40))
    T = round.(T, digits = 5)
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
    N = [16, 32, 64, 128]
    fig = Figure()
    ax = Axis(fig[1, 1])
    ylims!(ax, (-0.2, 1.2))
    xlims!(ax, (0, 5))
    for n in N
        plotear_MvsT!(ax, "N = $n" , deserialize("Tarea2/output/$(n)E_MvsT.out"))
    end
    axislegend(ax, position = :rt)
    return fig
end

function plotear_EvsTvsnzoom()
    N = [16, 32, 64, 128]
    fig = Figure()
    ax = Axis(fig[1, 1])

    ylims!(ax, (-2.2, -0.2))
    xlims!(ax, (1.8, 3.1))

    for n in N
        plotear_EvsT!(ax, "N = $n" , deserialize("Tarea2/output/$(n)E_MvsTzoom.out"))
    end
    axislegend(ax, position = :rb)
    return fig
end
function plotear_MvsTvsnzoom()
    N = [16, 32, 64, 128]
    fig = Figure()
    ax = Axis(fig[1, 1])
    ylims!(ax, (-0.2, 1.2))
    xlims!(ax, (1.8, 3.1))
    for n in N
        plotear_MvsT!(ax, "N = $n" , deserialize("Tarea2/output/$(n)E_MvsTzoom.out"))
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

function ajuste_MvsT(filename::String)
    data = deserialize(filename)
    x = getindex.(data[1:20], 1)
    y = getindex.(data[1:20], 3)
    @show x
    @show y
    function sigmoid(β, x)
        1/(1 .+ exp.((x .- β[1]) .* 50))
    end

    odr_data = SciPy.odr.RealData(x, y)
    odr_model = SciPy.odr.Model(sigmoid)
    odr_job = SciPy.odr.ODR(odr_data, odr_model, [2.65])

    odr_output = odr_job.run()
    
    return (odr_output.beta, odr_output.sd_beta)
end

#=
function ajuste(x, y, xe, ye; model = line, seed)
    odr_data = SciPy.odr.RealData(x, y, xe, ye)
    odr_model = SciPy.odr.Model(model)
    odr_job = SciPy.odr.ODR(odr_data, odr_model, seed)
    odr_output = odr_job.run()
    ajuste = DataFrame()
    ajuste.β = odr_output.beta
    ajuste.βe = odr_output.sd_beta
    ajuste.β_ = ajuste.β .± ajuste.βe
    return ajuste
end
=#