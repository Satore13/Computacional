struct Frame
   time::Float64
   y::Vector{Float64}
end

#=Todas las unidades temporales se asumen a priori en segundos, time_scale sirve sólo para el guardado de frames. Los únicos datos que se dan en la unidad 
    ya escalada son length en la función loop. 
=#
struct Simulation
    header::Vector{Symbol}
    video::Vector{Frame}
    dyn::Function
    fps::Float64
    time_scale::Float64
    parameters::Dict{Symbol, Float64}
    function Simulation(header::Vector{Symbol}, in_cond::Vector{Float64}, dyn::Function, fps::Float64; time_scale = 1.0, parameters = Dict{Symbol, Float64})
        if length(header) != length(in_cond)
            error("El cabecero y las condiciones iniciales han de tener la misma longitud")
        end
        y = dyn(in_cond, 0.0, parameters)
        length(y) != length(in_cond) && error("La función no devuelve un vector válido")

        initial_frame = Frame(0.0, in_cond)
        
        return new(header, Frame[initial_frame], dyn, fps, time_scale, parameters)
    end
end

function pushFrame!(s::Simulation, f::Frame)
    length(f.y) != length(s.header) && error("Frame no válido en esta simulación")
    push!(s.video, f)
end

function stepFrame(f::Frame, foo::Function, h::Float64, parameters)
    y = f.y
    t = f.time
    k = Vector{Vector{Float64}}(undef, 4)

    k[1] = h .* foo(y, t, parameters)
    k[2] = h .* foo(y .+ 0.5 .* k[1], t + 0.5 * h, parameters)
    k[3] = h .* foo(y .+ 0.5 .* k[2], t + 0.5 * h, parameters)
    k[4] = h .* foo(y .+ k[3], t + h, parameters)

    return Frame(t+ h, y .+ (1.0/6.0) .* (k[1] .+ 2 .*k[2] .+ 2 .*k[3] .+ k[4]))
end

function loop!(sim::Simulation, length::Float64, h0::Float64)
    length = length * sim.time_scale
    initial_time = last(sim.video).time
    final_time = initial_time + length
    @show final_time
    current_time = initial_time
    last_saved_frame = initial_time
    time_between_frames = sim.fps^-1 * sim.time_scale
    
    current_frame = last(sim.video)
    last_percentage = 0.0

    h = h0
    #El error que toleraremos será ε = h^5
    ε_tolerado = h0^5
    #Ejecutamos el bucle mientas que estemos por detrás del fin del tiempo que nos piden simular
    while current_time < final_time
        #Este bucle se ejecutará hasta que tengamos un paso apropiado
        while true
            #Primero tenemos que simular un paso con el h dado
            current_frame_with_current_h = stepFrame(current_frame, sim.dyn, h, sim.parameters)
            #y otro con el h medios
            current_frame_with_halved_h = stepFrame(current_frame, sim.dyn, h / 2.0, sim.parameters)

            #Obtenemos el mayor error entre todas las coordenadas de y calculado con h y de y calculado con h medios
            ε = maximum([ 16.0/15.0 * abs(y_h - y_h_halved) for (y_h, y_h_halved) in zip( current_frame_with_current_h.y, current_frame_with_halved_h.y)])

            s = max((ε / ε_tolerado)^0.2, 1e-8)
            h_max = h/s

            #Si h < h_max hemos calculado todo con más precisión de la necesaria con lo que duplicamos su valor
            if h < h_max
                h = 2 * h
                current_frame = current_frame_with_current_h
                break
            end
            if s < 2
                current_frame = current_frame_with_halved_h
                break
            end
            #Si s > 2 ninguno de los dos cálculos es lo suficientemente preciso así que hay que recalcular todo con h = h_max
            h = h_max
        end
        #Guardamos el tiempo del frame más nuevo
        current_time = current_frame.time
        #Sólo guardamos un frame cuando vemos que hay una distancia superior a la que estipulamos que hemos de tener para tener los fps pedidos
        if current_time - last_saved_frame > time_between_frames
            pushFrame!(sim, current_frame)
            last_saved_frame = current_time
        end

        #Algoritmo para imprimir el porcentaje de progreso
        current_percentage = (current_time - initial_time) / length * 100.0
        if current_percentage - last_percentage > 0.1
            print("Progreso: $(round(current_percentage, digits = 2))%, paso = $h\r")
        end

    end
    return nothing
end