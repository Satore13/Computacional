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
    h::Float64
    fps::Float64
    time_scale::Float64
    function Simulation(header::Vector{Symbol}, in_cond::Vector{Float64}, dyn::Function, h::Float64, fps::Float64; time_scale = 1.0)
        if length(header) != length(in_cond)
            error("El cabecero y las condiciones iniciales han de tener la misma longitud")
        end
        try 
            y = dyn(in_cond, 0.0)
            length(y) != length(in_cond) && error("La función no devuelve un vector válido")

        catch e
            println("Error al comprobar la función: ", e)
        end
        initial_frame = Frame(0.0, in_cond)
        
        return new(header, Frame[initial_frame], dyn, h, fps, time_scale)
    end
end

function pushFrame!(s::Simulation, f::Frame)
    length(f.y) != length(s.header) && error("Frame no válido en esta simulación")
    push!(s.video, f)
end

function stepFrame(f::Frame, foo::Function, h::Float64)
    y = f.y
    t = f.time
    k = Vector{Vector{Float64}}(undef, 4)
    k[1] = h .* foo(y, t)
    k[2] = h .* foo(y .+ 0.5 .* k[1], t + 0.5 * h)
    k[3] = h .* foo(y .+ 0.5 .* k[2], t + 0.5 * h)
    k[4] = h .* foo(y .+ k[3], t + h)

    return Frame(t+ h, y .+ (1.0/6.0) .* (k[1] .+ 2 .*k[2] .+ 2 .*k[3] .+ k[4]))
end

function loop!(s::Simulation, length::Float64)
    length = length * s.time_scale
    initial_time = last(s.video).time
    final_time = initial_time + length
    @show final_time
    current_time = initial_time
    last_saved_frame = initial_time
    time_between_frames = s.fps^-1 * s.time_scale

    current_frame = last(s.video)
    last_percentage = 0.0

    #Ejecutamos el bucle mientas que estemos por detrás del fin del tiempo que nos piden simular
    while current_time < final_time
        #Simulamos un frame con el paso que venga determinado por la simulación
        current_frame = stepFrame(current_frame, s.dyn, s.h)
        #Guardamos el tiempo del frame más nuevo
        current_time = current_frame.time
        #Sólo guardamos un frame cuando vemos que hay una distancia superior a la que estipulamos que hemos de tener para tener los fps pedidos
        if current_time - last_saved_frame > time_between_frames
            pushFrame!(s, current_frame)
            last_saved_frame = current_time
        end

        #Algoritmo para imprimir el porcentaje de progreso
        current_percentage = (current_time - initial_time) / length * 100.0
        if current_percentage - last_percentage > 0.5
            print("Progreso: $(round(current_percentage, digits = 2))%\r")
        end

    end
    return nothing
end