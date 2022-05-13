struct Frame
   time::Float64
   y::Vector{Float64}
end

struct Simulation
    header::Vector{Symbol}
    video::Vector{Frame}
    dyn::Function
    h::Float64
    function Simulation(header::Vector{Symbol}, in_cond::Vector{Float64}, dyn::Function, h::Float64)
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
        
        return new(header, Frame[initial_frame], dyn, h)
    end
end

function pushFrame!(s::Simulation, f::Frame)
    length(f.y) != length(header) && error("Frame no válido en esta simulación")
    push!(s.video, f)
end

function stepFrame!(f::Frame, foo::Function)
    y = f.y
    h = f.h
    t = f.time
    k = Vector{Vector{Float64}}(undef, 4)
    k[1] = h .* dyn(y, t)
    k[2] = h .* dyn(y .+ 0.5 .* k[1], t + 0.5 * h)
    k[3] = h .* dyn(y .+ 0.5 .* k[2], t + 0.5 * h)
    k[4] = h .* dyn(y .+ k[3], t + h)

    for i in eachindex(y)
        f.y[i] = y[i] + (1.0/6.0) * (k[1][i] + 2*k[2][i] + 2*k[3][i] + k[4][i])
    end
    return y
end

function stepFrame(f::Frame, foo::Function)
    y = f.y
    h = f.h
    t = f.time
    k = Vector{Vector{Float64}}(undef, 4)
    k[1] = h .* dyn(y, t)
    k[2] = h .* dyn(y .+ 0.5 .* k[1], t + 0.5 * h)
    k[3] = h .* dyn(y .+ 0.5 .* k[2], t + 0.5 * h)
    k[4] = h .* dyn(y .+ k[3], t + h)

    return y .+ (1.0/6.0) .* (k[1] .+ 2 .*k[2] .+ 2 .*k[3] .+ k[4])
end