include("RungeKutta.jl")
using Serialization

function dyn(y, t, parameters)
    #En orden θ1, θ2, p1, p2 = y1, y2, y3, y4

    θ1 = y[1]
    θ2 = y[2]
    p1 = y[3]
    p2 = y[4]
    #https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf
    l1 = parameters[:l1]
    l2 = parameters[:l2]
    m1 = parameters[:m1]
    m2 = parameters[:m2]
    g = parameters[:g]
    y_ = zero(y)
    h1 = (p1*p2 * sin(θ1 - θ2)) / (l1 * l2 * (m1 + m2 * sin(θ1 - θ2)^2))
    h2 = (m2 * l2^2 * p1^2 + (m1 + m2) * l1^2 * p2^2 - 2 * m2* l1*l2*p1*p2*cos(θ1 - θ2)) / (2 * (l1^2 * l2^2 * (m1 + m2 *sin(θ1 - θ2)^2))^2)

    
    
    y_[1] = (l2*y[3] - l1*y[4]*cos(y[1] - y[2]))/ (l1^2*l2*(m1 + m2 * sin(y[1] - y[2])^2)) 
    y_[2] = (-m2*l2*y[3]*cos(y[1] - y[2]) +(m1 + m2) *l1*y[4] ) / (m2 * l1 * l2^2 * (m1 + m2*sin(y[1] - y[2])^2))
    y_[3] = -(m1 + m2)*g*l1*sin(y[1])  - h1 + h2 * sin(2*(y[1]-y[2]))
    y_[4] = -m2*g*l2*sin(y[2]) + h1 - h2 * sin(2(y[1] - y[2]))

    
    return y_
end

function crear_simulacion(;θ1::Float64, θ2::Float64, l1::Float64 = 1.0, l2::Float64 = 1.0, m1::Float64 = 1.0, m2::Float64 = 1.0, g::Float64 = 10.0, h::Float64 = 1e-4, fps::Float64 = 30)
    Simulation([:θ1, :θ2, :p1, :p2], [θ1, θ2, 0.0, 0.0], dyn, fps, parameters = Dict([:m1 => m1, :m2 => m2, :l1 => l1, :l2 => l2, :g => g]), h = h)
end

#Funcion para pasar de las posiciones generalizadas a coordenadas cartesianas
function r1r2_de_polares(frame::Frame, parameters)

    θ1, θ2 = frame.y[1], frame.y[2]
    x1 = sin(θ1) * parameters[:l1]
    y1 = -cos(θ1) * parameters[:l1]

    x2 = sin(θ2) * parameters[:l2] + x1
    y2 = -cos(θ2) * parameters[:l2] + y1

    return Point2f(x1, y1), Point2f(x2, y2)
end

function ejecutar_simulaciones_similares()
    θ_range = range(start = π/4, stop = 3*π/4, length = 100)
    for (index, θ_inicial) in enumerate(θ_range)
        filename = "Tarea4/outputPD/simulacion_liap$index.out"
        sim = crear_simulacion(θ1 = θ_inicial, θ2 = θ_inicial)
        println("Ejecutando simulación número $index con θ inicial $(round(θ_inicial, digits = 6))")
        loop!(sim, 1000.0, 0.1)
        serialize(filename, sim)
        #Ahora ejecutamos con valores iniciales muy próximos
        
        filename = "Tarea4/outputPD/simulacion_liap$(index)alter.out"
        sim = crear_simulacion(θ1 = θ_inicial + 1e-6 , θ2 = θ_inicial + 1e-6)
        println("Ejecutando simulación número $index con θ inicial $(round(θ_inicial + 1e-6, digits = 6))")
        loop!(sim, 1000.0, 0.1)
        serialize(filename, sim)
    end
end

#batch_analizar_tiempo = [crear_simulacion(θ1 = π*0.5, θ2 = π*0.5, h = x) for x in [1e-3, 1e-4, 1e-5, 1e-6]]

#batch_analizar_angulos = [crear_simulacion(θ1 = x, θ2 = x ,h = 1e-4) for x in deg2rad.(30:30:150)]

function ejecutar_conjunto_simulaciones(batch::Vector{Simulation}, filename::String, duration::Float64 = 10000.0)
    filetree = "Tarea4/outputPD/"
    for s in batch
        loopfxts!(s, wrap_angles!, duration)
    end
    serialize(filetree * filename, batch)
end


function energia(f::Frame, p::Dict)
    θ1 = f.y[1]
    θ2 = f.y[2]
    p1 = f.y[3]
    p2 = f.y[4]

    l1 = p[:l1]
    l2 = p[:l2]

    m1 = p[:m1]
    m2 = p[:m2]

    g = p[:g]

    a = m2 * l2^2 * p1^2 + (m1 + m2) * l1^2 * p2^2 - 2 * m2 * l1 * l2 * p1 * p2 * cos(θ1 - θ2)
    b = 2 * m2 * l1^2 * l2^2 * (m1 + m2 * (sin(θ1-θ2))^2)
    c = (m1 + m2) * g * l1 * cos(θ1) + m2 * g * l2 * cos(θ2)     
    return  a / b - c  +  (m1 + m2) * g*l1 + m2 *g*l2
end

function wrap_angles!(frame::Frame, p::Dict{Symbol, Float64})
    frame.y[1] = mod2pi.(frame.y[1] .+ π) .- π
    frame.y[2] = mod2pi.(frame.y[2] .+ π) .- π
end