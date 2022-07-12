include("RungeKutta.jl")
using Serialization
using ForwardDiff
using GeometryBasics

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

function crear_simulacion(;θ1::Float64, θ2::Float64, l1::Float64 = 1.0, l2::Float64 = 1.0, m1::Float64 = 1.0, m2::Float64 = 1.0, g::Float64 = 10.0, h::Float64 = 1e-4, fps::Float64 = 30.0)
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

function angles_velocity(f::Frame, parameters::Dict)
    y = f.y
    θ1 = y[1]
    θ2 = y[2]
    p1 = y[3]
    p2 = y[4]
    l1 = parameters[:l1]
    l2 = parameters[:l2]
    m1 = parameters[:m1]
    m2 = parameters[:m2]
    g = parameters[:g]
    
    
    vθ1 = (l2*y[3] - l1*y[4]*cos(y[1] - y[2]))/ (l1^2*l2*(m1 + m2 * sin(y[1] - y[2])^2)) 
    vθ2 = (-m2*l2*y[3]*cos(y[1] - y[2]) +(m1 + m2) *l1*y[4] ) / (m2 * l1 * l2^2 * (m1 + m2*sin(y[1] - y[2])^2))
    
    return (vθ1, vθ2)
end
function energia_1(f::Frame, parameters::Dict)
    y = f.y
    θ1 = y[1]
    θ2 = y[2]
    p1 = y[3]
    p2 = y[4]
    l1 = parameters[:l1]
    l2 = parameters[:l2]
    m1 = parameters[:m1]
    m2 = parameters[:m2]
    g = parameters[:g]
    
    
    vθ1 = (l2*y[3] - l1*y[4]*cos(y[1] - y[2]))/ (l1^2*l2*(m1 + m2 * sin(y[1] - y[2])^2)) 
    vθ2 = (-m2*l2*y[3]*cos(y[1] - y[2]) +(m1 + m2) *l1*y[4] ) / (m2 * l1 * l2^2 * (m1 + m2*sin(y[1] - y[2])^2))
    
    return 0.5 * m1 * l1^2 * vθ1^2 - m1 *g * l1 * cos(θ1)

end

function energia_2(f::Frame, parameters::Dict)
    y = f.y
    θ1 = y[1]
    θ2 = y[2]
    p1 = y[3]
    p2 = y[4]
    l1 = parameters[:l1]
    l2 = parameters[:l2]
    m1 = parameters[:m1]
    m2 = parameters[:m2]
    g = parameters[:g]
    
    
    vθ1 = (l2*y[3] - l1*y[4]*cos(y[1] - y[2]))/ (l1^2*l2*(m1 + m2 * sin(y[1] - y[2])^2)) 
    vθ2 = (-m2*l2*y[3]*cos(y[1] - y[2]) +(m1 + m2) *l1*y[4] ) / (m2 * l1 * l2^2 * (m1 + m2*sin(y[1] - y[2])^2))
    
    return 0.5 * m2 *  (l1^2 * vθ1^2 + l2^2 * vθ2^2 + 2 * l1 * l2 * vθ1 * vθ2 * cos(θ1 - θ2)) - m2 * g * (l1 * cos(θ1) + l2 * cos(θ2))

end

function wrap_angles!(frame::Frame, p::Dict{Symbol, Float64})
    frame.y[1] = mod2pi.(frame.y[1] .+ π) .- π
    frame.y[2] = mod2pi.(frame.y[2] .+ π) .- π
    return frame
end


function batch_lyap()
    ranges = [ 5:5:60, 65:5:120, 125:5:180]
    
    for (i, range) in enumerate(ranges)
        range = collect(range)
        @show range
        batch = [crear_simulacion(θ1 = x, θ2 = x) for x in deg2rad.(range)]
        ejecutar_conjunto_simulaciones(batch, "lyapunov/batch$(i).out", 2000.0)
        println("Ejecutado batch $(i)")
    end
    for (i, range) in enumerate(ranges)
        range = collect(range)
        @show range
        batch = [crear_simulacion(θ1 = x , θ2 = x) for x in deg2rad.(range .+ 0.001)]
        ejecutar_conjunto_simulaciones(batch, "lyapunov/batch$(i)alter.out", 2000.0)
        println("Ejecutado batch $(i)")
    end
end

#https://test-sprott.physics.wisc.edu/chaos/lyapexp.htm
function calculate_lyapunov(sim::Simulation, steps::Integer)
    d0 = 10^-8 #Raíz cuadrada de la precisión

    time_local_exp = Float64[]

    trayA = sim.video[1]
    trayB = Frame(0.0, trayA.y .+ d0 * sqrt(1/4))

    diff_tray(A, B) = sqrt((A.y[1] - B.y[1]) ^ 2  + 
                                (A.y[2] - B.y[2])^2 + 
                                (A.y[3] - B.y[3])^2 + 
                                (A.y[4] - B.y[4])^2)

    porc_previo = 0.0
    suma = 0.0
    lyap_vs_time = Point2f[]
    for i in 1:steps
        trayA = stepFrame(trayA, sim.dyn, sim.h, sim.parameters)
        trayB = stepFrame(trayB, sim.dyn, sim.h, sim.parameters)

        d1 = diff_tray(trayA, trayB)
        local_exp = log2(abs(d1 / d0)) / sim.h

        #Actualizar trayectorias
        By1 = trayA.y[1] + d0/d1 * (trayB.y[1] - trayA.y[1])
        By2 = trayA.y[2] + d0/d1 * (trayB.y[2] - trayA.y[2])
        By3 = trayA.y[3] + d0/d1 * (trayB.y[3] - trayA.y[3])
        By4 = trayA.y[4] + d0/d1 * (trayB.y[4] - trayA.y[4])
        trayB = Frame(trayA.time, [By1, By2, By3, By4])

        suma += local_exp
        porc_act = i / steps * 100
        if porc_act - porc_previo >= 0.1
            print("Porcentaje actual: $(round(porc_act, digits = 1)) -- Media actual: $(suma / i) \t\t\r")
            push!(lyap_vs_time, Point2f(trayA.time, suma / i))
            porc_previo = porc_act
        end
    end
    println("Porcentaje actual: $(round(100.0, digits = 1)) -- Media actual: $(suma / steps) \t\t")
    push!(lyap_vs_time, Point2f(trayA.time, suma / steps))
    return lyap_vs_time
end

function batches_lyapunov()
    for θ₀ in 5:5:175
        x = deg2rad(θ₀)
        sim = crear_simulacion(θ1 = x, θ2 = x)
        println("Ejecutando simulación θ1 = θ2 = $(θ₀)")
        a = calculate_lyapunov(sim, 1_000 * 100)
        filename = "Tarea4/outputPD/truelyap/same_angle/$(θ₀).out"
        serialize(filename, a)
        println(prod(["-" for _ in 1:100]))
        sim = crear_simulacion(θ1 = 0.0, θ2 = x)
        println("Ejecutando simulación θ1 = 0, θ2 = $(θ₀)")
        a = calculate_lyapunov(sim, 1_000 * 100)
        filename = "Tarea4/outputPD/truelyap/theta1=0/$(θ₀).out"
        serialize(filename, a)
        println(prod(["-" for _ in 1:100]))
    end
end

function calculate_frecuencias(sim::Simulation, steps::Integer)
    tray = sim.video[1]
    frecuencias = Tuple{Float64, Float64}[]
    signo_ant1 = 1
    signo_ant2 = 1

    tiempo_ant = 0.0
    porc_previo = 0.0
    @show sim.h
    for i in 1:steps
        tray = stepFrame(tray, sim.dyn, sim.h, sim.parameters)
        wrap_angles!(tray, sim.parameters)
        w = angles_velocity(tray, sim.parameters)[2]
        p2 = tray.y[4]
        θ1 = tray.y[1]
        θ2 = tray.y[2]
        if signo_ant1 != sign(θ1)  #signo_ant2 != sign(θ2)
            push!(frecuencias, (abs(θ2), p2))#1/(2 * (tray.time - tiempo_ant)))
            signo_ant1 = sign(θ1)
        end
        porc_act = i / steps * 100
        if porc_act - porc_previo >= 5
            print("Porcentaje actual: $(round(porc_act, digits = 1)) -- Encontradas $(length(frecuencias)) frecuencias \t\t\r")
            porc_previo = porc_act
        end
    end
    sort!(frecuencias, by = (v) -> v[1])
    if length(frecuencias) < 20
        return getindex.(frecuencias, 2)
    end
    frecuencias = frecuencias[1:20]
    return getindex.(frecuencias, 2)
end

function batch_frecuencias_same_angle()
    freq_vs_ang = Tuple{Float64, Vector{Float64}}[]
    for θ₀ in 1:0.1:180
        x = deg2rad(θ₀)
        println("Frecuencias calculandose para θ₀ = $θ₀")
        ω = calculate_frecuencias(crear_simulacion(θ1 = x, θ2 = x, h = 1e-3), 100_000)
        println()
        println("Frecuencias calculada para θ₀ = $θ₀")
        push!(freq_vs_ang, (float(θ₀), ω))
    end
    serialize("Tarea4/outputPD/frecuencias/simetricas.out", freq_vs_ang)
    return freq_vs_ang
end

function batch_frecuencias_θ1__0()
    freq_vs_ang = Tuple{Float64, Vector{Float64}}[]
    for θ₀ in 1:0.1:180
        x = deg2rad(θ₀)
        println("Frecuencias calculandose para θ₀ = $θ₀")
        ω = calculate_frecuencias(crear_simulacion(θ1 = 0.0, θ2 = x, h = 1e-3), 100_000)
        println()
        println("Frecuencias calculada para θ₀ = $θ₀")
        push!(freq_vs_ang, (float(θ₀), ω))
    end
    serialize("Tarea4/outputPD/frecuencias/simetricas.out", freq_vs_ang)
    return freq_vs_ang
end


#Creamos carpetas necesarias si no existen
mkpath("Tarea4/outputPD")
mkpath("Tarea4/outputPD/frecuencias")
mkpath("Tarea4/outputPD/poincare")
mkpath("Tarea4/outputPD/lyapunov")
mkpath("Tarea4/outputPD/truelyap")
mkpath("Tarea4/outputPD/posiciones")
mkpath("Tarea4/analisis_error")


