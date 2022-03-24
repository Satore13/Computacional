include("Verlet.jl")
Particle, Frame = Verlet.BParticle, Verlet.BFrame
using CSV, CairoMakie, Serialization
using DataStructures

#=
Funciones auxiliares para ayudar con las conversiones de unidades - Se asume que las unidades externas son del SIU - 
Las unidades de masa son masas solares, las de distancia son la distancia tierra sol y las de tiempo corresponden a ≈ 58.1 días
=#

toOwnUnitsMass(m::Float64) = m/1.989e30 
toOwnUnitsTime(t::Float64) = t * √(1.989e30*6.67384e-11/((1.496e11)^3))
toOwnUnitsMeters(d::Float64) = d/1.496e11 
toOwnUnitsSpeed(v::Float64) = v/(1.496e11 *√(1.989e30*6.67384e-11/((1.496e11)^3)))#toOwnUnitsMeters(fromOwnUnitsTime(v))
fromOwnUnitsMass(m::Float64) = m * 1.989e30 
fromOwnUnitsTime(t::Float64) = t / √(1.989e30*6.67384e-11/((1.496e11)^3))
fromOwnUnitsMeters(d::Float64) = d * 1.496e11 
fromOwnUnitsSpeed(v::Float64) = fromOwnUnitsMeters(toOwnUnitsTime(v))

function norm(x::Vector{Float64})::Float64
    normsq = 0.0
    for xi in x
        normsq += xi^2
    end
    return √(normsq)
end
function gravitationalAcc(p::Particle, f::Frame)::Vector{Float64}
    m_sol = 1
    
    #Interacciones entre los planetas
    total_f = zeros(length(p.x))
    total_f = (sum(f.particles) do other_p
        if(other_p ≢ p)
            return (1 / norm(p.x - other_p.x)^3) * (-p.x) * other_p.mass
        end
        return zeros(length(p.x))
    end) 
    #Interacción planeta-sol
    total_f += (1 / norm(p.x)^3) .* (-p.x)
    
    return total_f
end
    
# Parametros de la simulacion
begin
    FPS = 10
    length_of_sim = 100
    step = 0.01
end
global PlanetsIndex = DataStructures.OrderedDict{Symbol, Int64}(:Mercurio => 1, :Venus=> 2, :Tierra => 3, :Marte => 4, :Jupiter => 5, :Saturno => 6, :Urano => 7, :Neptuno => 8, :Pluton => 9)

begin
    # Masas de los planetas en unidades de 10^24 metros que pasamos por nuestra función de conversión
    global PlanetsMasses = toOwnUnitsMass.([0.330, 4.87, 5.97, 0.073, 0.642, 1899, 568, 86.8, 102, 0.0125].*1e24)

    #Distancias planeta-sol
    global PlanetsDistance = [0.387, 0.723, 1, 1.52, 5.20, 9.57, 19.17, 30.18, 39.48]

    #Velocidad orbital (km/s)
    global PlanetsSpeed = toOwnUnitsSpeed.([47.9, 35.0, 29.8, 24.1, 13.1, 9.7, 6.8, 5.4, 4.7].*1e3)

    particleFromPlanet(i::Int64) = Particle([PlanetsDistance[i] ; 0.0],[0.0; PlanetsSpeed[i]], PlanetsMasses[i])
    global list_of_planets = [particleFromPlanet(index) for index in values(PlanetsIndex)]::Vector{Particle}
end

function runSimulation()
    f = Frame(0.0, list_of_planets)
    
    time_between_frames = FPS^-1
    number_of_frames = length_of_sim/step

    frames = Frame[f]

    println("Simulación del sistema solar con parámetros: ")
    println("FPS: ", FPS)
    println("Paso: ", step)
    println("Duración: ", length_of_sim)
    println("Frames reales a simular: ", number_of_frames)
    println("Frames a guardar: ", FPS * length_of_sim)
    println("-------------------")
    last_percentage = 0.0
    for i in 1:number_of_frames
        f = Verlet.stepFrame(f, step, gravitationalAcc)
        current_percentage = (i/number_of_frames*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
        if f.t - last(frames).t > time_between_frames
            push!(frames, f)
        end
    end
    println("Progreso: ", 100.00, "%")
    println("Completado")
    serialize("Tarea1/SistemaSolar.out", frames)
end

#Mal hecho el tener dos funciones así
function buildAnimationGeoc(filename::String = "Tarea1/SistemaSolarGeoc.out")
    #Leer los datos y guardarlos en un DataFrame
    println("Cargando archivo...")
    data = deserialize(filename)
    println("Archivo cargado")
    current_frame_n = Observable{Int64}(1)

    figure = Figure()
    ax = Axis(figure[1,1], title = @lift("t = " * string(round(58.1*data[$current_frame_n].t, digits = 2)) * " días"), aspect = 1)
    ax_full = Axis(figure[2,1], aspect = 1)
    xlims!(ax, (-2.5, 2.5))
    ylims!(ax, (-2.5, 2.5))
    xlims!(ax_full, (-50, 50))
    ylims!(ax_full, (-50, 50))

    function fromPlanetSymbolToX1X2(planet_symbol::Symbol, row_n::Int64)
        return (data[row_n].particles[PlanetsIndex[planet_symbol]].x[1],
                    data[row_n].particles[PlanetsIndex[planet_symbol]].x[2])
    end


    
    for plabel in PlanetsIndex
        coord = @lift(fromPlanetSymbolToX1X2(plabel.first, $current_frame_n))
        point = @lift(Point2f[$(coord)])
        scatter!(ax, point , label = string(plabel.first))
        scatter!(ax_full, point , label = string(plabel.first))          
    end 
    sol_coord = @lift(( data[$current_frame_n].particles[10].x[1], data[$current_frame_n].particles[10].x[2]))
    sol_point = @lift(Point2f[$sol_coord])
    scatter!(ax, sol_point , label = "Sol")
    scatter!(ax_full, sol_point , label = "Sol")



    Legend(figure[:,2], ax)

    last_percentage = 0.0
    number_of_rows = size(data)[1]
    record(figure, "Tarea1/SistemaSolarGeoc.mp4", range(1, number_of_rows); framerate = FPS) do row_number
        current_frame_n[] = row_number
        current_percentage = (row_number/number_of_rows*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
    end
end

function buildAnimation(filename::String = "Tarea1/SistemaSolar.out")
    #Leer los datos y guardarlos en un DataFrame
    println("Cargando archivo...")
    data = deserialize(filename)
    println("Archivo cargado")
    current_frame_n = Observable{Int64}(1)

    figure = Figure()
    ax = Axis(figure[1,1], title = @lift("t = " * string(round(data[$current_frame_n].t *58.1, digits = 2)) * " días"), aspect = 1)
    ax_full = Axis(figure[2,1], aspect = 1)
    xlims!(ax, (-2.5, 2.5))
    ylims!(ax, (-2.5, 2.5))
    xlims!(ax_full, (-50, 50))
    ylims!(ax_full, (-50, 50))

    function fromPlanetSymbolToX1X2(planet_symbol::Symbol, row_n::Int64)
        return (data[row_n].particles[PlanetsIndex[planet_symbol]].x[1],
                    data[row_n].particles[PlanetsIndex[planet_symbol]].x[2])
    end


    sol_scatter = scatter!(ax, [Point2f(0,0)], label = "Sol")
    sol_scatter = scatter!(ax_full, [Point2f(0,0)], label = "Sol")


    for plabel in PlanetsIndex
        coord = @lift(fromPlanetSymbolToX1X2(plabel.first, $current_frame_n))
        point = @lift(Point2f[$(coord)])
        scatter!(ax, point , label = string(plabel.first))
        scatter!(ax_full, point , label = string(plabel.first))          
    end 
   
    Legend(figure[:,2], ax)

    last_percentage = 0.0
    number_of_rows = size(data)[1]
    record(figure, "Tarea1/SistemaSolar.mp4", range(1, number_of_rows); framerate = FPS) do row_number
        current_frame_n[] = row_number
        current_percentage = (row_number/number_of_rows*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
    end
end

function totalEnergyFrame(f::Frame)
    total_k = sum([kEnergy(p) for p in f.particles])
    total_p = sum([uEnergy(p) for p in f.particles])

    return total_k + total_p
end

function plotEnergyvsTime()
    data = deserialize("Tarea1/SistemaSolar.out")

    t_ = getfield.(data, :t)
    E_ = [totalEnergyFrame(f) for f in data]

    relErr_ = abs.( (first(E_) .- E_) ./ first(E_))

    fig = Figure()
    ax = Axis(fig[1,1], title = "Error relativo en la energía")
    lines!(ax, t_, relErr_)

    return fig
end

function plotAngularMomentum()
    data = deserialize("Tarea1/SistemaSolar.out")

    t_ = getfield.(data, :t)
    L_ = [AngularMomentum(f) for f in data]

    relErr_ = abs.( (first(L_) .- L_) ./ first(L_))

    fig = Figure()
    ax = Axis(fig[1,1], title = "Error relativo en el Momento Angular")
    lines!(ax, t_, relErr_)

    return fig
end

function AngularMomentum(f::Frame)
    return sum([AngularMomentum(p) for p in f.particles])
end

function AngularMomentum(p::Particle)
    crossP = abs( p.x[1] * p.v[2] - p.x[2]*p.v[1])
    return p.mass * crossP
end

function kEnergy(p::Particle)::Float64
    0.5 * norm(p.v)^2 * p.mass
end
function uEnergy(p::Particle)::Float64  
    -p.mass/norm(p.x)
end

function createGeocView()
    data = deserialize("Tarea1/SistemaSolar.out")

    new_data = Frame[]
    for f in data
        new_t = f.t
        new_particles = Particle[]
        for p in f.particles
            new_v = p.v .- f.particles[PlanetsIndex[:Tierra]].v
            new_mass = p.mass
            new_x = p.x .- f.particles[PlanetsIndex[:Tierra]].x
            push!(new_particles, Particle(new_x, new_v, new_mass))
        end
        sol_v = (-1) .* f.particles[PlanetsIndex[:Tierra]].v
        sol_x = (-1) .* f.particles[PlanetsIndex[:Tierra]].x
        push!(new_particles, Particle(sol_x, sol_v, 1.0))
        push!(new_data, Frame(new_t, new_particles))
    end

    serialize("Tarea1/SistemaSolarGeoc.out", new_data)
end

function runAndRender()
    runSimulation()
    buildAnimation()
    createGeocView()
    buildAnimationGeoc()
end