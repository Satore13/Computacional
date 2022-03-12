include("Verlet.jl")
Particle, Frame = Verlet.Particle, Verlet.Frame
import Base.string
using CSV, CairoMakie, DataFrames
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
    length_of_sim = 10
    step = 0.0001
end
global PlanetsIndex = DataStructures.OrderedDict{Symbol, Int64}(:Mercurio => 1, :Venus=> 2, :Tierra => 3, :Marte => 4, :Jupiter => 5, :Saturno => 6, :Urano => 7, :Neptuno => 8)

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
    initial_frame = Frame(0.0, list_of_planets)
    simulation_data = Verlet.dataFrameRowFromFrame(initial_frame)
    
    time_between_frames = FPS^-1
    time_of_last_frame = 0.0
    number_of_frames = length_of_sim/step

    f = initial_frame
    
    println("Simulación del sistema solar con parámetros: ")
    println("FPS: ", FPS)
    println("Paso: ", step)
    println("Duración: ", length_of_sim)
    println("Frames totales a simular: ", number_of_frames)
    println("-------------------")
    last_percentage = 0.0
    for i in 1:number_of_frames
        f = Verlet.stepFrame(f, step, gravitationalAcc)
        current_percentage = (i/number_of_frames*100)
        if current_percentage - last_percentage > 1
            print("Progreso: ", round(current_percentage, digits = 2), "%\r")
            last_percentage = current_percentage
        end
        if f.t - time_of_last_frame > time_between_frames
            Verlet.addFrameRowToDataFrame!(simulation_data, f)
            time_of_last_frame = f.t
        end
    end
    println("Progreso: ", 100.00, "%")
    println("Completado")
    CSV.write("Tarea1/SistemaSolar.out", simulation_data)
end

begin
    #Parámetros de la animación
    width = (-2, 2)
    height = (-2, 2)
end
function buildAnimation()
    df = DataFrame(CSV.File("Tarea1/SistemaSolar.out"))
    
end
