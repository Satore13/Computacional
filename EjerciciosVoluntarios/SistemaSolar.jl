include("Verlet.jl")
Particle, Frame = Verlet.Particle, Verlet.Frame
import Base.string
using CSV
using DataStructures

#=
Funciones auxiliares para ayudar con las conversiones de unidades - Se asume que las unidades externas son el SIU - 
Las unidades de masa son masas solares, las de distancia son la distancia tierra sol y las de tiempo corresponden a ≈ 58.1 días
=#

toOwnUnitsMass(m::Float64) = m/1.989e30 
toOwnUnitsTime(t::Float64) = t * √(1.989e30*6.67384e-11/((1.496e11)^3))
toOwnUnitsMeters(d::Float64) = d/1.496e11 
fromOwnUnitsMass(m::Float64) = m * 1.989e30 
fromOwnUnitsTime(t::Float64) = t / √(1.989e30*6.67384e-11/((1.496e11)^3))
fromOwnUnitsMeters(d::Float64) = d * 1.496e11 



function norm(x::Vector{Float64})::Float64
    normsq = 0.0
    for xi in x
        normsq += xi^2
    end
    return √(normsq)
end
function gravitationalAcc(p::Particle, f::Frame)
    m_sol = 1
    
    #Interacciones entre los planetas
    total_f = 0.0
    total_f = (sum(f.particles) do other_p
        if(other_p ≢ p)
            return (1 / norm(p.x - other_p.x)^3) * (-p.x)
        end
        return zeros(length(p.x))
    end)
    #Interacción planeta-sol
    total_f += (1 / norm(p.x)^3) .* (-p.x)
    
    return total_f
end

# Parametros de la simulacion
FPS = 24

length_of_sim = 5
step = 0.001

global PlanetsIndex = DataStructures.OrderedDict{Symbol, Int64}(:Mercurio => 1, :Venus=> 2, :Tierra => 3, :Marte => 4, :Jupiter => 5, :Saturno => 6, :Urano => 7, :Neptuno => 8)

# Masas de los planetas en unidades de 10^24 metros que pasamos por nuestra función de conversión
global PlanetsMasses = toOwnUnitsMass.([0.330, 4.87, 5.97, 0.073, 0.642, 1899, 568, 86.8, 102, 0.0125].*1e24)
    
particleFromPlanet(i::Int64) = Particle([1.0*i, 0.0],[0.0, i*1.0], PlanetsMasses[i])
global list_of_planets = [particleFromPlanet(index) for index in values(PlanetsIndex)]::Vector{Particle}




function run()
    initial_frame = Frame(0.0, list_of_planets)
    simulation_data = Verlet.dataFrameRowFromFrame(initial_frame)
    
    time_between_frames = FPS^-1
    time_of_last_frame = 0.0
    number_of_frames = length_of_sim/step

    f = initial_frame
    
    for i in 1:number_of_frames
        f = Verlet.stepFrame(f, step, gravitationalAcc)
        if f.t - time_of_last_frame > time_between_frames
            Verlet.addFrameRowToDataFrame!(simulation_data, f)
            time_of_last_frame = f.t
        end
    end
    CSV.write("datos_sistema_solar.out", simulation_data)
end


