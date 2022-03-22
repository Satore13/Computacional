include("Verlet.jl")
import Base.rand
Particle, Frame = Verlet.Particle, Verlet.IndistFrame
using Serialization, CairoMakie
#=
Funciones auxiliares para ayudar con las conversiones de unidades - Se asume que las unidades externas son del SIU - 
Las unidades de masa son masas solares, las de distancia son la distancia tierra sol y las de tiempo corresponden a ≈ 58.1 días
En este sistema G = 1
=#

toOwnUnitsMass(m::Float64) = m/1.989e30 
toOwnUnitsTime(t::Float64) = t * √(1.989e30*6.67384e-11/((1.496e11)^3))
toOwnUnitsMeters(d::Float64) = d/1.496e11 
toOwnUnitsSpeed(v::Float64) = v/(1.496e11 *√(1.989e30*6.67384e-11/((1.496e11)^3)))#toOwnUnitsMeters(fromOwnUnitsTime(v))
fromOwnUnitsMass(m::Float64) = m * 1.989e30 
fromOwnUnitsTime(t::Float64) = t / √(1.989e30*6.67384e-11/((1.496e11)^3))
fromOwnUnitsMeters(d::Float64) = d * 1.496e11 
fromOwnUnitsSpeed(v::Float64) = fromOwnUnitsMeters(toOwnUnitsTime(v))

Range2f = Tuple{Float64, Float64}

struct ParticleDistribution
    energy_range::Range2f
    r_range::Range2f
    θ_range::Range2f
    mass_range::Range2f
    function ParticleDistribution(;energy_range::Range2f, r_range::Range2f, mass::Float64)
        if energy_range[1] == -Inf 
            return new((-mass/r_range[2], energy_range[2]), r_range, (0, 2*π), (mass, mass))
        end
        if energy_range[1] < -mass/r_range[2] 
            error("El rango inferior de energía no puede ser menor que -min(masa)/max(r)")
        end
        return new(energy_range, r_range, (0, 2*π), (mass, mass))
    end
end

function getRandBetween((a, b)::Range2f)::AbstractFloat
    length = b - a
    return rand()*length + a
end

function generateRandomParticle(dist::ParticleDistribution, type::Symbol)::Particle
    energy = getRandBetween(dist.energy_range)
    r = getRandBetween(dist.r_range)
    θ = getRandBetween(dist.θ_range)
    mass = getRandBetween(dist.mass_range)
    v = sqrt(2 * energy/ mass + 2/r)
    
    v_vector = v .* [-sin(θ), cos(θ)]
    x_vector = r .* [cos(θ), sin(θ)]

    return Particle(x_vector, v_vector, mass, Dict(:Tipo => type))
end

function norm²(x::Vector{Float64})::Float64
    normsq = 0.0
    for xi in x
        normsq += xi^2
    end
    return normsq
end
function norm(x::Vector{Float64})::Float64
    return √(norm²(x))
end
function gravitationalAcc(p::Particle, f::Frame)::Vector{Float64}
    #Interacción planeta-sol, para evitar inestabilidad toda interacción a menos de 0.01 del sol es capada a 0.01
    total_f = (1 / norm(p.x)^3) .* (-p.x)
    return total_f
end

function totalEnergyOfParticle(p::Particle)::Float64
    kEnergy(p) + uEnergy(p)
end

function energyOfFrame(f::Frame)::Float64
    sum(f.particles) do p
        totalEnergyOfParticle(p)
    end
end

function kEnergy(p::Particle)::Float64
    0.5 * norm²(p.v) * p.mass
end
function uEnergy(p::Particle)::Float64  
    -p.mass/norm(p.x)
end

#sframes son los frames a grabar y rframes son los frames simulados -no se graban, sólo se usan para aumentar la precisión-
begin
    global FPS = 10
    global rframes_per_sframes = 100
    global length = 1000
    global N_of_particles = 100
end


function run()
    pd = ParticleDistribution(energy_range = (-Inf, 0.0), r_range = (0.5, 1.0), mass = 1e-6)
    planetesimals = Particle[generateRandomParticle(pd, :Sólido) for _ in 1:N_of_particles]
    initial_frame = Frame(0.0, planetesimals)
    step = 1/(FPS*rframes_per_sframes)
    f = initial_frame
    frames = Frame[initial_frame]
    println("Duración: ", length * 58.1 /365, " años")
    
    percentage = Observable(0.0)
    on(percentage) do percentage
        print("Porcentaje: ", round(percentage, digits = 1), "%\r")
    end
    for nsf in 1:(FPS*length)
        for n in 1:rframes_per_sframes
            f = Verlet.stepFrame(f, step, gravitationalAcc)
        end
        push!(frames, f)
        if nsf/(FPS*length)*100 > percentage[] + 0.1
            percentage[] = nsf/(FPS*length)*100
        end
    end
    percentage[] = 100.0

    serialize("Tarea1/CreacionSistemaSolar/sim.out", frames)
    return nothing
end

function plotError()

    data = deserialize("Tarea1/CreacionSistemaSolar/sim.out")
    initial_energy = energyOfFrame(first(data))
    figure = Figure()
    calcError(f::Frame) = abs((energyOfFrame(f) - initial_energy) / initial_energy)
    errors = [calcError(f) for f in data]
    times_in_years = [f.t * 58.1/365 for f in data]   
    ax = Axis(figure[1,1], title = "Error relativo en la energía")

    lines!(ax, times_in_years, errors)
    return figure
end