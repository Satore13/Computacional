include("../Verlet.jl")
Frame, Particle = Verlet.BFrame, Verlet.BParticle
#= F(r) = 4ε/r * [12 (σ/ r)^12 - 6(σ / r)^6], tomamos σ = ε = 1

=#
using Serialization, CairoMakie
const L = 4 #Tamaño del receptáculo

function update(f::Frame)
    new_particles = Particle[]
    function wrap_pos(x)
        if x <= 0
            return x + L
        elseif x >= L
            return x - L
        else
            return x
        end
    end

    for P in f.particles
        
        new_x = [wrap_pos(xi) for xi in p.x]

        push!(new_particles, Particle(new_x, p.v, p.mass))
    end
end

function random_particle()::Particle
    mass = 1.0
    x = rand(Float64, 2) .* L
    v = [1.0, 0]
    return Particle(x, v, mass)
end

function run()
    step = 0.002
    population_size = 10
    simulation_length = 10
    number_of_real_frames = simulation_length / step

    initial_particles = Particle[]
    for _ in 1:population_size
        push!(initial_particles, random_particle())
    end
    return initial_particles
end