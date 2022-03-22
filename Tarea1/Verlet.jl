module Verlet
    #= AbstractFrame debería heredar un array de partículas que se consiga mediante getParticles así como una variable t =#
    abstract type AbstractFrame end
    #= Ha de heredar vectoresf64 x,v,a así como massf64 y tagsDict(Symbol => Symbol) =#
    abstract type AbstractParticle end
    using DataFrames, DataStructures
    #Estructura que almacena toda la información necesaria de una partícula
    struct Particle <: AbstractParticle
        x::Vector{Float64}
        v::Vector{Float64}
        a::Vector{Float64}
        mass::Float64
        tags::Dict{Symbol, Symbol}
        Particle(x::Vector{Float64}, v::Vector{Float64}, mass::Float64) = new(x, v,  zeros(Real, length(x)) ,mass, Dict())
        Particle(x::Vector{Float64}, v::Vector{Float64}, mass::Float64, tags::Dict{Symbol, Symbol}) = new(x, v,  zeros(Real, length(x)) ,mass, tags)
        Particle(x, v, a, mass) = new(x, v, a, mass)
    end

    #Microestado del sistema en un tiempo t, frame más básico
    struct BFrame <: AbstractFrame
        t::Float64
        particles::Vector{AbstractParticle}
    end

    getParticles(f::AbstractFrame)::Vector{AbstractParticle} = f.particles

    function getParticleWithTag(particles::Vector{P}, identifier::Pair{Symbol, Symbol})::Vector{P} where P <: AbstractParticle
        found_particles = P[]
        for particle in particles
            if identifier in particle.tags
                push!(found_particles, particle)
            end
        end
        return found_particles
    end

    function getParticleWithTag(f::T, identifier::Pair{Symbol, Symbol})::Vector{AbstractParticle} where T <: AbstractFrame
        return getParticleWithTag(f.particles, identifier)
    end

    function dataFrameRowFromFrame(frame::T)::DataFrame where T <: AbstractFrame
        df = DataFrame()
        insertcols!(df, "time" => [frame.t])
        for (id, p) in enumerate(getParticles(frame))
            for (i, x) in enumerate(p.x)
                insertcols!(df, ("p"*string(id)*"x"*string(i)) => [x])
            end
            for (i, v) in enumerate(p.v)
                insertcols!(df, ("p"*string(id)*"v"*string(i)) => [v])
            end
        end
        return df
    end
    function addFrameRowToDataFrame!(df::DataFrame, frame::T)::DataFrame where T <: AbstractFrame
        append!(df, dataFrameRowFromFrame(frame))
    end

    function stepFrame(frame::T, step::Float64, acceleration::Function, update::Function)::T where {T <: AbstractFrame, P<:AbstractParticle}
        aux_particles = Vector{P}(undef, 0)
        new_particles = Vector{P}(undef, 0)

        for p in getParticles(frame)
            initial_a = acceleration(p, frame)
            aux_w = p.v .+ step/2 .* initial_a
            new_x = p.x .+ step .* aux_w
            push!(aux_particles, P(new_x, p.v, initial_a,p.mass))
        end

        if update ≢ nothing
            aux_frame = update(frame)::T
        end
        aux_frame = T(frame.t + step, aux_particles)
        for p in getParticles(aux_frame)
            new_a = acceleration(p, aux_frame)::Vector{Float64}
            aux_w = p.v .+ step/2 .* p.a
            new_v = aux_w .+ step/2 .* new_a
            push!(new_particles, Particle(p.x, new_v, p.mass))
        end

        return update(T(frame.t + step, new_particles))::T
    end

    function stepFrame(frame::T, step::Float64, acceleration::Function)::T where T <: AbstractFrame
        stepFrame(frame, step, acceleration, (x) -> x)
    end

    #Funciones mecánicas adicionales
    calculateMomentum(p::Particle)::Vector{Float64} = p.v * p.mass
    calculateMomentum(frame::AbstractFrame)::Vector{Float64} = sum(calculateMomentum, getParticles(frame)) 
end