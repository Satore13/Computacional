module Verlet
    using DataFrames
    struct Particle
        x::Vector{Float64}
        v::Vector{Float64}
        a::Vector{Float64}
        mass::Float64
        Particle(x, v, mass) = new(x, v,  zeros(Real, length(x)) ,mass)
        Particle(x, v, a, mass) = new(x, v, a, mass)

    end

    struct Frame
        t::Float64
        particles::Vector{Particle}
    end

    function dataFrameRowFromFrame(frame::Frame)::DataFrame
        df = DataFrame()
        insertcols!(df, "time" => [frame.t])
        for (id, p) in enumerate(frame.particles)
            for (i, x) in enumerate(p.x)
                insertcols!(df, ("p"*string(id)*"x"*string(i)) => [x])
            end
            for (i, v) in enumerate(p.v)
                insertcols!(df, ("p"*string(id)*"v"*string(i)) => [v])
            end
        end
        return df
    end
    function addFrameRowToDataFrame!(df::DataFrame, frame::Frame)::DataFrame
        append!(df, dataFrameRowFromFrame(frame))
    end

    function stepFrame(frame::Frame, step::Float64, force::Function, update::Function)::Frame
        aux_particles = Vector{Particle}(undef, 0)
        new_particles = Vector{Particle}(undef, 0)

        for p in frame.particles
            initial_a = force(p, frame) ./ p.mass
            aux_w = p.v .+ step/2 .* initial_a
            new_x = p.x .+ step .* aux_w
            push!(aux_particles, Particle(new_x, p.v, initial_a,p.mass))
        end

        aux_frame = update(frame)::Frame

        aux_frame = Frame(frame.t + step, aux_particles)
        for p in aux_frame.particles
            new_a = force(p, aux_frame)::Vector{Float64} ./ p.mass
            aux_w = p.v .+ step/2 .* p.a
            new_v = aux_w .+ step/2 .* new_a
            push!(new_particles, Particle(p.x, new_v, p.mass))
        end

        return update(Frame(frame.t + step, new_particles))::Frame
    end

    function stepFrame(frame::Frame, step::Float64, force::Function)::Frame
        stepFrame(frame, step, force, (x) -> x)
    end


    function calculateForceOnParticle(p::Particle, f::Frame)::Vector{Float64}
        total_f = zeros(Real, length(p.x))
        for other_p in f.particles
            #Skip calculating force with itself
            other_p === p && continue
        end
        total_f += [-p.mass*9.8*sin(p.x[1])]
        #total_f += [9.8] * p.mass
        return total_f
    end
end