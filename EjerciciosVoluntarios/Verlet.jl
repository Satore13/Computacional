using DataFrames

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

struct Particle
    x::Vector{Real}
    v::Vector{Real}
    a::Vector{Real}
    mass::Real
    Particle(x, v, mass) = new(x, v,  zeros(Real, length(x)) ,mass)
    Particle(x, v, a, mass) = new(x, v, a, mass)

end

struct Frame
    t::Real
    particles::Vector{Particle}
end

function stepFrame(frame::Frame, step::Real)::Frame
    aux_particles = Vector{Particle}(undef, 0)
    new_particles = Vector{Particle}(undef, 0)

    for p in frame.particles
        initial_a = calculateForceOnParticle(p, frame) ./ p.mass
        aux_w = p.v .+ step/2 .* initial_a
        new_x = p.x .+ step .* aux_w
        push!(aux_particles, Particle(new_x, p.v, initial_a,p.mass))
    end
    aux_frame = Frame(frame.t + step, aux_particles)
    for p in aux_frame.particles
        new_a = calculateForceOnParticle(p, aux_frame) ./ p.mass
        aux_w = p.v .+ step/2 .* p.a
        new_v = aux_w .+ step/2 .* new_a
        push!(new_particles, Particle(p.x, new_v, p.mass))
    end

    return Frame(frame.t + step, new_particles)
end



function calculateForceOnParticle(p::Particle, f::Frame)::Vector{Real}
    total_f = zeros(Real, length(p.x))
    for other_p in f.particles
        #Skip calculating force with itself
        other_p === p && continue

        total_f = (other_p.x - p.x)
    end
    #total_f += [9.8] * p.mass
    return total_f
end

function run()
    f = Frame(0.0, [Particle([10.0*i], [0.0],1.0) for i in 1:2])
    simulation_data = dataFrameRowFromFrame(f)
    for i in 1:1000
        f = stepFrame(f, 0.01)
        addFrameRowToDataFrame!(simulation_data, f)
    end
    show(simulation_data)
end
