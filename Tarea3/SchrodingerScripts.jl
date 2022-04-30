include("Schrodinger.jl")
using GLMakie
using Serialization

function plotear_probabilidad!(ax, φ::WaveFunction, color = :blue)
    N = size(φ) - 1
    x = collect(0:N)
    y = [probability_density(φ, xi) for xi in x]
    @show y
    lines!(ax, x, y, color = color)
end

function plotear_probabilidad(φ::WaveFunction)
    ax = Axis(Figure()[1, 1])
    plotear_probabilidad!(ax, φ)
    current_figure()
end

function plotear_anim(sim::Simulation, filename::String="Tarea3/videos/sim.mp4")
    N = sim.N
    x = collect(0:N)
    index = Observable{Int64}(1)
    φn = @lift(sim.φn[$index])    
    energy_op = build_energy_op(build_potential_op(sim.V))

    fig = Figure()

    ax_wf = Axis(fig[1, 1:3], xlabel = "x", yaxisposition = :left)
    ax_potential = Axis(fig[1, 1:3], yaxisposition = :right)
    hidespines!(ax_potential)
    hidexdecorations!(ax_potential)
    ρ = @lift(Vector0([probability_density($φn, i) for i in 0:N]))
    density_dots = @lift([Point2f(i, $ρ[i]) for i in 0:N])
    plt_density = lines!(ax_wf, density_dots)
    mean_x = @lift(mean_value($φn, position_op))
    plt_mean_x = vlines!(ax_wf, mean_x, color = :black)
    potential_dots = [Point2f(j, sim.V[j]) for j in 0:N]
    plt_potential = lines!(ax_potential, potential_dots, color = :red)
    plt_energy = hlines!(ax_potential, mean_value(sim.φn[0], energy_op), color = :orange)
    xlims!(ax_wf, (0, N))
    xlims!(ax_potential, (0, N))
    Legend(fig[1, 4], [plt_density, plt_mean_x, plt_potential, plt_energy], ["Densidad de probabilidad", "<x>", "V(x)", "<H>"])

    ax_parts = Axis(fig[2, 1:3], xlabel = "x")
    imag_dots = @lift([Point2f(i, imag($φn.val[i])) for i in 0:N])
    real_dots = @lift([Point2f(i, real($φn.val[i])) for i in 0:N])
    lines!(ax_parts, imag_dots, label = "Parte imaginaria")
    lines!(ax_parts, real_dots, label = "Parte real")
    xlims!(ax_parts, (0, N))
    Legend(fig[2, 4], ax_parts)

    ax_norm = Axis(fig[4, 1:3], xlabel = "t")
    norm_dots = Observable(Point2f[(-1,0)])
    lines!(ax_norm, norm_dots, label = "|φ|²")
    ylims!(0, 1.1)
    xlims!(0, eachindex(sim.φn)[end])
    Legend(fig[4, 4], ax_norm)

    ax_momentum = Axis(fig[3, 1:3], xlabel = "t")
    momentum_dots = Observable(Point2f[(0, sim.k0)])
    energy_dots = Observable(Point2f[(-1, 0)])
    lines!(ax_momentum, momentum_dots, label = "<p>")
    lines!(ax_momentum, energy_dots, label = "<H>")
    xlims!(ax_momentum, (0, eachindex(sim.φn)[end]))
    ylims!(ax_momentum, (-sim.k0, sim.k0))
    Legend(fig[3, 4], ax_momentum)


    record(fig, filename, eachindex(sim.φn),framerate = 20) do i
        index[] = i
        norm_dots[] = push!(norm_dots[], Point2f(i, norm(sim.φn[i])))
        momentum_dots[] = push!(momentum_dots[], Point2f(i, mean_value(sim.φn[i], momentum_op)))
        energy_dots[] = push!(energy_dots[], Point2f(i, mean_value(sim.φn[i], energy_op)))
    end
end