include("RungeKutta.jl")
using GLMakie

function dyn(y, t, parameters)
    #En orden θ1, θ2, p1, p2 = y1, y2, y3, y4

    #https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf
    l1 = parameters[:l1]
    l2 = parameters[:l2]
    m1 = parameters[:m1]
    m2 = parameters[:m2]
    g = parameters[:g]
    y_ = zero(y)
    h1 = (y[3]*y[4] * sin(y[1] - y[2])) / (l1 * l2 * (m1 + m2 * sin(y[1] - y[2])^2))
    h2 = (m2 * l2^2 * y[3]^2 + (m1 + m2) * l1^2 * y[4]^2 - 2 * m2* l1*l2*y[3]*y[4]*cos(y[1] - y[2])) / (2 * l1^2 * l2^2 * (m1 + m2 *sin(y[1] - y[2])^2))^2

    y_[1] = (l2*y[3] - l1*y[4]*cos(y[1] - y[2]))/ (l1^2*l2*(m1 + m2 * sin(y[1] - y[2])^2)) 
    y_[2] = (-m2*l2*y[3]*cos(y[1] - y[2]) +(m1 + m2) *l1*y[4] ) / (m2 * l1 * l2^2 * (m1 + m2*sin(y[1] - y[2])^2))
    y_[3] = -(m1 + m2)*g*l1*sin(y[1])  - h1 + h2 * sin(2*(y[1]-y[2]))
    y_[4] = -m2*g*l2*sin(y[2]) + h1 - h2 * sin(2(y[1] - y[2]))

    return y_
end

function crear_simulacion(;θ1::Float64, θ2::Float64, l1::Float64 = 1.0, l2::Float64 = 1.0, m1::Float64 = 1.0, m2::Float64 = 1.0, g::Float64 = 9.81)
    Simulation([:θ1, :θ2, :p1, :p2], [θ1, θ2, 0.0, 0.0], dyn, 30.0, parameters = Dict([:m1 => m1, :m2 => m2, :l1 => l1, :l2 => l2, :g => g]))
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

function animar_simulacion(sim::Simulation, filename::String = "out.mp4")
    filename = "Tarea4/videosPD/" * filename



    r1 = Observable{Point2f}(r1r2_de_polares(sim.video[1], sim.parameters)[1])
    r2 = Observable{Point2f}(r1r2_de_polares(sim.video[1], sim.parameters)[2])
    cuerdas = @lift [Point2f(0.0, 0.0), $r1, $r2]

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1)
    scatter!(ax, r1, color = :black)
    scatter!(ax, r2, color = :black)
    lines!(ax, cuerdas, color = :black)

    side = sim.parameters[:l1] + sim.parameters[:l2]
    xlims!(ax, (-side * 1.2, side * 1.2))
    ylims!(ax, (-side * 1.2, side * 1.2))

    @show sim.fps
    @show sim.video[end].time
    record(fig, filename, eachindex(sim.video), framerate = Int64(sim.fps)) do i
        current_frame = sim.video[i]
        @show current_frame.time
        r1[], r2[] = r1r2_de_polares(current_frame, sim.parameters)
    end
end