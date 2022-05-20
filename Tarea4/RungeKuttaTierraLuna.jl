include("RungeKutta.jl")
using GLMakie

const ω = 2.6617e-6
const masa_nave = 1.0
const masa_tierra = 5.9736e24
const masa_luna = 0.07349e24
const distancia_tierra_luna = 3.844e8
const constante_gravitacional = 6.67e-11
const radio_tierra = 6.378160e6
const radio_luna = 1.7374e6

function dyn(y, t)
    dy = copy(y)

    Δ = constante_gravitacional * masa_tierra / distancia_tierra_luna^3
    μ = masa_luna / masa_tierra     
    r_aux = sqrt(1 + y[1] ^2 - 2*y[1]*cos(y[2] - ω * t))

    dy[1] = y[3]
    dy[2] = y[4]/y[1]^2
    dy[3] = y[4]^2/y[1]^3 - Δ * (1 / y[1] ^ 2 + μ / r_aux^3 * (y[1]  -  cos(y[2] - ω * t)))
    dy[4] = -Δ * μ * y[1] / r_aux^3 * sin(y[2] - ω * t) 


    return dy

end


#Velocidad de escape ≈ 11_190 m /s
#La velocidad inicial se espera en m/s y los ángulos en grados
function run(;ρ0::Float64 = radio_tierra, v0::Float64, θ::Float64, φ::Float64, length::Float64, time_scale::Float64, h::Float64 = 1.0, fps = 30.0)
    v_norm = v0 / distancia_tierra_luna

    #Momento lineal inicial en la coordenada polar (en coordenadas normalizadas)
    p_ρ0 = v_norm * cosd(θ - φ)

    ρ_norm = ρ0 / distancia_tierra_luna

    #Momento angular inicial en la coordenada azimutal (en coordenadas normalizadas)
    p_φ0 = ρ_norm * v_norm * sind(θ - φ)

    y0 = [ρ0 / distancia_tierra_luna, deg2rad(φ) , p_ρ0, p_φ0]
    s = Simulation([:ρ, :φ, :pρ, :pφ], y0, dyn, h, fps, time_scale = time_scale)
    loop!(s, length)
    return s
end

#Algunos ejemplos
run_orbita_geosincrona() = run(ρ0 = 35_786_000.0+radio_tierra, v0 = 3_070.0, φ = 90.0, θ = 0.0, length = 48.0, time_scale = 3600.0, h = 1.0)


function animar_simulacion(sim::Simulation, filename::String = "out.mp4")
    filename = "Tarea4/videos/"*filename

    fig = Figure()
    time_scale = sim.time_scale
    h = floor(Int,  time_scale / 3600)
    m = floor(Int, (time_scale % 3600) / 60) 
    s = time_scale % 60


    current_time_str = Observable{String}("")
    ax = Axis(fig[1, 1], aspect = 1, title  = "1s animación = $(h)h $(m)m $(s)s real", xlabel = current_time_str)
    xlims!(ax, -1.1, 1.1)
    ylims!(ax, -1.1, 1.1)


    current_pos_luna = Observable{Point2f}((0,0))
    current_pos = Observable{Point2f}((0,0))
    scatter!(ax, current_pos)
    scatter!(ax, current_pos_luna)

    record(fig, filename, eachindex(sim.video), framerate = Int64(sim.fps)) do i
        current_frame = sim.video[i]

        time = current_frame.time
        
        h = floor(Int,  time / 3600)
        m = floor(Int, (time % 3600) / 60) 
        s = time_scale % 60

        current_pos_luna[] = (cos(ω * time), sin(ω * time))
        current_time_str[] = "Δt real = $(h)h $(m)m $(s)s"
        ρ = current_frame.y[1]
        φ = current_frame.y[2]
        @show ρ
        @show φ
        @show current_frame.time
        current_pos[] = (ρ * cos(φ), ρ * sin(φ))
    end
end