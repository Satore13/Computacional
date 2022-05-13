using OffsetArrays
#= 
    Forma rápida de crear un vector con origen en el 0 a partir del paquete OffsetArrays

=#
Vector0{T} = OffsetVector{T, Vector{T}}
function Vector0(T::Type)
    OffsetVector(T[], OffsetArrays.Origin(0))
end
function Vector0(T::Type, size)
    OffsetVector(Vector{T}(undef, size), OffsetArrays.Origin(0))
end
function Vector0(v::Vector)
    OffsetVector(v, OffsetArrays.Origin(0))
end
struct WaveFunction
    val::Vector0{ComplexF64}
    h::Float64
    WaveFunction(val::Vector0{ComplexF64}) = new(val, 1.0)
end
function (φ::WaveFunction)()
    return φ.val
end
function (φ::WaveFunction)(j::Int64)
    return φ.val[j]
end


function Base.size(φ::WaveFunction)
    return size(φ(), 1)
end

function norm(φ::WaveFunction)::Float64
    sqrt(sum(abs2.(φ()))) 
end
function normalize(φ::WaveFunction)::WaveFunction
    normφ = norm(φ)
    WaveFunction( [ φ(j) / normφ for j in eachindex(φ())] )
end
function probability_density(φ::WaveFunction, j::Int64)::Float64
    abs2.(φ(j))
end

function generar_onda_gaussiana(j0, k0, σ, N)::WaveFunction
    normalize(WaveFunction(Vector0([cis( (j * k0)) * exp(-8*(4*j - N)^2 / N^2) for j in 0:N])))
end

struct QOperator
    f::Function
end
function Base.:+(φ, ψ)
    if size(φ) ≠ size(ψ)
        @error("Las funciones de onda tienen que estar definidas sobre el mismo intervalo en el espacio de posiciones")
    end
    WaveFunction([φ(j) + ψ(j) for j in eachindex(φ())])
end
function (f::QOperator)(φ::WaveFunction)::WaveFunction
    return f.f(φ)
end

function mean_value(φ::WaveFunction, A::QOperator)
    ψ = A(φ)
    real(sum( [ adjoint(phi_n) * psi_n for (phi_n, psi_n) in zip(φ(), ψ())] ))
end

position_op = QOperator(
    (φ) -> WaveFunction([ i*φ(i) for i in eachindex(φ())])
)
momentum_op = QOperator(
    (φ) -> begin
        N = size(φ) - 1
        ψ_val = Vector0(ComplexF64, N + 1)
        ψ_val[N] = 0.0

        for j in 0:(N-1)
            ψ_val[j] = (φ(j+1) - φ(j)) * (-im)
        end
        WaveFunction(ψ_val)
    end
)
function build_potential_op(V::Vector0{ComplexF64})::QOperator
    QOperator((φ) -> WaveFunction([ V[i]*φ(i) for i in eachindex(φ())]))
end
function build_energy_op(V::QOperator)::QOperator
    return QOperator(
        (φ) -> begin
            momentum_op(momentum_op(φ)) + V(φ)
        end
    )
end


struct Simulation
    N::Int64
    n_ciclos::Int64
    s::Float64
    
    φn::Vector0{WaveFunction}
    k0::ComplexF64
    α::Vector0{ComplexF64}
    V::Vector0{ComplexF64}
    A::Vector0{ComplexF64}
end

function calcular_α(N::Int64, A::Vector0{ComplexF64})
    α = Vector0(ComplexF64, N)
    α[N-1] = 0.0
    for j in (N-1):-1:1
        γj = (A[j] + α[j])^-1
        α[j-1] = -γj
    end
    return α
end

function calcular_βn(N::Int64, A::Vector0{ComplexF64}, bn::Vector0{ComplexF64}, α::Vector0{ComplexF64})
    βn = Vector0(ComplexF64, N)
    βn[N-1] = 0.0
    for j in (N-1):-1:1
        γj = (A[j] + α[j])^-1
        βn[j - 1] = γj * (bn[j] - βn[j])
    end
    return βn
end

function calcular_bn(φ::WaveFunction, s::Float64)
    return [4 * im * φj / s for φj in φ()]
end

function calcular_A(V::Vector0{ComplexF64}, s::Float64)
    return [-2 + 2 * im / s - Vj for Vj in V]
end

function calcular_χn(N::Int64, A::Vector0{ComplexF64}, bn::Vector0{ComplexF64}, α::Vector0{ComplexF64}, βn::Vector0{ComplexF64})::Vector0{ComplexF64}
    #Lo primero es calcular las A
    χ = Vector0(ComplexF64, N + 1)
    χ[0] = χ[N] = 0.0
    for j in 1:(N-1)
        χ[j] = - 1/ (A[j] + α[j]) * χ[j - 1] + (bn[j] - βn[j])/(A[j] + α[j])
    end
    return χ
end

function generar_simulacion(N, n_ciclos, x0 = N/4, σ = N/16, choice_of_potential::Symbol = :none, Λ::Vector = [], φi::Union{WaveFunction, Nothing} = nothing)::Simulation
    #Valores independientes del tiempo
    k0 = 2π * n_ciclos / N
    if choice_of_potential == :none
        V = Vector0([Complex(0.0) for _ in 0:N])
    elseif choice_of_potential == :square
        V = Vector0([begin
            if 2 * N / 5 ≤ j ≤ 3*N/5
                Complex(Λ[1]*k0^2)
            else 
                Complex(0.0)
            end
        end for j in 0:N])
        @show V
    end
    
    #Calcular configuraciones iniciales
    if isnothing(φi)
        φ0 = generar_onda_gaussiana(x0, k0, σ, N)
    end
    s = 1/(4*abs2(k0))
    @show s
    A = calcular_A(V, s)
    α = calcular_α(N, A)

    return Simulation(N, n_ciclos, s, Vector0(WaveFunction[φ0]), k0, α, V, A)
end

function paso_simulacion!(simulacion::Simulation)
    current_wf = simulacion.φn[end]
    s = simulacion.s
    N = simulacion.N
    A = simulacion.A
    α = simulacion.α
    current_b = calcular_bn(current_wf, s)
    current_β = calcular_βn(N, A, current_b, α)
    current_χ = calcular_χn(N, A, current_b, α, current_β)
    new_wf = [current_χ[j] - current_wf(j) for j in eachindex(current_wf())]

    push!(simulacion.φn, WaveFunction(new_wf))
end

function bucle_simulacion!(simulation::Simulation, pasos::Integer)
    for i in 1:pasos
        paso_simulacion!(simulation)
    end
end
function bucle_simulacion(N, n_ciclos, pasos; x0 = N/4, σ = N/16, choice_of_potential::Symbol = :none, Λ::Vector = [])
    sim = generar_simulacion(N, n_ciclos, x0, σ, choice_of_potential, Λ)
    bucle_simulacion!(sim, pasos)
    return sim
end