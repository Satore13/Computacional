import Base.size
import Base.rand
import Base.show
struct Spin
    val::Integer
    function Spin(val::Integer)
        if !(val in [-1, 1])
            return error("El spin sólo puede valer ±1")
        else 
            return new(val)
        end
    end
end
function rand(::Type{Spin})
    return Spin(rand([-1, 1]))
end
function inverse_spin(s::Spin)
    if s === Spin(1)
        return Spin(-1)
    else 
        return Spin(1)
    end
end

struct Red
    Nudos::Matrix{Spin}
    Temperatura::Float64
    Red(Nudos, Temp) = size(Nudos)[1] != size(Nudos)[2] ? error("La red tiene que ser cuadrada") : new(Nudos, Temp)
    Red(Nudos) = Red(Nudos, 0.0)
    Red(s::Spin, N::Integer) = Red([s for _ in 1:N, _ in 1:N])
end
function rand(::Type{Red}, N::Integer)
    s = [rand(Spin) for _ in 1:N, _ in 1:N]
    return Red(s)
end
function rand(::Type{Red}, N::Integer, temp::Float64)
    s = [rand(Spin) for _ in 1:N, _ in 1:N]
    return Red(s, temp)
end
function size(r::Red)
    return size(r.Nudos)[1]
end
#Generar una red igual salvo inversión de un spin aleatorio
function generar_red_similar(r::Red)
    N = size(r)
    i, j = rand(1:N, 2)

    new_spins = deepcopy(r.Nudos)
    new_spins[i, j] = inverse_spin(new_spins[i ,j])

    return Red(new_spins, r.Temperatura)
end

function energia(r::Red)

    #Conseguimos la matriz de spines y el tamaño de la red
    s = getfield.(r.Nudos, :val)
    N = size(r)

    #Hacemos equivalentes los spines 0 y N - Condiciones periódicas
    function wrap(i)
        if i > N
            i -= N
        elseif i < 1
            i += N
        end
        return i
    end
    
    #Elementos de la sumatoria
    elms = [ 
        begin
            s[i, j] * (s[i, wrap(j + 1)] + s[i, wrap(j - 1)] + s[wrap(i + 1), j] + s[wrap(i - 1), j])
        end
        for i = 1:N, j in 1:N]


    return -0.5 * sum(elms) / (2*N)
end

function magnetizacion(r::Red)
    s = getfield.(r.Nudos, :val)
    N = size(r)

    return abs(sum(sij for sij in s)/N^2)
end

function probabilidad_transicion(from::Red, to::Red)
    if !(from.Temperatura ≈ to.Temperatura)
        error("La temperatura ha de ser la misma para ambas redes")
    end
    Δenergia = energia(to) - energia(from)

    β = 1/from.Temperatura

    return min(1, exp(-β * Δenergia))
end

#= Paso del algoritmo de Metropolis:
        Generamos una red similar y calculamos la probabilidad de transición, si un número generado aleatoria y uniformemente es menor que esa probabilidad 
    pasamos a esa nueva red. En caso contrario, nos quedamos con la configuración previa
=#
function paso(r::Red)
    #Generamos una red igual salvo un spin invertido
    nueva_r = generar_red_similar(r)

    p = probabilidad_transicion(r, nueva_r)

    ξ = rand()
 
    if ξ < p
        return nueva_r
    else 
        return r
    end
end

function bucle_simulacion(red_inicial::Red, pasos_de_MC::Integer, guardar_cada::Integer = 1)::Vector{Red}
    N = size(red_inicial)
    r = red_inicial
    datos = Red[]
    push!(datos, red_inicial)

    porcentaje_previo = 0.0

    for n_pMC in 1:pasos_de_MC
        for _ in 1:N^2 
            r = paso(r)
        end
        porcentaje_actual = n_pMC / pasos_de_MC * 100
        if porcentaje_actual - porcentaje_previo > 0.01
            print("Progreso: ", round(porcentaje_actual, digits = 2), "%\r")
            porcentaje_previo = porcentaje_actual
        end
        if n_pMC % guardar_cada == 0
            push!(datos, r)
        end
    end

    return datos
end

