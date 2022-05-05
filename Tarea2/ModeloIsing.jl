struct Spin
    val::Int64
    function Spin(val::Int64)
        if !(val in [-1, 1])
            return error("El spin sólo puede valer ±1")
        else 
            return new(val)
        end
    end
end
function Base.rand(::Type{Spin})
    return Spin(rand([-1, 1]))
end
function inverse_spin(s::Spin)::Spin
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
    Red(s::Spin, N::Int64, T::Float64) = Red([s for _ in 1:N, _ in 1:N], T)
    Red(s::Spin, N::Int64) = Red(s, N, 0.0)
end
function Base.rand(::Type{Red}, N::Integer)
    s = [rand(Spin) for _ in 1:N, _ in 1:N]
    return Red(s)
end
function Base.rand(::Type{Red}, N::Integer, temp::Float64)
    s = [rand(Spin) for _ in 1:N, _ in 1:N]
    return Red(s, temp)
end
function Base.size(r::Red)::Int64
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

function energia_norm(r::Red)::Float64

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


    return -0.5 * sum(elms) / size(r)^2 
end

function Δenergia(from::Red, to::Red)
    #Conseguimos la matriz de spines y el tamaño de la red
    to_s = getfield.(to.Nudos, :val)
    from_s = getfield.(from.Nudos, :val)

    N = size(from)

    #Hacemos equivalentes los spines 0 y N - Condiciones periódicas
    function wrap(i)
        if i > N
            i -= N
        elseif i < 1
            i += N
        end
        return i
    end
    
    function elm(s, i, j)
        s[i, j] * (s[i, wrap(j + 1)] + s[i, wrap(j - 1)] + s[wrap(i + 1), j] + s[wrap(i - 1), j])
    end
    #Elementos de la sumatoria
    elms = [ 
        begin
            elm(to_s, i, j) - elm(from_s, i, j) 
        end
        for i = 1:N, j in 1:N]


    return -0.5 * sum(elms)
end

function magnetizacion_norm(r::Red)
    s = getfield.(r.Nudos, :val)
    return sum(sij for sij in s)/ size(r)^2
end

function probabilidad_transicion(from::Red, to::Red)
    if !(from.Temperatura ≈ to.Temperatura)
        error("La temperatura ha de ser la misma para ambas redes")
    end

    β = 1/from.Temperatura

    return min(1, exp(-β * Δenergia(from, to)))
end

#= Paso del algoritmo de Metropolis:
        Generamos una red similar y calculamos la probabilidad de transición, si un número generado aleatoria y uniformemente es menor que esa probabilidad 
    pasamos a esa nueva red. En caso contrario, nos quedamos con la configuración previa
=#

function energia_entre_similares(r::Red, x::Int64, y::Int64)::Float64
    #Forma rápida de calcular la diferencia de energía con una configuración similar
    N = size(r)
    function wrap(i::Int64)::Int64
        if i > N
            i -= N
        elseif i < 1
            i += N
        end
        return i
    end
    s0, s1, s2, s3, s4 = r.Nudos[x, y].val, r.Nudos[wrap(x - 1), y].val, r.Nudos[wrap(x+1), y].val, r.Nudos[x, wrap(y+1)].val, r.Nudos[x, wrap(y-1)].val

    return float(2*s0*(s1 + s2 + s3 + s4))

    #=return float(s[x, y] * (s[wrap(x - 1), y] + s[wrap(x + 1), y] + s[x, wrap(y + 1)] + s[x, wrap(y - 1)]) +
            s[x, y] * s[wrap(x - 1), y] +
            s[x, y] * s[wrap(x + 1), y] +
            s[x, y] * s[x, wrap(y - 1)] +
            s[x, y] * s[x, wrap(y + 1)])=#
end

function paso!(r::Red, x::Int64, y::Int64)
    β = 1/r.Temperatura::Float64
    ΔE = energia_entre_similares(r, x, y)
    p = exp(-β * ΔE)::Float64
    if rand() <= p
        r.Nudos[x, y] = inverse_spin(r.Nudos[x, y]) 
    end
    return
end

function pasoMC!(r::Red)
    N = size(r)::Int64
    for _ in 1:N
        directions = Tuple{Int64, Int64}[]
        
        for _ in 1:N
            notok = true
            while notok
                new_d = (rand(1:N), rand(1:N))
                if !(in(directions)(new_d))
                    push!(directions, new_d)
                    notok = false
                end
            end
        end
        Threads.@threads for (x, y) in directions
            paso!(r, x, y)
        end
    end
end

function bucle_simulacion(red_inicial::Red, pasos_de_MC::Int64, guardar_cada::Int64 = 1)::Vector{Red}
    N = size(red_inicial)
    r = deepcopy(red_inicial)
    datos = Red[]
    push!(datos, red_inicial)
    porcentaje_previo = 0.0
    for n_pMC in 1:pasos_de_MC
        pasoMC!(r)
        porcentaje_actual = n_pMC / pasos_de_MC * 100
        if porcentaje_actual - porcentaje_previo > 0.1
            print("Progreso: ", round(porcentaje_actual, digits = 2), "%\r")
            porcentaje_previo = porcentaje_actual
        end
        if n_pMC % guardar_cada == 0
            push!(datos, deepcopy(r))
        end
    end

    return datos
end

function local_correlation(datos::Vector{Red}, x::Integer, y::Integer, r::Integer)
    length_of_sim = size(datos, 1)
    traslations = [(r, 0), (-r, 0), (0, r), (0, -r)]
    N = size(datos[1])

    #Hacemos equivalentes los spines 0 y N - Condiciones periódicas
    function wrap(i)
        if i > N
            i -= N
        elseif i < 1
            i += N
        end
        return i
    end
    
    suma = 0.0
    for r in datos
        my_spin = r.Nudos[x, y].val
        other_positions = [(wrap(x + tx), wrap(y + ty)) for (tx, ty) in traslations]
        suma += sum([ my_spin * r.Nudos[other_x, other_y].val for (other_x, other_y) in other_positions]) / 4
    end
    return suma / length_of_sim
end

function global_correlation(red::Vector{Red}, r::Integer)
    N = size(red[1])
    suma = 0.0
    for i in 1:N, j in 1:N
        suma += local_correlation(red, i, j, r)
    end

    return suma / N^2
end

function calcular_correlacion(datos::Vector{Red})::Vector{NTuple{2, Float64}}
    max_val = ceil(Int64, size(datos[1])/2)
    r = collect(0:max_val)
    return [(ri, global_correlation(datos, ri)) for ri in r]
end

function calor_especifico(datos::Vector{Red})::Float64
    number_of_states = size(datos, 1)
    T = datos[1].Temperatura
    mean_energy = 0.0
    mean_energysq = 0.0
    for r in datos
        mean_energy += energia_norm(r)
        mean_energysq += energia_norm(r)
    end
    mean_energy_and_then_sq  = (mean_energy / number_of_states)^2
    mean_energysq = mean_energysq / number_of_states^2
    return (1/T) * (mean_energysq - mean_energy_and_then_sq)
end