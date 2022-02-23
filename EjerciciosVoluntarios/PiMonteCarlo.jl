using Random


function generateCoord()
    x = (rand())
    y = (rand())

    return (x,y)
end

function checkIfInCircle((x, y))
    return x^2 + y^2 <=  1
end

begin
    N = 1000
    local inside::BigInt = 0
    for i in 1:N
        v = generateCoord()
        if(checkIfInCircle(v))
            inside += 1
        end
    end
    
    println("π ≈ ", 4*inside/(N))
end