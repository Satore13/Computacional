using Random


function generateCoord()
    x = (rand() - 0.5)
    y = (rand() - 0.5)

    return (x,y)
end

function checkIfInCircle((x, y))
    return x^2 + y^2 <=  0.25
end

begin
    iterations = 1000
    local inside::BigInt = 0
    local outside::BigInt = 0
    for i in 1:iterations
        v = generateCoord()
        if(checkIfInCircle(v))
            inside += 1
        else
            outside += 1
        end 
        if ((i*100) ÷ iterations) % 10  == 0
            println((i*100) ÷ iterations)
        end
    end
    
    print("π ≈ ", 4*inside/(inside + outside))
    png(plt, "pipi") 
end