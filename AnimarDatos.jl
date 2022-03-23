using Plots, DataFrames, CSV

function run()
    df = DataFrame(CSV.File("Tarea1/PenduloReal/sim_out.csv"))

    animation = Animation()
    FPS = convert(Int64, round((df.time[2] - df.time[1])^-1))
    println("FPS = " * string(FPS))
    for row in eachrow(df)
        frame_plot = scatter([row.time], [row.p1x1], xlims=(-5,5), ylim = (-5,5), label="")
        frame(animation, frame_plot)
    end
    return gif(animation, fps = FPS)
end
@time run()