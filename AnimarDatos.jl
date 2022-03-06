using Plots, DataFrames, CSV

function run()
    df = DataFrame(CSV.File("datos_sistema_solar.out"))

    animation = Animation()
    FPS = convert(Int64, round((df.time[2] - df.time[1])^-1))
    println("FPS = " * string(FPS))
    for row in eachrow(df)
        frame_plot = scatter([row.p1x1], [row.p1x2], label="Tierra", xlims=(-5,5), ylim = (-5,5))
        scatter!(frame_plot, [0.0], [0.0], label="Sol")
        frame(animation, frame_plot)
    end
    return gif(animation, fps = FPS)
end
@time run()