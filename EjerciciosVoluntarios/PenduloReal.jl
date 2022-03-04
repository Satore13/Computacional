include("Verlet.jl")
Particle = Verlet.Particle
Frame = Verlet.Frame
begin
    N = 1
    global f = Frame(0.0, [Particle([0.0], [float(2*π)], 1.0)])
    simulation_data = Verlet.dataFrameRowFromFrame(f)
    for i in 1:10000
        global f = Verlet.stepFrame(f, 0.001, Verlet.calculateForceOnParticle)
        Verlet.addFrameRowToDataFrame!(simulation_data, f)
    end
    show(simulation_data)
end 