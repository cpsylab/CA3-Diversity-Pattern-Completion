using Base.Threads
using CSV, StableRNGs, Plots, LaTeXStrings, StatsPlots
include("ca3.jl")

# Set parameter ranges for experiments
const ms = [5, 10, 25, 50]
const as = [0.01, 0.05, 0.1, 0.2]
const cs = 0.25:0.25:0.75
const μgs = [0.1, 0.2, 0.3, 0.4]
const σgs = [0., 0.05, 0.1, 0.2]

const rng = Dict(m => StableRNG(326) for m in ms)

# Function to create iterator over parameters 
param_iter(pars) = Iterators.product(pars...) |> collect |> vec 

# Run the simulation for a given set of parameters
println("Starting Simulation")
@threads for m ∈ ms
    println("Running m = $m")
    df = @pipe [
                run_sim(run, NetworkParams(n=3000, nT=10, m=m, a=a, c=c, g1=g1, g1_std=g1_std); rng=rng[m]) 
            for (run, (a, c, g1, g1_std)) in enumerate(param_iter([as, cs, μgs, σgs]))] |> 
        vcat(_...)

    # Write the dataframe to file 
    CSV.write("results/data_m$m.csv", df)
end
