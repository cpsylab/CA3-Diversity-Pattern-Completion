using CSV, StableRNGs, Plots, LaTeXStrings, StatsPlots
include("ca3.jl")

const rng = StableRNG(326)

# Set parameter ranges for experiments
const ms = 5:5:70
const as = [0.01, 0.05, 0.1, 0.2]
const cs = 0.25:0.25:0.75
const μgs = [0.1, 0.2, 0.25, 0.3, 0.35, 0.4]
const σgs = [0., 0.05, 0.1, 0.2]

# Function to create iterator over parameters 
param_iter(ms, μgs, σgs) = Iterators.product(ms, μgs, σgs) |> collect |> vec 

# Run the simulation for a given set of parameters
for m ∈ ms
    println("Running m = $m")
    df = @pipe [
            run_sim(run, NetworkParams(n=1000, m=m, a=a, c=c, g1=g1, g1_std=g1_std); rng=rng) 
            for (run, (a, c, g1, g1_std)) in enumerate(param_iter(as, cs, μgs, σgs))] |> 
        vcat(_...)

    # Write the dataframe to file 
    CSV.write("results/data_m$m.csv", df)
end

# df_act = combine(groupby(df, [:run, :m, :g1, :g1_std, :timestep]), 
#     :act_corrs => mean => :act_corrs, 
#     :act_corrs => std => :act_corrs_std, 
#     :hit_rates => mean => :hit_rates,
#     :hit_rates => std => :hit_rates_std, 
#     :false_alarm_rates => mean => :false_alarm_rates, 
#     :false_alarm_rates => std => :false_alarm_rates_std)

# figs = [
#     @df filter([:m, :g1] => (x,y) -> (x == m) & (y == g1), df_act) plot(
#         :timestep, :act_corrs, yerror=:act_corrs_std, group=:g1_std, ylim=[0,1],
#         title="$m, $g1", xlabel="Timestep", ylabel="Correlation",
#         markerstrokecolor=:auto, legend=false, titlefont = font(5), 
#         xtickfontsize=5, ytickfontsize=5, xlabelfontsize=6, ylabelfontsize=6,
#     ) for (g1, m) in (Iterators.product(μgs, ms) |> collect |> vec)
# ] 

# plot(figs..., layout=(4, 6))

# plot(df_act.timestep, df_act.act_corrs, yerror=df_act.act_corrs_std, legend=false, c=:black, ylim=[0,1])
# plot!(df_act.timestep, df_act.hit_rates, yerror=df_act.hit_rates_std, legend=false, c=:blue, ylim=[0,1], linecolor=:blue)
# plot!(df_act.timestep, df_act.false_alarm_rates, yerror=df_act.false_alarm_rates_std, ylim=[0,1], c=:red, legend=false, linecolor=:red)

# plot(act_corrs, legend=false, c=:black, alpha=0.5, ylim=[0,1])
# plot!(hit_rates, legend=false, c=:blue, alpha=0.5, ylim=[0,1])
# plot!(false_alarm_rates, ylim=[0,1], c=:red, alpha=0.5)
