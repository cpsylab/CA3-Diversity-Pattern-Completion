using CSV, DataFrames, Plots, StatsPlots, Pipe, Statistics
using LaTeXStrings, Distributions

# Create rainbow color palette
rainbow(n) = palette([:red, :orange, :gold, :green, :blue, :indigo, :violet], n)

# Set parameter ranges for experiments
const ms = 5:5:70
const as = [0.01, 0.05, 0.1, 0.2]
const cs = 0.25:0.25:0.75
const μgs = [0.1, 0.2, 0.25, 0.3, 0.35, 0.4]
const σgs = [0., 0.05, 0.1, 0.2]

ci95(x) = abs(quantile(Normal(0,1), (1-0.95)/2))*std(x)/sqrt(length(x))

# Establish lineplot properties 
lineplot_prop = [

]

df = vcat([@pipe CSV.read("results/data_m$m.csv", DataFrame) |> 
        filter(:timestep => x -> x == 10, _) |> 
        filter(:act_corrs => x -> !isnan(x), _) for m ∈ [5, 10, 15, 20, 25]]...)

meancis(df, var, group) = combine(groupby(df, [group, var]), 
        :act_corrs => mean => :act_corrs_mean, 
        :act_corrs => ci95 => :act_corrs_std,
        :hit_rates => mean => :hit_rates_mean, 
        :hit_rates => ci95 => :hit_rates_std,
        :false_alarm_rates => mean => :false_alarm_rates_mean,
        :false_alarm_rates => ci95 => :false_alarm_rates_std,) 


# Make the main figure 
tickfontsize = 5
p1 = @df meancis(df, :g1_std, :m) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:m, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(5), ylabel="Correlation", title=L"m", legend=:topright, legend_font_pointsize=6, legend_columns=3, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p2 = @df meancis(df, :g1_std, :a) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:a, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(4), title=L"a", legend=:topright, legend_font_pointsize=6, legend_columns=2, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p3 = @df meancis(df, :g1_std, :g1) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:g1, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(6), title=L"g^I",  legend=:topright, legend_font_pointsize=6, legend_columns=2, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p4 = @df meancis(df, :g1_std, :c) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:c, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(3), title=L"c^\ast",  legend=:topright, legend_font_pointsize=6, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)

p5 = @df meancis(df, :g1_std, :m) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:m, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(5), ylabel="Hit Rate", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p6 = @df meancis(df, :g1_std, :a) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:a, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(4), legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p7 = @df meancis(df, :g1_std, :g1) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:g1, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(6),  legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p8 = @df meancis(df, :g1_std, :c) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:c, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(3), legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)

p9 = @df meancis(df, :g1_std, :m) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:m, markerstrokecolor=:auto, ylims=(0,0.75), lw=2, palette=rainbow(5), ylabel="False Alarm Rate", xlabel=L"\sigma_g", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p10 = @df meancis(df, :g1_std, :a) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:a, markerstrokecolor=:auto, ylims=(0,0.75), lw=2, palette=rainbow(4), xlabel=L"\sigma_g", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p11 = @df meancis(df, :g1_std, :g1) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:g1, markerstrokecolor=:auto, ylims=(0,0.75), lw=2, palette=rainbow(6), xlabel=L"\sigma_g",  legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p12 = @df meancis(df, :g1_std, :c) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:c, markerstrokecolor=:auto, ylims=(0,0.75), lw=2, palette=rainbow(3), xlabel=L"\sigma_g", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)


mainfig = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, layout=(3,4), dpi=350)
savefig(mainfig, "results/mainfig.png")