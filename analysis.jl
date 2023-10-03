using CSV, DataFrames, Plots, StatsPlots, Pipe, Statistics
using LaTeXStrings, Distributions, GLM

# Create rainbow color palette
rainbow(n) = palette([:red, :orange, :gold, :green, :blue, :indigo, :violet], n)

# Set parameter ranges for experiments
const ms = 5:5:70
const as = [0.01, 0.05, 0.1, 0.2]
const cs = 0.25:0.25:0.75
const μgs = [0.1, 0.2, 0.25, 0.3, 0.35, 0.4]
const σgs = [0., 0.05, 0.1, 0.2]

ci95(x) = abs(quantile(Normal(0,1), (1-0.95)/2))*std(x)/sqrt(length(x))
compute_sparsity(df) = 1 .- ((1 .- df.a) .* df.false_alarm_rates + df.a .* df.hit_rates)


# LINEAR REGRESSION 
# df = vcat([@pipe CSV.read("results/data_m$m.csv", DataFrame) |> 
#         filter(:timestep => x -> x == 10, _) |> 
#         filter(:act_corrs => x -> !isnan(x), _) for m ∈ 5:5:65]...)
# df.sparsity = compute_sparsity(df)

# res = lm(@formula(act_corrs ~ m + a + c + g1 + g1_std), df)

# lm(@formula(act_corrs ~ m + a*g1*g1_std), df)


# PLOTS 
df = vcat([@pipe CSV.read("results/data_m$m.csv", DataFrame) |> 
        filter(:timestep => x -> x == 10, _) |> 
        filter(:act_corrs => x -> !isnan(x), _) for m ∈ [5,15, 25, 50, 60]]...)

df.sparsity = compute_sparsity(df)


meancis(df, var, group) = combine(groupby(df, [group, var]), 
        :act_corrs => mean => :act_corrs_mean, 
        :act_corrs => ci95 => :act_corrs_std,
        :hit_rates => mean => :hit_rates_mean, 
        :hit_rates => ci95 => :hit_rates_std,
        :false_alarm_rates => mean => :false_alarm_rates_mean,
        :false_alarm_rates => ci95 => :false_alarm_rates_std,
        :sparsity => mean => :sparsity_mean,
        :sparsity => ci95 => :sparsity_std) 


# Make the main figure 
Nm = 5
tickfontsize = 5; ylabelfontsize=8
p1 = @df meancis(df, :g1_std, :m) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:m, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(Nm), ylabel="Correlation", title=L"m", legend=:topright, legend_font_pointsize=6, legend_columns=3, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize, ylabelfontsize=ylabelfontsize)
p2 = @df meancis(df, :g1_std, :a) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:a, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(4), title=L"a", legend=:topright, legend_font_pointsize=6, legend_columns=2, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p3 = @df meancis(df, :g1_std, :g1) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:g1, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(6), title=L"\mu_g",  legend=:topright, legend_font_pointsize=6, legend_columns=2, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p4 = @df meancis(df, :g1_std, :c) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, group=:c, markerstrokecolor=:auto, ylims=(0,1.25), lw=2, palette=rainbow(3), title=L"c^\ast",  legend=:topright, legend_font_pointsize=6, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)

p5 = @df meancis(df, :g1_std, :m) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:m, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(Nm), ylabel="Hit Rate", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize, ylabelfontsize=ylabelfontsize)
p6 = @df meancis(df, :g1_std, :a) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:a, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(4), legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p7 = @df meancis(df, :g1_std, :g1) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:g1, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(6),  legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p8 = @df meancis(df, :g1_std, :c) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, group=:c, markerstrokecolor=:auto, ylims=(0,1.1), lw=2, palette=rainbow(3), legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)

p9 = @df meancis(df, :g1_std, :m) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:m, markerstrokecolor=:auto, ylims=(0,1), lw=2, palette=rainbow(Nm), ylabel="False Alarm Rate", xlabel=L"\sigma_g", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize, ylabelfontsize=ylabelfontsize)
p10 = @df meancis(df, :g1_std, :a) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:a, markerstrokecolor=:auto, ylims=(0,1), lw=2, palette=rainbow(4), xlabel=L"\sigma_g", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p11 = @df meancis(df, :g1_std, :g1) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:g1, markerstrokecolor=:auto, ylims=(0,1), lw=2, palette=rainbow(6), xlabel=L"\sigma_g",  legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)
p12 = @df meancis(df, :g1_std, :c) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:c, markerstrokecolor=:auto, ylims=(0,1), lw=2, palette=rainbow(3), xlabel=L"\sigma_g", legend=false, ytickfontsize=tickfontsize, xtickfontsize=tickfontsize)


mainfig = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, layout=(3,4), dpi=500)
savefig(mainfig, "results/mainfig.png")

# Activity correlation crossfig
plts = []
for (a, μg) ∈ Iterators.product(as, μgs)
        plt = @pipe df |> 
        filter(:a => x -> x == a, _) |>
        filter(:g1 => x -> x == μg, _) |> 
        meancis(_, :g1_std, :m) |> (
        @df sort!(_, :g1_std) plot(:g1_std, :act_corrs_mean, yerr=:act_corrs_std, 
        group=:m, markerstrokecolor=:auto, 
        ylims=(0,1), yticks=0:1,
        xticks=0:0.1:0.2,
        legend=:topright, 
        size=(200,200), xlabel=L"\sigma_g", ylabel="ρ", 
        title="a=$a, μg=$μg", titlefontsize=8,
        xlabelfontsize=8, ylabelfontsize=8, tickfontsize=8,
        palette=rainbow(length(ms))))
        push!(plts, plt)
end
crossfig = plot(plts..., layout=(6,4), size=(900,1200), legend=false, dpi=500)
savefig(crossfig, "figures/crossfig.png")

# Hit rate correlation crossfig
plts = []
for (a, μg) ∈ Iterators.product(as, μgs)
        plt = @pipe df |> 
        filter(:a => x -> x == a, _) |>
        filter(:g1 => x -> x == μg, _) |> 
        meancis(_, :g1_std, :m) |> (
        @df sort!(_, :g1_std) plot(:g1_std, :hit_rates_mean, yerr=:hit_rates_std, 
        group=:m, markerstrokecolor=:auto, 
        ylims=(0,1), yticks=0:0.2:1,
        xticks=0:0.1:0.2,
        legend=:topright, 
        size=(200,200), xlabel=L"\sigma_g", ylabel="H", 
        title="a=$a, μg=$μg", titlefontsize=8,
        xlabelfontsize=8, ylabelfontsize=8, tickfontsize=8,
        palette=rainbow(length(ms))))
        push!(plts, plt)
end
crossfig_hits = plot(plts..., layout=(6,4), size=(900,1200), legend=false, dpi=500)
savefig(crossfig_hits, "figures/crossfig_hits.png")

# FA rate correlation crossfig
plts = []
for (a, μg) ∈ Iterators.product(as, μgs)
        plt = @pipe df |> 
        filter(:a => x -> x == a, _) |>
        filter(:g1 => x -> x == μg, _) |> 
        meancis(_, :g1_std, :m) |> (
        @df sort!(_, :g1_std) plot(:g1_std, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, 
        group=:m, markerstrokecolor=:auto, 
        ylims=(0,1), yticks=0:0.2:1,
        xticks=0:0.1:0.2,
        legend=:topright, 
        size=(200,200), xlabel=L"\sigma_g", ylabel="FA", 
        title="a=$a, μg=$μg", titlefontsize=8,
        xlabelfontsize=8, ylabelfontsize=8, tickfontsize=8,
        palette=rainbow(length(ms))))
        push!(plts, plt)
end
crossfig_fa = plot(plts..., layout=(6,4), size=(900,1200), legend=false, dpi=500)
savefig(crossfig_fa, "figures/crossfig_fa.png")

# FA rate correlation crossfig
plts = []
for (a, μg) ∈ Iterators.product(as, μgs)
        plt = @pipe df |> 
        filter(:a => x -> x == a, _) |>
        filter(:g1 => x -> x == μg, _) |> 
        meancis(_, :g1_std, :m) |> (
        @df sort!(_, :g1_std) plot(:g1_std, :sparsity_mean, yerr=:sparsity_std, 
        group=:m, markerstrokecolor=:auto, 
        ylims=(0,1), yticks=0:0.2:1,
        xticks=0:0.1:0.2,
        legend=:topright, 
        size=(200,200), xlabel=L"\sigma_g", ylabel="FA", 
        title="a=$a, μg=$μg", titlefontsize=8,
        xlabelfontsize=8, ylabelfontsize=8, tickfontsize=8,
        palette=rainbow(length(ms))))
        push!(plts, plt)
end
crossfig_sp = plot(plts..., layout=(6,4), size=(900,1200), legend=false, dpi=500)
savefig(crossfig_p, "figures/crossfig_p.png")

# PLOT DYNAMICS OVER TIME 
meancis_t(df, group) = combine(groupby(df, [group, :timestep]), 
        :act_corrs => mean => :act_corrs_mean, 
        :act_corrs => ci95 => :act_corrs_std,
        :hit_rates => mean => :hit_rates_mean, 
        :hit_rates => ci95 => :hit_rates_std,
        :false_alarm_rates => mean => :false_alarm_rates_mean,
        :false_alarm_rates => ci95 => :false_alarm_rates_std,
        :sparsity => mean => :sparsity_mean,
        :sparsity => ci95 => :sparsity_std) 



df_t = vcat([@pipe CSV.read("results/data_m$m.csv", DataFrame) |> 
        filter(:act_corrs => x -> !isnan(x), _) for m ∈ [5, 10, 15, 20, 25, 30, 35, 45, 60]]...)
df_t.sparsity = compute_sparsity(df_t)



val = Symbol(:g1_std)
@df meancis_t(df_t, :m) plot(:timestep, :sparsity_mean, yerr=:sparsity_std, group=:m, linewidth=2, markerstrokecolor=:auto, palette=rainbow(length(unique(:m))), ylims=(0,1))
@df meancis_t(df_t, :g1_std) plot(:timestep, :hit_rates_mean, yerr=:hit_rates_std, group=:g1_std, linewidth=2, markerstrokecolor=:auto, palette=rainbow(length(unique(:g1_std))), ylims=(0,1), legend=:bottomright)
@df meancis_t(df_t, :g1_std) plot!(:timestep, :false_alarm_rates_mean, yerr=:false_alarm_rates_std, group=:g1_std, linewidth=2, linestyle=:dash, markerstrokecolor=:auto, palette=rainbow(length(unique(:g1_std))), ylims=(0,1))


@df meancis_t(df_t, :m) scatter(:sparsity_mean, :act_corrs_mean, group=:m, linewidth=2, markerstrokecolor=:auto, palette=rainbow(length(unique(:m))), xlims=(0.5,1), ylims=(0,1), legend=:bottomright)