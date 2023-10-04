using CSV, StableRNGs, Plots, LaTeXStrings, StatsPlots
using DataFrames
include("ca3.jl")

################################################################################
# Functions for extracting patterns and calculating sensitivity and specificity
################################################################################
extract_pattern_series(X, pat_index) = mapreduce(permutedims, vcat, [x[pat_index, :] for x ∈ X])
function sensspec(A, B)
    H = (sum((A .== 1) .& (B .== 1), dims=2) ./ sum(A .== 1, dims=2)) |> vec
    FA = (sum((A .== 0) .& (B .== 1), dims=2) ./ sum(A .== 0, dims=2)) |> vec
    return H, FA
end

################################################################################
# Simulation 
################################################################################

serr(x) = std(x) ./ sqrt(length(x))
function hyperhypo_errors(params, rng)
    T, Z = create_patterns(params; rng=rng)
    Z = ceil.(Z)
    WJ = learn_weights(params, T; rng=rng) 
    X, Tprime = init_partial_pattern(params, T, Z; rng=rng) 
    act_corrs,temp_corrs,hit_rates,false_alarm_rates, d_primes, X_series, Tprime_series, g_inhib = run_dynamics(params, Z, T, WJ, X, Tprime; rng=rng, return_patterns=true)
    hyperex = [g <  mean(g_inhib)  for g ∈ g_inhib]
    p_err = mean(Z .!= X_series[end], dims=1) |> vec
    mu_err_h, sd_err_h = @pipe p_err[hyperex .== 1] |> (mean(_), 1.96*serr(_))
    mu_err_l, sd_err_l = @pipe p_err[hyperex .== 0] |> (mean(_), 1.96*serr(_))
    return [mu_err_h, sd_err_h, mu_err_l, sd_err_l]
end


rng = StableRNG(326)
hherrres = []
for g1 ∈ [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    println("Running g1=$g1")
    g1smax = sqrt(g1*(1-g1))
    res_ = hcat([(
        @pipe NetworkParams(n=3000, nT=10, m=10, a=0.1, c=0.5, g1=g1, g1_std=g1std) |> 
            hyperhypo_errors(_, rng) ) for g1std ∈ 0.01:0.02:g1smax ]...)
    push!(hherrres, res_)
end

ps = []
for (i,g1) ∈ enumerate([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    hres = hherrres[i]
    g1smax = sqrt(g1*(1-g1))
    p_ = plot(ylims=[0, 1.0], titlefontsize=11, annotations=(0.3*g1smax, 0.6, text(L"\mu_g="*string(g1), 10)), legend=:topright)
    p_ = plot!(0.01:0.02:g1smax, hres[3,:], ribbon=hres[4,:], label="Hypo", xlabel=L"\sigma_g", ylabel="P(Error)", c=:blue, xrotation=45)
    p_ = plot!(0.01:0.02:g1smax, hres[1,:], ribbon=hres[2,:], label="Hyper", xlabel=L"\sigma_g", ylabel="P(Error)", c=:red, xrotation=45)
    push!(ps, p_)
end
hherrfig = plot(ps...,layout=(2,3), dpi=750)
savefig(hherrfig, "figures/hherrfig.png")