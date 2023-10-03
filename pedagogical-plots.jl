using Plots, LaTeXStrings
using Distributions

# Plot the relationship between T_ij and Z_ij 
zfig = plot(-2:0.01:2, x->exp(-x^2), 
    c=:black, legend=false, 
    xlabel=L"T_{ij}", ylabel=L"Z_{ij} = e^{-T_{ij}^2}", 
    yticks=0:0.2:1, lw=2,
    size=(200, 200))
savefig(zfig, "figures/TZ.pdf")

# stdp plot 
stdp_fig = plot(
    -2:0.01:2, x->exp(-abs(x)/1), 
    c=:black, legend=false, lw=2,
    xlabel=L"T_{ki} - T_{kj}", ylabel=L" e^{-|T_{ki} - T_{kj}|/\tau_{pot}}", 
    size=(200, 200))
savefig(stdp_fig, "figures/stdp.pdf")

# Beta distribution plot
"""
    BetaMS(m::Float64, s::Float64)

    Create a Beta distribution with mean `m` and standard deviation `s`.
"""
function BetaMS(m::Float64, s::Float64)
    @assert s < sqrt(m*(1-m))

    α = (m^2 - m^3 - m*s^2)/ s^2
    β = ((m-1)*(m^2 - m + s^2))/s^2
    return Beta(α, β)
end

colors = [:red, :orange, :green, :blue, :purple]
fig = plot() 
for s ∈ [0.1, 0.2, 0.3, 0.4]
    fig = plot!(0:0.01:1, x -> pdf(BetaMS(0.4, s), x), 
        label=L"\sigma_g="*"$s", lw=2,
        c=colors[findfirst(x->x==s, [0.1, 0.2, 0.3, 0.4])],
        xlabel=L"g_i^I", ylabel="PDF", 
        xticks=0:0.2:1, 
        size=(200, 200))
end
fig
savefig(fig, "figures/beta.pdf")

plot(fig, stdp_fig, layout=(2,1), size=(200, 400))
