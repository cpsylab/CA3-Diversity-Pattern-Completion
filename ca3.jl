################################################################################
# ca3.jl 
#     Contains functions defining the operation of the CA3 network 
# 
#    Functions:
#       PART 1: Utility functions
#           BetaMS(m, s)
#           corrcoef(x, y)
#       PART 2: Encoding phase functions 
#           create_patterns(params; rng)
#           learn_weights(params, T; rng=rng)
#       PART 3: Retrieval phase functions
#           init_partial_pattern(params, T, Z; rng)
#           run_dynamics(params, Z, T, WJ, X, Tprime; rng)
#
# Written by Abraham Nunes, 2023 
# Extended from Matlab and C code from 
#   [1] Mishra RK, Kim S, Guzman SJ, Jonas P. Symmetric spike timing-dependent 
#       plasticity at CA3–CA3 synapses optimizes storage and recall in 
#       autoassociative networks. Nat Commun. 2016;7(1):11552. 
#   [2] Guzman SJ, Schlögl A, Frotscher M, Jonas P. Synaptic mechanisms of 
#       pattern completion in the hippocampal CA3 network. Science. 
#       2016;353(6304):1117-1123.
#   [3] Bennett MR, Gibson WG, Robinson J. Dynamics of the CA3 Pyramidal Neuron 
#       Autoassociative Memory Network in the Hippocampus.
#   [4] Gibson WG, Robinson J. Statistical analysis of the dynamics of a 
#       sparse associative memory. Neural Networks. 1992;5(4):645-661.
################################################################################
using Parameters, DataFrames
using Tullio, Pipe
using Distributions, Random, Statistics

@with_kw struct NetworkParams 
    m::Int = 10         # n patterns stored 
    n::Int = 1000       # n neurons 
    c::Float64 = 0.5    # prob connection between neurons
    a::Float64 = 0.1    # prob activation of each neuron
    b1::Float64 = 0.5   # Fraction of correct active pattern 
    bn::Float64 = 1.    # Fraction of spurious active pattern 
    σt::Float64 = 0.2   # Standard deviation of spike times
    g0::Float64 = 0.    # base threshold value
    g1::Float64 = 0.2   # inhibitory neuron gain (mean in the case of "beta" inhibition)
    g1_std::Float64 = 0. # std of inhibitory neuron gain (must be < sqrt(mean*(1-mean)))
    τpot::Float64 = 1.  # STDP time constant 
    nT::Int = 5         # number of time steps 
    tau3::Float64 = 1.  # Neuronal membrane potential time constant 
    stdp_type::String = "symmetric" # "symmetric" or "asymmetric
    connectivity_algorithm:: String = "random" # "random" or "sonet" 
    αrecip::Float64=1. 
    αconv::Float64=1. 
    αdiv::Float64=1. 
    αchain::Float64=1.
end

# -----------------------------------------------------------------------------
# PART 1: Utility functions 
# -----------------------------------------------------------------------------

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

"""
    corrcoef(x::Vector, y::Vector)

    Calculate the correlation coefficient between finite elements of two vectors.
"""
function corrcoef(x::Vector, y::Vector)
    finite_indices = intersect(findall(isfinite, x), findall(isfinite, y))
    x, y = x[finite_indices], y[finite_indices]
    return mean((x .- mean(x)) .* (y .- mean(y))) / (std(x) * std(y))
end

"""
    sdt_analysis(x::Vector, y::Vector)

    Compute the hit rate, false alarm rate, and d' for two vectors `x` and `y`.
"""
norminv(x) = quantile(Normal(0, 1), x)
function sdt_analysis(x::Vector, y::Vector)
    # Find indices where x and y are nonzero, and indices where y is zero
    x_nonzero, y_nonzero = findall(!iszero, x), findall(!iszero, y)
    y_zero = findall(iszero, y)

    # Compute hit rate and false alarm rate
    hit_rate = length(intersect(x_nonzero, y_nonzero)) / length(y_nonzero)
    false_alarm_rate = length(intersect(x_nonzero, y_zero)) / length(y_zero)

    # Compute d'
    d_prime = norminv(hit_rate) - norminv(false_alarm_rate)
    return [hit_rate, false_alarm_rate, d_prime]
end

# -----------------------------------------------------------------------------
# PART 2: Functions for encoding phase 
# -----------------------------------------------------------------------------

""" 
    create_patterns(params::NetworkParams; rng=rng)

    Create patterns for the network.
"""
function create_patterns(params::NetworkParams; rng=rng)
    @unpack n,m,a,σt = params
    idz = randperm(rng, n*m)[1:round(Int, m*n*a)]   # indices of non-zero elements
    T, Z = fill(Inf, m, n), zeros(m, n)             # initialize 
    T[idz] .= rand(rng, Normal(0, σt), length(idz)) # fill non-zero elements with random values
    Z[idz] .= exp.(-(T[idz].^2))
    return T,Z
end

"""
    make_sonet_connectivity_matrix(params::NetworkParams; rng=rng)

    Create a connectivity matrix for the network using the SONET algorithm.

    Reference: 
        Zhao L, Beverlin B, Netoff T, Nykamp D. Synchronization from Second 
        Order Network Connectivity Statistics. Frontiers in Computational 
        Neuroscience. 2011;5.
"""
function make_sonet_connectivity_matrix(params::NetworkParams; rng=rng)
    return 0
end 


""" 
    make_connectivity_matrix(params::NetworkParams; rng=rng)

    Create a connectivity matrix for the network.
"""
function make_connectivity_matrix(params::NetworkParams; rng=rng)
    @unpack n,c,connectivity_algorithm = params
    if connectivity_algorithm == "random"
        return rand(rng, Bernoulli(c), n, n)
    elseif connectivity_algorithm == "sonet" 
        return make_sonet_connectivity_matrix(params; rng=rng)
    else 
        error("connectivity_algorithm must be 'random' or 'sonet'")
    end
end

"""
    learn_weights(params::NetworkParams, T::Matrix; rng=rng)

    Learn weights for the network.
""" 
function learn_weights(params::NetworkParams, T::Matrix; rng=rng)
    @unpack m, n, c, τpot, stdp_type = params
    W = make_connectivity_matrix(params; rng=rng)
    @tullio dT[i, j, k] := T[i, j] - T[i, k]

    J = zeros(n, n)
    for i ∈ 1:m
        if stdp_type == "asymmetric"
            jj = @pipe sign.(dT[i,:,:])*exp.(-abs.(dT[i,:,:])./τpot) |> clamp.(_, 0, 1)
        else 
            jj = @pipe exp.(-abs.(dT[i,:,:])./τpot) |> clamp.(_, 0, 1)
        end
        jj[isnan.(jj)] .= 0
        J .+= jj
    end

    return @tullio WJ[i,j] := J[i,j]*W[i,j]
end

# -----------------------------------------------------------------------------
# PART 3: Functions for retrieval phase
# -----------------------------------------------------------------------------

"""
    init_partial_pattern(params::NetworkParams, T::Matrix, Z::Matrix; rng=rng)

    Initialize the network with a partial pattern.
"""
function init_partial_pattern(params::NetworkParams, T::Matrix, Z::Matrix; rng=rng)
    @unpack m, n, b1, bn = params

    Tprime, X = fill(Inf, m, n), zeros(m, n) # initialize
    b1_idx = randperm(rng, m*n)[1:round(Int, m*n*b1)] # indices of correct active pattern
    X[b1_idx] = Z[b1_idx]
    Tprime[b1_idx] = T[b1_idx]

    bn_idx = randperm(rng, m*n)[1:round(Int, m*n*(1-bn))] # indices of spurious active pattern
    bn_idx = setdiff(bn_idx, findall(!iszero, Z[1:length(Z)]))

    Tprime[bn_idx] = rand(rng, Uniform(0, 1), length(bn_idx)) # Why not Gaussian???
    X[bn_idx] = exp.(-Tprime[bn_idx].^2)
    return X,Tprime
end

"""
    run_dynamics(params, Z, T, WJ, X, Tprime; rng)

    Run the dynamics of the network.

    Arguments: 
        - params: NetworkParams. Parameters of the network.
        - Z: Matrix. Ground-truth patterns of the network.
        - T: Matrix. Ground-truth firing times of the network.
        - WJ: Matrix. Weight matrix of the network.
        - X: Matrix. Initial state of the network.
        - Tprime: Matrix. Initial firing times of the network.
        - rng: Random number generator.

    Returns: 
        - act_corrs: Matrix. Correlation between the network's activity and the ground-truth patterns.
        - temp_corrs: Matrix. Correlation between the network's firing times and those in the ground-truth patterns.
        - hit_rates: Matrix. Hit rates of the network over time (rows) and patterns (columns)
        - false_alarm_rates: Matrix. False alarm rates of the network over time (rows) and patterns (columns)
        - d_primes: Matrix. d' values of the network over time (rows) and patterns (columns)
"""


function run_dynamics(params::NetworkParams, Z::Matrix, T::Matrix, WJ::Matrix, X::Matrix, Tprime::Matrix; rng=rng)
    @unpack m, n, nT, g0, g1, g1_std, tau3 = params
    
    # Initialize inhibition scaling based on whether it is fixed or variable
    if g1_std == 0.
        g_inhib = g1 * ones(n)
    else
        g_inhib = rand(rng, BetaMS(g1, g1_std), n)
    end

    # Define functions for calculating correlations, hit rates, false alarm rates, and d'
    getcorr(x,y) = [corrcoef(x[i,:], y[i,:]) for i in 1:m]
    getsdt(x,y) = mapreduce(permutedims, vcat, [sdt_analysis(x[i,:], y[i,:]) for i in 1:m])

    # Initialize matrices to store results
    tmpT = Tprime
    tmpT[findall(!isfinite, tmpT)] .= NaN
    temp_corrs=zeros(nT+1,m)
    act_corrs=zeros(nT+1,m)
    hit_rates=zeros(nT+1,m)
    false_alarm_rates=zeros(nT+1,m)
    d_primes=zeros(nT+1,m)

    # Get initial values of correlations, hit rates, false alarm rates, and d'
    temp_corrs[1,:] = getcorr(tmpT, T)
    act_corrs[1,:] = getcorr(X, Z)
    sdtres = getsdt(X, Z)
    hit_rates[1,:], false_alarm_rates[1,:], d_primes[1,:] = sdtres[:,1], sdtres[:,2], sdtres[:,3]


    δ = 1 
    for k ∈ 1:nT 
        tmpT = fill(Inf, size(X))
        total_inhibition = sum(X, dims=2) # recalculate threshold

        f = 1:m
        for i ∈ 1:m
            for j ∈ 1:n
                begin
                    ix1 = @views findall(!iszero, X[i,:])                  # find indices of active neurons
                    ix2 = @views sortperm(Tprime[i, ix1] .+ δ)             # find indices of active neuron spike times (sorted)
                    sortedT = @view(Tprime[i, ix1])[ix2]                   # sorted spike times 
                    expDsortedT = [1.; exp.(.-diff(sortedT)./tau3)]        # exp of difference between consecutive elements
                    ix3 = @views ix1[ix2]
                end
                
                    
                # Calculate the activity of neuron j at time-step k 
                f = 0.
                for l ∈ 1:length(sortedT)
                        f = f*expDsortedT[l] + WJ[ix3[l], j]
                        if f > (n .* g0 .+ g_inhib[j] .* total_inhibition[i])
                            tmpT[i, j] = sortedT[l] # first time f > θ
                            break
                        end
                    # end
                end
            end 
        end 

        # Update activity of neurons
        X = (tmpT .< Inf) .|> Float64

        # Correlation of before vs after iteration
        temp_corrs[k+1,:] = getcorr(tmpT, T) 
        act_corrs[k+1,:] = getcorr(X, Z)

        # Calculate hit rates, false alarm rates, and d'
        sdtres = getsdt(X, Z)
        hit_rates[k+1,:], false_alarm_rates[k+1,:], d_primes[k+1,:] = sdtres[:,1], sdtres[:,2], sdtres[:,3]
        
        # Update firing times
        Tprime = tmpT
    end 

    return act_corrs, temp_corrs, hit_rates, false_alarm_rates, d_primes
end


# -----------------------------------------------------------------------------
# SIMULATION FUNCTIONS
# -----------------------------------------------------------------------------

"""
    run_sim(run::Int64, params::NetworkParams; rng=rng)

    Run the simulation.

    Arguments: 
        - run: Int64. The run number.
        - params: NetworkParams. The parameters of the network.
        - rng: StableRNG. The random number generator.

    Returns: 
        - DataFrame. The results of the simulation.
"""
function run_sim(run::Int64, params::NetworkParams; rng=rng)
    T, Z = create_patterns(params; rng=rng)
    WJ = learn_weights(params, T; rng=rng)
    X, Tprime = init_partial_pattern(params, T, Z; rng=rng)
    act_corrs,temp_corrs,hit_rates,false_alarm_rates, d_primes = run_dynamics(params, Z, T, WJ, X, Tprime; rng=rng)
    return DataFrame(
        :run => fill(run, (1 + params.nT)*params.m),
        :m => fill(params.m, (1 + params.nT)*params.m), 
        :n => fill(params.n, (1 + params.nT)*params.m),
        :nT => fill(params.nT, (1 + params.nT)*params.m),
        :c => fill(params.c, (1 + params.nT)*params.m),
        :a => fill(params.a, (1 + params.nT)*params.m),
        :b1 => fill(params.b1, (1 + params.nT)*params.m),
        :bn => fill(params.bn, (1 + params.nT)*params.m),
        :σt => fill(params.σt, (1 + params.nT)*params.m),
        :g0 => fill(params.g0, (1 + params.nT)*params.m),
        :g1 => fill(params.g1, (1 + params.nT)*params.m),
        :g1_std => fill(params.g1_std, (1 + params.nT)*params.m),
        :τpot => fill(params.τpot, (1 + params.nT)*params.m),
        :tau3 => fill(params.tau3, (1 + params.nT)*params.m),
        :stdp_type => fill(params.stdp_type, (1 + params.nT)*params.m),
        :connectivity_algorithm => fill(params.connectivity_algorithm, (1 + params.nT)*params.m),
        :αrecip => fill(params.αrecip, (1 + params.nT)*params.m),
        :αconv => fill(params.αconv, (1 + params.nT)*params.m),
        :αdiv => fill(params.αdiv, (1 + params.nT)*params.m),
        :αchain => fill(params.αchain, (1 + params.nT)*params.m),
        :timestep => repeat(0:params.nT, params.m),
        :pattern => repeat(1:params.m, params.nT+1),
        :act_corrs => act_corrs[:],
        :temp_corrs => temp_corrs[:],
        :hit_rates => hit_rates[:],
        :false_alarm_rates => false_alarm_rates[:],
        :d_primes => d_primes[:]
    )
end
