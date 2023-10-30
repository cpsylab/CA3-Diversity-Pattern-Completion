# CA3-Diversity-Pattern-Completion

Code for the paper *The Impact of Electrophysiological Diversity on Pattern Completion in Lithium Nonresponsive Bipolar Disorder: A Computational Modelling Approach* 

## Directory Structure 

- `main.jl`: Primary script for running the experiments 
- `ca3.jl`: Functions defining network, encoding, and retrieval operations
- `pedagogical-plots.jl`: Plots of the STDP function and other transformations for figures
- `hypervshyposim.jl`: Script for analysis of whether hyper or hypoexctable variation drives effects of variability on pattern completion 
- `analysis.jl`: Scripts for analyzing data 

## Running 

``` 
julia --project pedagogical-plots.jl
julia --project --threads NTHREADS main.jl
julia --project --threads NTHREADS hypervshyposim.jl
julia --project analysis.jl
```

To obtain the figures in the paper, chag

