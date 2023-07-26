# Ising

This is the reference implementation of the Metropolis Markov-Chain Monte Carlo algorithm for the 2D Ising model using the [Carlo](https://github.com/lukas-weber/Carlo.jl) framework. 

The file `src/Ising.jl` contains the implementation of the `Carlo.AbstractMC` interface, which is all that is needed to implement a Carlo Monte Carlo algorithm.

## Installation
```julia
julia> using Pkg; Pkg.add("Ising"); Pkg.add("Carlo")
```

## Example
There is an example job script `example/job` to show how to run simulations using Carlo. Simply execute it to access the Carlo command-line interface.

```bash
cd example
./job --help

./job run
```

If everything went well, you should find the file `example/job.results.json` containing the means and errorbars of your calculation.