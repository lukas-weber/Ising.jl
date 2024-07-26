# Ising
[![DOI](https://zenodo.org/badge/558064408.svg)](https://zenodo.org/doi/10.5281/zenodo.12981007)


This is the reference implementation of the Metropolis Markov-Chain Monte Carlo algorithm for the 2D Ising model using the [Carlo.jl](https://github.com/lukas-weber/Carlo.jl) framework.

The file `src/Ising.jl` contains the implementation of the `Carlo.AbstractMC` interface, which is all that is needed to implement a Carlo Monte Carlo algorithm.

## Installation
```julia
julia> using Pkg; Pkg.add("Ising"); Pkg.add("Carlo")
```

## Example
There is an example job script `example/job` to show how to run simulations using Carlo. Simply execute it to access the Carlo command-line interface.

```bash
cd example

julia --project example_job.jl run
```

If everything went well, you should find the file `example/job.results.json` containing the means and errorbars of your calculation.

For unfinished jobs, you can at any point use

```bash
julia --project example_job.jl status
```

to retrieve progress information and 

```bash
julia --project example_job.jl merge
```

to merge the data already collected into a results file.

To start over from scratch, remember to `julia --project example_job.jl delete` the existing data or use the `run -r` flag. For more information, run `julia --project example_job.jl --help`.
