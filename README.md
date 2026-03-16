# SC2Epistasis
This repository contains the code for reproducing the results of the manuscript: [Epistasis and the changing fitness landscapes of SARS-CoV-2](https://doi.org/10.64898/2026.03.12.711354).

## Overview
The manuscript builds on the clade specific fitness effects of amino acid mutations obtained from [SARS2-mut-fitness-v2](https://github.com/neherlab/SARS2-mut-fitness-v2/tree/main), based on the mutation counts from the [UShER tree](https://github.com/jbloomlab/SARS2-mut-fitness/tree/main/results_gisaid_2024-04-24/expected_vs_actual_mut_counts). A minimal epistatic model is inferred by minimization of an objective function defined from the discrepancies of fitness effects across clades.

## Repository structure
- `data/`: contains input files such as mutational effects and clade founder mutations.
- `results/`: contains numerical output files and in the `results/figures/` subdirectory, plots appearing in the manuscript.
- `script/`: scripts to perform the analysis and generate the figures.
- `src/`: source code for major functions definition.
- `Project.toml` and `Manifest.toml`: files containing dependencies for environment reproducibility.
- 
### Main output files
- `results/delta_fit.csv`: dataframe reporting the fitness effects of mutations in pairs of clades whose founders differ at least for one mismatch.
- `results/clade_diff.csv`: dataframe summarizing the mismatches between clade founders.
- `fit_zscore_epi.jld2`: JLD file quantifying the raw epistatic signal.
- `results/jcoup_l1_1em3_3d.csv`: dataframe storing the inferred epistatic couplings. Columns include:
  - $\sigma^{wt}_i$: wildtype amino acid at site i.
  - $i$: residue index of site i.
  - $\sigma_i$: mutant amino acid at site i.
  - $j$: residue index of site j.
  - $\sigma_j$: background amino acid at site j.
  - $J$: inferred coupling $J_{ij}(\sigma_i, \sigma_j; \sigma^{wt}_i)$.

### Main scripts
- `script/delta_fit.jl`: script to generate the dataframe with fitness effects between pairs of clades.
- `script/init_inf.jl`: script to initialize the structs used for inference.
- `script/infer_couplings.jl`: script to infer the epistatic couplings by minimizing the objective 

### Main source code
- `src/types.jl`: definition of tailored structs used for inference.
- `src/numerical_optim.jl` : definition of the optimization algorithm.
  
## Code
The code for performing the analysis relies on the [Julia](https://julialang.org/) programming language. To run the code, follow the steps below.

### Clone the repository

    git clone git@github.com:neherlab/SC2Epistasis.git

### Set up Julia environment
Navigate to the repository directory and open the Julia REPL with the number of desired threads from the shell:

    julia -t <number_of_threads>

Switch to the package manager by pressing `]`. Then, you need to activate and instantiate the environment:

    (@v1.11) pkg> activate .
    (SC2Epistasis) pkg> instantiate

If you use Julia for the first time, you may also need to build the PyCall package in order to use the matplolib plotting library:

    julia> ENV["PYTHON"]=""
    (SC2Epistasis) pkg> build PyCall

You can now import the module `SC2Epistasis` in the Julia REPL:

    julia> using SC2Epistasis

Which contains the functions defined in the `src/` directory. After the first installation, you can run code from the REPL simply by:

    shell> julia -t <number_of_threads>
    (SC2Epistasis) pkg> activate .
    julia> using SC2Epistasis

### Running scripts
You can run the scripts in the `script/` directory from the Julia REPL using the `include` function, i.e.:

    julia> include("script/<script_name>.jl")

Every script imports the required packages and defines additional necessary functions.

## Infer coupling parameters
In order to infer the coupling parameters, run the script:

    julia> include("script/infer_couplings.jl")

This implicitly runs the script `script/init_inf.jl`, which initializes the tailored structs used for inference:
- `data::Vector{Data}`: structs array containing the input data.
- `optx::Vector{Optx}`: structs array containing the optimization parameters.
- `qform::Vector{QForm}`: structs array containing the parameters of the quadratic form.

Then, it runs the inference function `learn(optx, qform, l1_vec; kwargs...)` which returns the optimal value of the objective and modify in place the coupling parameters stored in `optx`.

## References
- Manuscript this repository refers to [Sesta and Neher (2026)](https://doi.org/10.64898/2026.03.12.711354).
- Context dependent estimates of SARS-CoV-2 mutation rates [Haddox et al. (2025)](https://doi.org/10.1093/nar/gkaf503).
- Original approach for estimating fitness effects of mutations [Bloom & Neher (2023)](https://doi.org/10.1093/ve/veae026).
- Potts model for epistasis inference in the context of MSAs of protein sequences [Morcos et al. (2011)](https://doi.org/10.1073/pnas.1111471108).
