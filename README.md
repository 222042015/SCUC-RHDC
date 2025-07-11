# Relax-and-Cut for Temporal SCUC Decomposition

This repository contains the implementation for the manuscript titled: "Relax-and-Cut for Temporal SCUC Decomposition".

To reproduce the results presented in the manuscript, execute
```bash
bash scripts.sh
```

To run with a single instance, execute
```
julia --threads=1 test_cb.jl -dataset $instance -nInt 6 -nCont 6 -stepsize 6 -threads 1
```

The code for visualizing the generated results is located in the Jupyter Notebook: `load_result.ipynb`.

This work builds upon the "UnitCommitment.jl" package.
* **Alinson S. Xavier, Aleksandr M. Kazachkov, Ogün Yurdakul, Jun He, Feng Qiu**. "UnitCommitment.jl: A Julia/JuMP Optimization Package for Security-Constrained Unit Commitment (Version 0.4)". Zenodo (2024). [DOI: 10.5281/zenodo.4269874](https://doi.org/10.5281/zenodo.4269874).
