# Coupled Boolean Networks:

`coupledboolnet` allows users to create and simulate single or multicellular-coupled Boolean networks with different interaction rules (lienar threshold or Ising model).   


## How to use

Users supply connectivty value ($k$) from command line:

    python -m coupleboolnet [k-value]

By default, if `importsavedata = True`, `coupledboolnet` searches specfied directory path under `coupledboolnet/data/` directory.  

Otherwise, the default setting is a tissue of 10-by-10 cells with 10-genes with perturbation, and Ising interaction.


## Contents 

- Edit `bnsettings.py` for extensive paramter of control with boolean networks. Options include perturbation, temperature $T$, interaction strength $h$. 

- Run methods from `runanalysis2.py` separately for visualization tools of the simulation outcomes.

- `coupledboolent/kamiak` contains scripts for running on SLURM supported clusters.


## Work to be done:

1. Run on Kamiak: Done
2. Visualization script: Done
3. Python Docstring
