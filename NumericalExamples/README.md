# Numerical Examples

This folder contains the MATLAB scripts needed to reproduce the numerical experiments included in the paper.

- `small_example.m` produces the example on the Sioux Falls Road Network for enforcing Katz Centrality Measure,
- `heatmaptest.m` produces the *preliiminary results* for the Katz centrality measure and saves them as a `.mat` file that can then be used to produce the figures in the paper with the `plot_heatmaps.m` script,
- `Test_Katz_Forcing.m` produces the Scenarios 1 and 2 test for the Katz centrality measure
- `Test_PR_Forcing.m` produces the Scenarios 1 and 2 test for the PageRank centrality measure
- `comparison_with_TSDP.m` runs the comparison with the TSDP algorithm.


> [!IMPORTANT]  
> The tests use code that is contained in the `../enforcers` and `../Centrality_Forcing_IPM` folders. 
> Make sure you have added them to the MATLAB path (`addpath(<path>)`).

> [!WARNING]
> Comparison with the TSDP algorithm requires cloning the repository with the 
> `--recurse-submodules` option, see the information in the 
> [main README](https://github.com/Cirdans-Home/enforce-katz-and-pagerank/blob/main/README.md) of the repository.

## Toeplitz

The experiments were performed on one of the compute nodes of the Toeplitz cluster at the Green Data 
Center of the University of Pisa. That is, it is one of the four Nodes Intel(R) Xeon(R) CPU E5-2650 
v4 @ 2.20GHz with 2 threads per core, 12 cores per socket and 2 sockets with 256 GB of memory of which 
250 available.

### Batch scripts

The folder also contains some batch scripts (`.sh`) for the SLURM queue manager that are used to run 
the examples on the compute nodes of the Toeplitz cluster at the Gree Data Center of the University of 
Pisa. That is, they queue the tests by doing `sbatch run_<something>.sh`.


