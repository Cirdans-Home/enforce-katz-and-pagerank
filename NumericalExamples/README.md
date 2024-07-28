# Numerical Examples

This folder contains the MATLAB scripts needed to reproduce the numerical experiments included in the paper.

- `heatmaptest.m` produces the *preliiminary results* in Section 6.4.2 and saves them as a `.mat` file.


> [!IMPORTANT]  
> The tests use code that is contained in the `../enforcers` and `../Centrality_Forcing_IPM` folders. 
> Make sure you have added them to the MATLAB path (`addpath(<path>)`).

## Toeplitz

The experiments were performed on one of the compute nodes of the Toeplitz cluster at the Green Data 
Center of the University of Pisa. That is, it is one of the four Nodes Intel(R) Xeon(R) CPU E5-2650 
v4 @ 2.20GHz with 2 threads per core, 12 cores per socket and 2 sockets with 256 GB of memory of which 
250 available.

### Batch scripts

The folder also contains some batch scripts (`.sh`) for the SLURM queue manager that are used to run 
the examples on the compute nodes of the Toeplitz cluster at the Gree Data Center of the University of 
Pisa. That is, they queue the tests by doing `sbatch run_<something>.sh`.


