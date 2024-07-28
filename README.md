# Enforcing Katz and PageRank Centrality Measures in Complex Networks

We investigate the problem of **enforcing a desired centrality measure** in complex networks, 
while still keeping the original pattern of the network. Specifically, by representing the 
network as a graph with suitable nodes and weighted edges, we focus on computing the smallest 
perturbation on the weights required to obtain a prescribed **PageRank** or **Katz centrality index** 
for the nodes. By quantifying the least perturbation amount, we gain insights into the properties 
of the network and can make targeted modifications to achieve specific centrality objectives. 
Our approach relies on optimization procedures that scale with the size of the network, 
enabling the exploration of different scenarios and the alteration of network structure 
and dynamics.

If you use this code or any of the ideas here, please cite the paper:
```bibtex
@misc{EnforcePRKATZ
   title={Enforcing Katz and PageRank Centrality Measures in Complex Networks},
   author={S. Cipolla, F. Durastante, B. Meini}
}
```

## Collaborators

- Stefano Cipolla [:email:](s.cipolla@soton.ac.uk)[:earth_africa:](https://stefanocipolla.github.io/)
- Fabio Durastante [:email:](fabio.durastante@unipi.it)[:earth_africa:](https://fdurastante.github.io)
- Beatrice Meini [:email:](beatrice.meini@unipi.it)[:earth_africa:](https://people.dm.unipi.it/meini/)

## Getting the code

To obtain the code, you can use Git with the following command:
```bash
git clone --recurse-submodule git@github.com:Cirdans-Home/enforce-katz-and-pagerank.git
```
This command will not only clone the main repository but also initialize and clone all its submodules; see the next section for some details.


### Target Stationary Distribution Problem

We make a comparison with the results from

> "Assigning Stationary Distributions to Sparse Stochastic Matrices", Nicolas Gillis and Paul Van Dooren, arXiv:2312.16011, 2023. [https://arxiv.org/abs/2312.16011](arxiv:2312.16011)

Their [code](https://gitlab.com/ngillis/TSDP) is included as a submodule of this repository, hence if you want to run the related comparison tests you need the clone command with the `--recurse-submodules` option.

### GUROBI (Optional)

Some versions of the solvers can use the [Gurobi](https://www.gurobi.com) solver for the optimization problem. It is possible to get an Academic License for it. If you don't want to use it, this is not required since implementations are written in pure Matlab and don't require external tools.

### Replicating numerical experiments

The `NumericalExamples` folder contains the necessary scripts and instructions to replicate the examples included in the paper.
