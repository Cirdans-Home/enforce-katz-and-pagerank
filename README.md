# Enforce Katz and PageRank

Code for enforcing the Katz and PageRank centrality scores.

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

