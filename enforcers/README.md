# Enforcers

This folder contains the individual routines that solve the different 
instances of the problem with the algorithms described in the Paper.

- `enforce_katz.m`:
   ```matlab
%ENFORCE_KATZ Given a matrix A and a parameter alpha finds a matrix Delta
%with pattern P such that:
%   (I - \alpha(A+\Delta))^{-1}1 = \hat{mu}
%and \beta \|\Delta\|_F^2 + (1-\beta) \|\Delta\|_1 minimal
%   INPUT A sparse matrix
%         alpha such that (I - \alpha A) is invertible
%         muhat desired vector of scores with entries >= 1
%         P pattern matrix
%         beta scalar in (0,1) regulating the objective function
%         tol tolerance for the IPM Solver
%   OUTPUT Delta perturbation matrix
%          stat structure containing statistics
%          mucheck vector of centralities computed with the perturbed
%          matrix, in principle should be equal to the muhat entry
```
- `enforce_pagerank_closedform.m`
   ```matlab
%ENFORCE_PAGERANK_CLOSEDFORM Produces the perturbation in closed form
%following the results given in Proposition 4.1 and Proposition 4.2 of the
%paper. This is a perturbation, having the same sparsity pattern of A plus
%the diagonal, and is not guaranteed to be minimal
%   INPUT: A adjacency matrix of the network
%          alpha PageRank parameter
%          v preference vector
%          pihat desired steady state
%   OUTPUT: Phat modified matrix P
%           Delta perturbation from Proposition 4.1
%           fronorm = ||G - Ghat||_F
```

> [!TIP]
> These versions are useful to have the code in a standalone way, to work 
> on the single parts of the algorithm and interact with the various parts 
> it can be more useful to consult the specific folders for the two 
> problems in which the code is separated into the various functions 
> and more readable.