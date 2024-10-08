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
- `enforce_katz_gurobi.m`
   ```matlab
   %ENFORCE_KATZ_GUROBI Given a matrix A and a parameter alpha finds a matrix 
   %Delta with pattern P such that:
   %   (I - \alpha(A+\Delta))^{-1}1 = \hat{mu}
   %and \beta \|\Delta\|_F^2 + (1-\beta) \|\Delta\|_1 minimal
   %   INPUT A sparse matrix
   %         alpha such that (I - \alpha A) is invertible
   %         muhat desired vector of scores with entries >= 1
   %         P pattern matrix
   %         beta scalar in (0,1) regulating the objective function
   %   OUTPUT Delta perturbation matrix
   %          stat structure containing statistics
   %          mucheck vector of centralities computed with the perturbed
   %          matrix, in principle should be equal to the muhat entry
   %
   ```
- `enforce_pagerank.m`
   ```matlab
   %ENFORCE PAGERANK Given a matrix A and a parameter alpha finds a matrix Delta
   %with pattern P such that:
   %    (I- alpha ( diag(A1)^{-1}(A+Delta))^T)pihat = (1-alpha)v
   %    Delta 1 = 0
   %    offdiagonal(A+Delta) >=0
   %    --> and  beta \|\Delta\|_F^2 + (1-beta) \|\Delta\|_1 minimal
   %   INPUT A sparse matrix
   %         alpha such that A is irreducible
   %         pihat desired vector of pagerank with e^Tpihat=1
   %         P pattern matrix
   %         beta scalar in (0,1) regulating the objective function
   %         tol tolerance for the IPM Solver
   %   OUTPUT Delta perturbation matrix
   %          stat structure containing statistics
   %          picheck vector of centralities computed with the perturbed
   %          matrix, in principle should be equal to the pihat entry
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
- `enforce_pagerank_gurobi.m`
  ```matlab
  %ENFORCE PAGERANK Given a matrix A and a parameter alpha finds a matrix Delta
  %with pattern P such that:
  %    (I- alpha ( diag(A1)^{-1}(A+Delta))^T)pihat = (1-alpha)v
  %    Delta 1 = 0
  %    offdiagonal(A+Delta) >=0
  %    --> and  beta \|\Delta\|_F^2 + (1-beta) \|\Delta\|_1 minimal
  %   INPUT A sparse matrix
  %         alpha such that A is irreducible
  %         pihat desired vector of pagerank with e^Tpihat=1
  %         P pattern matrix
  %           beta scalar in (0,1) regulating the objective function
  %   OUTPUT Delta perturbation matrix
  %          stat structure containing statistics
  %          picheck vector of centralities computed with the perturbed
  %          matrix, in principle should be equal to the pihat entry
  %          rhat shift parameter to enforce positive diagonal entries, may be zero
  %          Phat modified probability transition matrix
  %
  ```


> [!TIP]
> These versions are useful to have the code in a standalone way, to work 
> on the single parts of the algorithm and interact with the various parts 
> it can be more useful to consult the specific folders for the two 
> problems in which the code is separated into the various functions 
> and more readable.

> [!CAUTION]
> The `enforce_katz_gurobi.m` uses the [Gurobi](https://www.gurobi.com) solver for the optimization problem
> it is possible to get an Academic License for it. If you don't want to use, the 
> `enforce_katz.m` implementation is written in pure Matlab and doesn't require external 
> tools.
