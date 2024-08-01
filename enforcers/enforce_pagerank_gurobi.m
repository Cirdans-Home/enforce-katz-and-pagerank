function [Delta,varargout] = enforce_pagerank_gurobi(A,alpha,pihat,v,P,beta)
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
%   OUTPUT Delta perturbation matrix
%          stat structure containing statistics
%          picheck vector of centralities computed with the perturbed
%          matrix, in principle should be equal to the pihat entry
%          rhat shift parameter to enforce positive diagonal entries, may be zero
%          Phat modified probability transition matrix
%

if nargout >= 2
    varargout{1} = struct();
end

% General infos
n              = size(A,1);

% Build the projector
proj                = pattern_projector(P);  % Projector Onto the Pattern of A
reduced_size = size(proj,1);
K                   = commutation(n);
deg                 = A*ones(n);
% Indices for free variables, i.e., indices of diagonal elements after
% pattern projection
kk    = 1:n;
index = kk+(kk-1)*n;

[free_variables,~] = find(proj(:,index));

scaling                          = 1;
scaling_option              = 1; % Do not Change Scaling Options
scaling_direction          = 'l'; % Do not Change Scaling Options

if beta == 1
    % Solving only with Frobenius norm constraints
    % Building the GUROBI model
    model.modelname = 'PageRank-Frobenius';
    model.Q = speye(reduced_size);
    g       = reshape(A,n*n,1);
    g(free_variables) = 0;
    g        = proj*g;
    model.obj = -2*full(g);
    work = spdiags(1./deg,0,n,n)*pihat;
    L       = [kron(work.',speye(n))*(K*proj.');...
        kron(ones(n,1).',speye(n))*proj.'];
    b       = [ (1/alpha).*(pihat -(1-alpha).*v) - A.'*(spdiags(1./deg,0,n,n)*pihat)+kron(work.',speye(n))*(K*(proj.'*g) )  ;...
        kron(ones(n,1).',speye(n))*(proj.'*g) ];
    model.A   = L;
    model.rhs = full(b);
    model.lb  = zeros(size(model.A,2),1);
    model.lb(free_variables) = -Inf;
    model.sense = repmat('=',size(model.A,1),1);
    % Call GUROBI on the problem to solve
    result = gurobi(model);
    if nargout >= 2
        varargout{1}.time       = result.runtime;
        varargout{1}.iter       = result.itercount;
        varargout{1}.baritercount = result.baritercount;
    end

    % Recover Matrix and Desired Ranking
    [ival,jval,~] = find(P);
    Delta = result.x - g;
    Delta(abs(Delta)<1e-13)=0;
    Delta = sparse(ival,jval,Delta,n,n);
    % Compute the optimized PageRank centrality
    if nargout >= 3
        I              = speye(n,n);
        rhat         = min(diag( spdiags(1./deg,0,n,n)*(A+Delta) ));

        if rhat < 1e-8
            rhat          = 1- alpha*rhat;
            r               = rhat;
            alphahat   = 1- (1-alpha)/r;
            Phat         = 1/(r-1+alpha).*(alpha*(spdiags(1./deg,0,n,n)*(A+Delta))+(r-1).*speye(n) );
            picheck    =  (I - alphahat*Phat.')\((1-alphahat).*v) ;
            varargout{2} = picheck/sum(picheck);
	else
            rhat      = 0;
            picheck =  (I - alpha*(spdiags(1./deg,0,n,n)*(A+Delta)).')\((1-alpha).*v) ;
            varargout{2} = picheck/sum(picheck);
	end
        if nargout >= 4
            varargout{3} = rhat;
        end
        if nargout >= 5
            varargout{4} = Phat;
        end

    end

else
    % Solving with sparsity constraints
    % Building the GUROBI model
    model.name = 'PageRank-l1-fro';
    tau = (1-beta)/beta;
    Q   = 2*speye(reduced_size);
    work = spdiags(1./deg,0,n,n)*pihat;
    L       = [kron(work.',speye(n))*(K*proj.');...
        kron(ones(n,1).',speye(n))*proj.'];
    c       = reshape(A,n*n,1);
    c(free_variables) = 0;
    c        = proj*c;
    b       = [   (1/alpha).*(pihat -(1-alpha).*v) - A.'*(spdiags(1./deg,0,n,n)*pihat)+kron(work.',speye(n))*(K*(proj.'*c) )  ;...
        kron(ones(n,1).',speye(n))*(proj.'*c) ];
    b        = [b;-c];
    model.Q  = blkdiag(Q,sparse(reduced_size,reduced_size),sparse(reduced_size,reduced_size));
    g        = [-2.*c; tau.*ones(reduced_size,1); tau.*ones(reduced_size,1)];
    model.obj = full(g);
    L   = [L, sparse(2*n,reduced_size), sparse(2*n,reduced_size);...
        - speye(reduced_size), speye(reduced_size),      -speye(reduced_size)];
    if (scaling == 1)
        DD = Scale_the_problem(L,scaling_option,scaling_direction);
        L = spdiags(DD,0,size(L,1),size(L,1))*L;  % Apply the left scaling.
        b = b.*DD;
    end
    model.A     = L;
    model.rhs   = full(b);
    model.lb  = zeros(size(model.A,2),1);
    model.lb(free_variables) = -Inf;
    model.sense = repmat('=',size(model.A,1),1);

    % Call GUROBI on the problem to solve
    result = gurobi(model);
    if nargout >= 2
        varargout{1}.time       = result.runtime;
        varargout{1}.iter       = result.itercount;
        varargout{1}.baritercount = result.baritercount;
    end

    % Recover Matrix and Desired Ranking
    [ival,jval,~] = find(P);
    Delta    = result.x(1:reduced_size) -c;
    Delta(abs(Delta)<1e-13)=0;
    Delta    = sparse(ival,jval,Delta,n,n);
    % Compute the optimized PageRank centrality
    if nargout >= 3
        I              = speye(n,n);
        rhat         = min(diag( spdiags(1./deg,0,n,n)*(A+Delta) ));
        if rhat < 1e-8
            rhat          = 1- alpha*rhat;
            r               = rhat;
            alphahat   = 1- (1-alpha)/r;
            Phat         = 1/(r-1+alpha).*(alpha*(spdiags(1./deg,0,n,n)*(A+Delta))+(r-1).*speye(n) );
            picheck     =  (I - alphahat*Phat.')\((1-alphahat).*v) ;
            varargout{2} = picheck/sum(picheck);
	else
            rhat      = 0;
            picheck =  (I - alpha*(spdiags(1./deg,0,n,n)*(A+Delta)).')\((1-alpha).*v) ;
            varargout{2} = picheck/sum(picheck);
	end
        if nargout >= 4
            varargout{3} = rhat;
        end
        if nargout >= 5
            varargout{4} = Phat;
        end
    end
end

end

%% Computational Routines
% Commutation Matrix
function [K] = commutation(n)
I = reshape(1:n*n, [n, n]); % initialize a matrix of indices of size(A)
I = I';    % Transpose it
I = I(:); % vectorize the required indices
K = speye(n*n); % Initialize an identity matrix
K = K(I,:); % Re-arrange the rows of the identity matrix
end

function [proj] = pattern_projector(P)
k     = find(P);
N    = size(P,1);
proj = sparse(1:1:length(k),k,ones(length(k),1), length(k),N*N );
end

function [nr_res_p,nr_res_d,mu] = IPM_Res(c,A,A_tr,b,Q,x,y,z,pos_vars, ...
    num_of_pos_vars)
nr_res_p = b-A*x;                 % Non-regularized primal residual
nr_res_d = c-A_tr*y-z + Q*x;      % Non-regularized dual residual
if num_of_pos_vars > 0
    mu = (x(pos_vars)'*z(pos_vars))/num_of_pos_vars;
else
    mu = 0 ;
end
end

function [D,D_L] = Scale_the_problem(A,scale_option,direction)
% ==================================================================================================================== %
% [D] = Scale_the_problem(A):
% -------------------------------------------------------------------------------------------------------------------- %
% This function, takes as an input a sparse matrix A, representing the constraints of a quadratic progrmaming problem.
% It checks whether the matrix is well scaled, and if not, it applies some linear transformations to the matrix, in
% order to improve its numerical properties. The method return a diagonal matrix D, which the user should
% use in order to recover the solution of the initial problem (after solving the scaled one).
% The optional parameter: scale_option. This parameter can take 3 values:
%   (i)   scale_option = 0, no scaling will be applied.
%   (ii)  scale_option = 1, iterative geometric scaling will be used.
%   (iii) scale_option = 2, equilibrium scaling is employed.
%   (iv)  scale_option = 3, nearest power of 2 scaling is used.
%   (v)   scale_option = 4, mixed strategy, based on the properties of the respective row/column.
% The optional parameter: direction. This parameter can take 3 values:
%   (i)   direction = 'r', right scaling (default).
%   (ii)  direction = 'l', left scaling.
% For more information about these scaling choices, the reader is refered to:
%                                   https://en.wikibooks.org/wiki/GLPK/Scaling
%
% Author: Spyridon Pougkakiotis.
% ==================================================================================================================== %
if (nargin < 2 || isempty(scale_option))
    scale_option = 1; % Set geometric scaling as the default option.
end
if (nargin < 3 || isempty(direction))
    direction = 'r';
end
if (direction == 'l')
    A = A';
end
pos_A = abs(A);       % Need it to identify non-zero elements.
D = zeros(size(A,2),1);
D_L = [];
pos_ind = pos_A > 0;
% ================================================================================================================ %
% Based on the input parameters, build the desired scaling factor.
% ---------------------------------------------------------------------------------------------------------------- %
if (max(max(pos_A)) <= 10 && min(min(abs(A(pos_ind))))>= 0.1 || scale_option == 0)
    fprintf('No scaling necessary.\n'); % Well scaled or no scale.
    D = ones(size(A,2),1);
elseif (scale_option == 1) % Geometric scaling (applied on columns for computational efficiency).
    fprintf('The constraint matrix is scaled. Geometric scaling is employed.\n');
    for j = 1:size(A,2)
        rows = pos_A(:,j) > 0; % Find all non-zero elements for this column
        if (any(rows))
            %size(pos_A(rows,j))
            maximum = max(pos_A(rows,j));
            minimum = min(pos_A(rows,j));
            if (maximum*minimum > 10^(-12) && maximum*minimum < 10^(12))
                D(j) = 1/sqrt(maximum*minimum);
            else
                D(j) = 1;
            end
        else
            D(j) = 1; % Extreme case, where one column is all zeros.
        end
    end
elseif (scale_option == 2) % Equilibrium scaling (applied on columns for efficiency).
    fprintf('The constraint matrix is scaled. Equilibrium scaling is employed.\n');
    for j = 1:size(A,2)
        rows = pos_A(:,j) > 0; % Find all non-zero elements for this column.
        maximum = max(pos_A(rows,j));
        if (maximum > 10^(-6))
            D(j) = 1/maximum;
        else
            D(j) = 1;
        end
    end
elseif (scale_option == 3) % Nearest power of 2 scaling + geometric scaling (avoiding rounding errors).
    fprintf('The constraint matrix is scaled. Geometric scaling with nearest power of 2 is employed.\n');
    for j = 1:size(A,2)
        rows = pos_A(:,j) > 0; % Find all non-zero elements for this column.
        if (any(rows))
            maximum = max(pos_A(rows,j));
            minimum = min(pos_A(rows,j));
            p = nextpow2(sqrt(maximum*minimum));
            if (maximum*minimum > 10^(-12) && maximum*minimum < 10^(12))
                D(j) = 1/(2^(p-1));
            else
                D(j) = 1;
            end
        else
            D(j) = 1; % Extreme case, where one column is all zeros.
        end
    end
elseif (scale_option == 4)
    fprintf('The constraint matrix is scaled. Mixed scaling is employed.\n');
    for j = 1:size(A,2)
        rows = pos_A(:,j) > 0; % Find all non-zero elements for this column
        if (any(rows))
            %size(pos_A(rows,j))
            maximum = max(pos_A(rows,j));
            minimum = min(pos_A(rows,j));
            if (maximum > 10^3 && minimum < 10^(-3))
                p = nextpow2(sqrt(maximum*minimum));
                D(j) = 1/(2^(p-1));
            elseif (1/minimum > maximum && minimum > 10^(-6))
                p = nextpow2(minimum);
                D(j) = 1/(2^(p-1));
            elseif (maximum < 10^(6))
                p = nextpow2(maximum);
                D(j) = 1/(2^(p-1));
            else
                D(j) = 1;
            end
        else
            D(j) = 1; % Extreme case, where one column is all zeros.
        end
    end
end

% ================================================================================================================ %
end
