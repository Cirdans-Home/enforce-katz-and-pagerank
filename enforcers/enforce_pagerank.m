function [Delta,varargout] = enforce_pagerank(A,alpha,pihat,v,P,beta,tol)
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
%          rhat shift parameter to enforce positive diagonal entries, may be zero
%          Phat modified probability transition matrix
%

if nargout >= 2
    varargout{1} = struct();
end

% General infos
n              = size(A,1);
print_mode     = 3;

% Build the projector
proj                = pattern_projector(P);  % Projector Onto the Pattern of A
reduced_size = size(proj,1);
K                    = commutation(n);
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
    H      = 2*speye(reduced_size);
    work = spdiags(1./deg,0,n,n)*pihat;
    L       = [kron(work.',speye(n))*(K*proj.');...
        kron(ones(n,1).',speye(n))*proj.'];
    g       = reshape(A,n*n,1);
    g(free_variables) = 0;
    g        = proj*g;

    b       = [   (1/alpha).*(pihat -(1-alpha).*v) - A.'*(spdiags(1./deg,0,n,n)*pihat)+kron(work.',speye(n))*(K*(proj.'*g) )  ;...
        kron(ones(n,1).',speye(n))*(proj.'*g) ];


    if (scaling == 1)
        DD = Scale_the_problem(L,scaling_option,scaling_direction);
        L = spdiags(DD,0,size(L,1),size(L,1))*L;  % Apply the left scaling.
        b = b.*DD;
    end

    % Running the solver
    IterStruct     = struct();
    IterStruct.Fact                   = 'chol';
    rho               = 1e-9;
    delta            = rho;
    pc_mode     = 2;
    tic;
    [Delta,~,~,Info] = PPM_IPM(-2*g, ...
        L, ...
        b, ...
        H, ...
        free_variables, ...
        tol, ...
        200,...
        pc_mode, ...
        print_mode, ...
        IterStruct, ...
        rho, ...
        delta);
    elapstime = toc;
    if nargout >= 2
        varargout{1}.time           = elapstime;
        varargout{1}.opt            = Info.opt;
        varargout{1}.IPMiter      = Info.IPM_It;
        varargout{1}.primalres  = [Info.NatRes.primal];
        varargout{1}.dualres     = [Info.NatRes.dual];
        varargout{1}.compl       = [Info.NatRes.compl];
    end

    % Recover Matrix and Desired Ranking
    [ival,jval,~] = find(P);
    Delta(abs(Delta)<1e-13)=0;
    Delta = Delta -g;
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
            varargout{2}     =  (I - alphahat*Phat.')\((1-alphahat).*v) ;
        else
            rhat      = 0;
            varargout{2} =  (I - alpha*(spdiags(1./deg,0,n,n)*(A+Delta)).')\((1-alpha).*v) ;
        end
        if nargout >= 4
            varargout{3} = rhat;
        end
        if nargout >= 5
            varargout{4} = Phat;
        end

        % alphahat = 1- (1-alpha)/r;
        % Phat        =  1/(r-1+alpha).*(alpha*(spdiags(1./deg,0,n,n)*(A+Delta))+(r-1).*speye(n) );
        % varargout{2} = (I - alphahat*Phat.')\((1-alphahat).*v) ;
    end

else
    constr_var = setdiff(1:1:reduced_size, free_variables);
    % Solving with sparsity constraints
    tau = (1-beta)/beta;
    Q   = 2*speye(reduced_size);
    work = spdiags(1./deg,0,n,n)*pihat;
    L    = [kron(work.',speye(n))*(K*proj.');...
        kron(ones(n,1).',speye(n))*proj.'];


    c       = reshape(A,n*n,1);
    c(free_variables) = 0;
    c        = proj*c;

    b       = [   (1/alpha).*(pihat -(1-alpha).*v) - A.'*(spdiags(1./deg,0,n,n)*pihat)+kron(work.',speye(n))*(K*(proj.'*c) )  ;...
        kron(ones(n,1).',speye(n))*(proj.'*c) ];
    b        = [b;-c(constr_var)];
    H       = blkdiag(Q,sparse(length(constr_var),length(constr_var)),sparse(length(constr_var),length(constr_var)));
    g        = [-2.*c; tau.*ones(length(constr_var),1); tau.*ones(length(constr_var),1)];
    II        = speye(reduced_size);
    L        = [L, sparse(2*n,length(constr_var)), sparse(2*n,length(constr_var));...
        - II(constr_var,:),     speye(length(constr_var)),      -speye(length(constr_var))];


    if (scaling == 1)
        DD = Scale_the_problem(L,scaling_option,scaling_direction);
        L = spdiags(DD,0,size(L,1),size(L,1))*L;  % Apply the left scaling.
        b = b.*DD;
    end

    % Running the solver
    IterStruct      = struct();
    IterStruct.Fact = 'chol';
    rho             = 1e-9;
    delta           = rho;
    pc_mode         = 2;
    tic;
    [Delta,~,~,Info] = PPM_IPM(g,L,b,H,free_variables,tol,200,...
        pc_mode,print_mode,IterStruct,rho,delta);
    elapstime = toc;
    if nargout >= 2
        varargout{1}.time          = elapstime;
        varargout{1}.opt            = Info.opt;
        varargout{1}.iter            = Info.ExIt;
        varargout{1}.IPMiter     = Info.IPM_It;
        varargout{1}.primalres  = [Info.NatRes.primal];
        varargout{1}.dualres     = [Info.NatRes.dual];
        varargout{1}.compl       = [Info.NatRes.compl];
    end
    % Recover Matrix and Desired Ranking
    [ival,jval,~] = find(P);
    Delta = Delta(1:reduced_size);
    Delta(abs(Delta)<1e-13)=0;
    Delta    = Delta -c;
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
            varargout{2}     =  (I - alphahat*Phat.')\((1-alphahat).*v) ;
        else
            rhat      = 0;
            varargout{2} =  (I - alpha*(spdiags(1./deg,0,n,n)*(A+Delta)).')\((1-alpha).*v) ;
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


function [x,y,z,Info] = PPM_IPM(c,A,b,Q,free_variables,...
    tol,maxit,pc,printlevel,IterStruct,rho,delta)
%  IPM   Primal-dual Regularized interior-point method with decoupled variables.
%
%  This is the driver function of an IPM for solving the
%  quadratic programming problem
%
%    min c'x + 1/2*x'Qx  subject to A*x=b, x_C>=0.       (1)
%
%  printlevel options:
%  0: turn off iteration output
%  1: print primal and dual residual and duality measure
%  2: print centering parameter and step length
%  3: print residuals in the solution of the step equations
%  printlevel Default: 1.
% ==================================================================================================================== %
% Parameter filling and dimensionality testing.
% -------------------------------------------------------------------------------------------------------------------- %
[m, n] = size(A);
% Make sure that b and c are column vectors of dimension m and n.
if (size(b,2) > 1); b = (b)'; end
if (size(c,2) > 1); c = (c)'; end
if (~isequal(size(c),[n,1]) || ~isequal(size(b),[m,1]) )
    error('problem dimension incorrect');
end

% Make sure that A is sparse and b, c are full.
if (~issparse(A)); A = sparse(A); end
if (~issparse(Q)); Q = sparse(Q); end
if (issparse(b));  b = full(b);   end
if (issparse(c));  c = full(c);   end

% Set default values for missing parameters.
if (nargin < 5 || isempty(free_variables)); free_variables = []; end
if (nargin < 6 || isempty(tol));            tol = 1e-4;          end
if (nargin < 7 || isempty(maxit));          maxit = 100;         end
if (nargin < 8 || isempty(pc));             pc = 1;              end
if (nargin < 9 || isempty(printlevel));     printlevel = 1;      end
% ======================================================================= %
% Initialization of the structure for convergence hystory
Info = struct();

% ======================================================================= %
% Initialization - Mehrotra's Initial Point for QP:
% Choose an initial starting point (x,y,z). For that, we ignore the
% non-negativity constraints, as well as the regularization variables and
% solve the relaxed optimization problem (which has a closed form
% solution). Then, we shift the solution, to respect the non-negativity
% constraints. The point is expected to be well centered.
% ----------------------------------------------------------------------- %
A_tr = A';              % Store the transpose for computational efficiency.
pos_vars = setdiff((1:n)',free_variables);
num_of_pos_vars = size(pos_vars,1);
e_pos_vars = ones(num_of_pos_vars,1);

% Turn off Predictor-Corrector when PMM is only running.
if (num_of_pos_vars == 0 && pc ~= 1)
    pc = 1;
end

% ======================================================================= %
% Use PCG to solve two least-squares problems for efficiency (along with
% the Jacobi preconditioner).
% ----------------------------------------------------------------------- %
D = sum(A.^2,2) + 10;
Jacobi_Prec = @(x) (1./D).*x;
NE_fun = @(x) (A*(A_tr*x) + 10.*x);
[x,~] = pcg(NE_fun,b,10^(-8),min(1000,m),Jacobi_Prec);
x = A_tr*x;
[y,~] = pcg(NE_fun,A*(c+Q*x),10^(-8),min(1000,m),Jacobi_Prec);
z = c+ Q*x - A_tr*y;
% ======================================================================= %
if (norm(x(pos_vars)) <= 10^(-4))
    x(pos_vars) = 0.1.*ones(num_of_pos_vars,1); % 0.1 is chosen arbitrarily
end

if (norm(z(pos_vars)) <= 10^(-4))
    z(pos_vars) = 0.1.*ones(num_of_pos_vars,1); % 0.1 is chosen arbitrarily
end

delta_x = max(-1.5*min(x(pos_vars)),0);
delta_z = max(-1.5*min(z(pos_vars)), 0);
temp_product = (x(pos_vars) + (delta_x.*e_pos_vars))'*(z(pos_vars) ...
    + (delta_z.*e_pos_vars));
delta_x_bar = delta_x + (0.5*temp_product)/(sum(z(pos_vars),1) + ...
    num_of_pos_vars*delta_z);
delta_z_bar = delta_z + (0.5*temp_product)/(sum(x(pos_vars),1) + ...
    num_of_pos_vars*delta_x);

z(pos_vars) = z(pos_vars) + delta_z_bar.*e_pos_vars;
x(pos_vars) = x(pos_vars) + delta_x_bar.*e_pos_vars;
z(free_variables) = 0;

if (issparse(x));  x = full(x); end
if (issparse(z));  z = full(z); end
if (issparse(y));  y = full(y); end
% ======================================================================= %
% PPM parameters initialization

iter        = 0;
PPM_red     = 0.7;
IPM_maxit   = 70;
IPM_Tot_It = 0;

xk = x;
yk = y;
zk = z;

% ********************************************************************
% PPM-IPM Main-Loop
% =============
[nr_res_p,nr_res_d,mu] = IPM_Res(c,A,A_tr,b,Q,x,y,z,pos_vars, ...
    num_of_pos_vars);
tic
while (iter < maxit)
    iter = iter+1;
    % Check for termination. We have found a sufficiently accurate
    % solution if the natural residual is below tol.
    % Info.NatRes(iter) = r;
    Info.NatRes(iter).primal = norm(nr_res_p);
    Info.NatRes(iter).dual   = norm(nr_res_d);
    Info.NatRes(iter).compl  = mu;
    if printlevel>0
        fprintf('==================\n');
        fprintf('PPM iteration: %4d\n', iter);
        %       fprintf('Natural Residual: %8.2e\n', Info.NatRes(iter));
        fprintf('NR Primal: %8.2e, NR Dual: %8.2e, NR Comp: %8.2e\n',...
            Info.NatRes(iter).primal,Info.NatRes(iter).dual,Info.NatRes(iter).compl);
        fprintf('==================\n');
    end
    if  (norm(nr_res_p)/(max(10,norm(b))) < tol && norm(nr_res_d)/(max(10,norm(c))) < tol &&  mu < tol )
        if printlevel > 0
            fprintf('optimal solution found\n');
        end
        Info.opt    = 1;
        Info.ExIt   = iter-1;
        Info.IPM_It = IPM_Tot_It;
        Info.nnz     = last_nnz;
        break;
    end
    [x,y,z,Info.IPM(iter)] = prox_eval(c,A,A_tr,Q,b,xk,yk,zk,rho,delta,...
        free_variables,pos_vars,num_of_pos_vars,...
        PPM_red^iter,IPM_maxit,pc,printlevel,IterStruct);
    IPM_Tot_It = Info.IPM(iter).IPMIter+IPM_Tot_It;
    last_nnz    = Info.IPM(iter).nnz;
    if Info.IPM(iter).opt == 2
        Info.ExIt   = iter;
        Info.opt    = 2;
        Info.IPM_It = IPM_Tot_It;
        Info.nnz     = last_nnz;
        break;
    end
    %r = natural_PPM_res(c,A,b,Q,x,y,pos_vars);
    [nr_res_p,nr_res_d,mu] = IPM_Res(c,A,A_tr,b,Q,x,y,z,pos_vars, ...
        num_of_pos_vars);
    xk = x;
    yk = y;
    zk = z;
    if  iter == maxit
        Info.opt      = 0;
        Info.ExIt     = iter;
        Info.IPM_It = IPM_Tot_It;
        Info.nnz      = last_nnz;
    end
    % The PPM has terminated either because the solution accuracy is
    % reached or the maximum number of iterations is exceeded. Print
    % result.
end

IPMTT = toc;
if  printlevel>0
    fprintf('time: %g\n', IPMTT);
end
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

% S. Cipolla, J. Gondzio.

function [x,y,z,OInfo] = prox_eval(c,A,A_tr,Q,b,xk,yk,zk,rho,delta,...
    free_variables,pos_vars,num_of_pos_vars,...
    tol,maxit,pc,pl,IterStruct)
% ==================================================================================================================== %
% This function is an Interior Point-Proximal Method of Multipliers, suitable for solving linear and convex quadratic
% programming problems. The method takes as input a problem of the following form:
%
%                                    min   c^T x + (1/2)x^TQx +
%                                    (rho/2)|x-xk|^2 + (delta/2)|y|^2
%                                    s.t.  A x + delta(y-yk) = b,
%                                          x_C >= 0, for i in C \subset {1,...,n},
%                                          x_F free, for i in F = {1,...,n}\C.
%
% INPUT PARAMETERS:
% IP_PMM(c, A, Q, b): find the optimal solution of the problem, with an error tolerance of 10^(-6).
%                     Upon success, the method returns x (primal solution), y (Lagrange multipliers) and
%                     z >= 0 (dual optimal slack variables). If the run was unsuccessful, the method  either returns
%                     a certificate of infeasibility, or terminates after 100 iterations. By default, the method
%                     scales the constraint matrix.
% IP_PMM(c, A, Q, b, free_variables): The last parameter is a matrix of indices, pointing to the free variables of the
%                                     problem. If not given, it is assumed that there are no free variables.
% IP_PMM(c, A, Q, b, free_variables, tol): This way, the user can specify the tolerance to which the problem is solved.
% IP_PMM(c, A, Q, b, free_variables, tol, max_it): This way, the user can also specify the maximum number of iterations.
% IP_PMM(c, A, Q, b, free_variables, tol, maxit, pc): predictor-corrector option.
%     1 = No Predictor-Corrector, 2 Merhotra P-C, 3 Gondzio Multiple Corrections
% IP_PMM(c, A, Q, b, free_variables, tol, max_it,pc, printlevel): sets the printlevel.
%                                                              0: turn off iteration output
%                                                              1: print primal and dual residual and duality measure
%                                                              2: print centering parameter and step length
% OUTPUT: [x,y,z,opt,iter], where:
%         x: primal solution
%         y: Lagrange multiplier vector
%         z: dual slack variables
%         opt: true if problem was solved to optimality, false if problem not solved or found infeasible.
%         iter: numeber of iterations to termination.
%
% Authors: S. Cipolla, J. Gondzio.
% ==================================================================================================================== %
% Output Structure
OInfo = struct();
% ===============================================================
% Auxilliary Initializations
[m,n] = size(A);
x=xk;
y=yk;
z=zk;
e_pos_vars = ones(num_of_pos_vars,1);       % Vector of ones of dimension |C|.
% ==================================================================================================================== %
% Initialize parameters
% -------------------------------------------------------------------------------------------------------------------- %
iter = 0;
alpha_x = 0;     % Step-length for primal variables (initialization)
alpha_z = 0;     % Step-length for dual variables (initialization)
sigmamin = 0.05; % Heuristic value.
sigmamax = 0.95; % Heuristic value.
sigma = 0;
OInfo.opt = 0;

if (num_of_pos_vars > 0)                             % Defined only when non-negativity constraints are present.
    mu = (x(pos_vars)'*z(pos_vars))/num_of_pos_vars; % Initial value of mu.
    res_mu = zeros(n,1);
else
    mu = 0;     % Switch to a pure PMM method (no inequality constraints).
    res_mu = [];
end
header(pl);     % Set the printing choice.

if (pc == 1)
    retry = 0;  % Num of times a factorization is re-built (for different regularization values)
else
    retry_p = 0;
    retry_c = 0;
end
max_tries = 10; % Maximum number of times before exiting with an ill-conditioning message.
mu_prev = 0;
%reg_limit = 1e-1*max(5*tol*(1/max(norm(A,'inf')^2,norm(Q,'inf')^2)),5*10^(-10));
reg_limit  = 1e-1*rho;
% ==================================================================================================================== %
while (iter < maxit)
    % -------------------------------------------------------------------------------------------------------------------- %
    % IP-PMM Main Loop structure:
    % Until (PPM related residual < tol) do
    %   Choose sigma in [sigma_min, sigma_max] and solve:
    %
    % [ (Q + Theta^{-1} + rho I)   -A^T ](Delta x)   -(c + Qx_k - A^Ty_k -[z_k] + rho (x-xk))+ X_k^-1(sigma*mu-Z_kX_k)
    % [         A               delta I ](Delta y) = -(Ax_k + delta (y-yk) - b)
    %
    % with suitable modifications when the constarined set is smaller that the
    % whole ste of variables.
    %   Find two step-lengths a_x, a_z in (0,1] and update:
    %       x_{k+1} = x_k + a_x Delta x, y_{k+1} = y_k + a_z Delta y, z_{k+1} = z_k + a_z Delta z
    %   k = k + 1
    % End
    % -------------------------------------------------------------------------------------------------------------------- %
    if (iter > 1)
        res_p  = new_res_p;
        res_d  = new_res_d;
        res_n =  new_res_n;
    else
        F     = Primal_Dual_Res(x,y,z,pos_vars,c,A,b,Q,mu,xk,yk,rho, delta);
        res_d = -F{1};                   % Regularized dual residual.
        res_p = -F{2};                   % Regularized primal residual.
        res_n = natural_prox_res(c,A,b,Q,x,y,pos_vars,xk,yk,rho,delta);
    end
    % ================================================================================================================ %
    % Check termination criteria
    % ---------------------------------------------------------------------------------------------------------------- %

    if res_n < 10^4*tol*min(1,norm(x-xk)+norm(y-yk))
        OInfo.IPMIter= iter;
        OInfo.opt    = 1;
        OInfo.nnz    = NS.nnz;
        break;
    end
    % ================================================================================================================ %
    iter = iter+1;
    % ================================================================================================================ %
    % ================================================================================================================ %
    % Compute the Newton factorization.
    % ---------------------------------------------------------------------------------------------------------------- %
    pivot_thr = reg_limit;
    %pivot_thr    = 1e-8;
    NS = Newton_factorization(A,A_tr,Q,x,z,delta,rho,pos_vars,free_variables,pivot_thr,IterStruct);
    % ================================================================================================================ %
    switch pc
        case 1 % No predictor-corrector.
            % ============================================================================================================ %
            % Compute the parameter sigma and based on the current solution
            % ------------------------------------------------------------------------------------------------------------ %
            if (iter > 1)
                sigma = max(1-alpha_x,1-alpha_z)^5;
                %sigma  = 1e3* min(x).^2/mu;
            else
                sigma = 0.5;
            end

            sigma = min(sigma,sigmamax);
            sigma = max(sigma,sigmamin);
            % ============================================================================================================ %
            if (num_of_pos_vars > 0)
                res_mu(pos_vars) = (sigma*mu).*e_pos_vars - x(pos_vars).*z(pos_vars);
            end
            % ============================================================================================================ %
            % Solve the Newton system and calculate residuals.
            % ------------------------------------------------------------------------------------------------------------ %
            [dx,dy,dz,instability] = Newton_backsolve(NS,res_p,res_d,res_mu,pos_vars,free_variables,IterStruct);
            if (instability == true) % Checking if the matrix is too ill-conditioned. Mitigate it.
                if (retry < max_tries)
                    fprintf('The system is re-solved, due to bad conditioning.\n')
                    iter = iter -1;
                    retry = retry + 1;
                    reg_limit = min(reg_limit*100,tol);
                    continue;
                else
                    fprintf('The system matrix is too ill-conditioned.\n');
                    OInfo.opt = 2;
                    OInfo.IPMIter = iter;
                    break;
                end
            end
            % ============================================================================================================ %
        case 2 % Mehrotra predictor-corrector. ONLY when num_of_pos_vars > 0!!
            % ================================================================================================================ %
            % Predictor step: Set sigma = 0. Solve the Newton system and compute a centrality measure.
            % ---------------------------------------------------------------------------------------------------------------- %
            res_mu(pos_vars) = - x(pos_vars).*z(pos_vars);
            % ============================================================================================================ %
            % Solve the Newton system with the predictor right hand side -> Optimistic view, solve as if you wanted to
            %                                                               solve the original problem in 1 iteration.
            % ------------------------------------------------------------------------------------------------------------ %
            [dx,dy,dz,instability] = Newton_backsolve(NS,res_p,res_d,res_mu,pos_vars,free_variables,IterStruct);
            if (instability == true) % Checking if the matrix is too ill-conditioned. Mitigate it.
                if (retry_p < max_tries)
                    fprintf('The system is re-solved, due to bad conditioning  of predictor system.\n')
                    iter = iter -1;
                    retry_p = retry_p + 1;
                    reg_limit = min(reg_limit*100,tol);
                    continue;
                else
                    fprintf('The system matrix is too ill-conditioned.\n');
                    OInfo.opt = 2;
                    OInfo.IPMIter = iter;
                    break;
                end
            end
            retry = 0;
            % ============================================================================================================ %

            % ============================================================================================================ %
            % Step in the non-negativity orthant.
            % ------------------------------------------------------------------------------------------------------------ %
            idx = false(n,1);
            idz = false(n,1);
            idx(pos_vars) = dx(pos_vars) < 0; % Select all the negative dx's (dz's respectively)
            idz(pos_vars) = dz(pos_vars) < 0;
            alphamax_x = min([1;-x(idx)./dx(idx)]);
            alphamax_z = min([1;-z(idz)./dz(idz)]);
            tau = 0.995;
            alpha_x = tau*alphamax_x;
            alpha_z = tau*alphamax_z;
            % ============================================================================================================ %
            centrality_measure = (x(pos_vars) + alpha_x.*dx(pos_vars))'*(z(pos_vars) + alpha_z.*dz(pos_vars));
            mu = (centrality_measure/(num_of_pos_vars*mu))^2*(centrality_measure/num_of_pos_vars);
            % ================================================================================================================ %

            % ================================================================================================================ %
            % Corrector step: Solve Newton system with the corrector right hand side. Solve as if you wanted to direct the
            %                 method in the center of the central path.
            % ---------------------------------------------------------------------------------------------------------------- %
            res_mu(pos_vars) = mu.*e_pos_vars - dx(pos_vars).*dz(pos_vars);
            % ============================================================================================================ %
            % Solve the Newton system with the predictor right hand side -> Optimistic view, solve as if you wanted to
            %                                                               solve the original problem in 1 iteration.
            % ------------------------------------------------------------------------------------------------------------ %
            [dx_c,dy_c,dz_c,instability] = Newton_backsolve(NS,zeros(m,1),zeros(n,1),res_mu,pos_vars,free_variables,IterStruct);
            if (instability == true) % Checking if the matrix is too ill-conditioned. Mitigate it.
                if (retry_c < max_tries)
                    fprintf('The system is re-solved, due to bad conditioning of corrector.\n')
                    iter = iter -1;
                    retry_c = retry_c + 1;
                    mu = mu_prev;
                    reg_limit = min(reg_limit*100,tol);
                    continue;
                else
                    fprintf('The system matrix is too ill-conditioned, increase regularization.\n');
                    OInfo.opt = 2;
                    OInfo.IPMIter = iter;
                    break;
                end
            end
            % ============================================================================================================ %
            dx = dx + dx_c;
            dy = dy + dy_c;
            dz = dz + dz_c;

    end
    % ================================================================================================================ %
    % Compute the new iterate:
    % Determine primal and dual step length. Calculate "step to the boundary" alphamax_x and alphamax_z.
    % Then choose 0 < tau < 1 heuristically, and set step length = tau * step to the boundary.
    % ---------------------------------------------------------------------------------------------------------------- %
    if (num_of_pos_vars > 0)
        idx = false(n,1);
        idz = false(n,1);
        idx(pos_vars) = dx(pos_vars) < 0; % Select all the negative dx's (dz's respectively)
        idz(pos_vars) = dz(pos_vars) < 0;
        alphamax_x = min([1;-x(idx)./dx(idx)]);
        alphamax_z = min([1;-z(idz)./dz(idz)]);
        tau  = 0.995;
        alpha_x = tau*alphamax_x;
        alpha_z = tau*alphamax_z;
    else
        alpha_x = 1;         % If we have no inequality constraints, Newton method is exact -> Take full step.
        alpha_z = 1;
    end
    % ================================================================================================================ %

    % ================================================================================================================ %
    % Make the step.
    % ---------------------------------------------------------------------------------------------------------------- %
    x = x+alpha_x.*dx; y = y+alpha_z.*dy; z = z+alpha_z.*dz;
    if (num_of_pos_vars > 0) % Only if we have non-negativity constraints.
        mu_prev = mu;
        mu = (x(pos_vars)'*z(pos_vars))/num_of_pos_vars;
        % mu_rate = abs((mu-mu_prev)/max(mu,mu_prev));
    end
    % ================================================================================================================ %
    % Computing the new residuals.
    % ================================================================================================================ %
    F     = Primal_Dual_Res(x,y,z,pos_vars,c,A,b,Q,mu,xk,yk,rho, delta);
    new_res_d = -F{1};                   % Regularized dual residual.
    new_res_p = -F{2};                   % Regularized primal residual.
    new_res_n = natural_prox_res(c,A,b,Q,x,y,pos_vars,xk,yk,rho,delta);
    % ================================================================================================================ %
    % Print iteration output.
    % ---------------------------------------------------------------------------------------------------------------- %
    pres_inf = norm(new_res_p);
    dres_inf = norm(new_res_d);
    output(pl,iter,pres_inf,dres_inf,mu,sigma,alpha_x,alpha_z);
    % ================================================================================================================ %
end % while (iter < maxit)

if iter == maxit
    OInfo.IPMIter=maxit;
    OInfo.opt = 0;
end


% The IPM has terminated because the solution accuracy is reached or the maximum number
% of iterations is exceeded. Print result.

if (pl >0 )
    fprintf('iterations: %4d\n', iter);
    fprintf('primal feasibility: %8.2e\n', norm(res_p));
    fprintf('dual feasibility: %8.2e\n', norm(res_d));
    fprintf('complementarity: %8.2e\n', full(dot(x,z)/n));
end
end


% ==================================================================================================================== %
% header + output printing functions:
% pl = 1: primal-dual infeasibility and mu is printed at each iteration k
% pl = 2: primal-dual infeasibility, mu, sigma, and step-lengths are printed at each iteration k
% -------------------------------------------------------------------------------------------------------------------- %
function header(pl)
if (pl >= 1)
    fprintf(' ');
    fprintf('%4s    ', 'iter');
    fprintf('%8s  ', 'pr feas');
    fprintf('%8s  ', 'dl feas');
    fprintf('%8s  ', 'mu');
end
if (pl >= 2)
    fprintf('  ');
    fprintf('%8s  ', 'sigma');
    fprintf('%8s  ', 'alpha_x');
    fprintf('%8s  ', 'alpha_z');
end
if (pl >= 1)
    fprintf('\n ====    ========  ========  ========');
end
if (pl >= 2)
    fprintf('    ========  ========  ========');
end
if (pl >= 1); fprintf('\n'); end
end


function output(pl,it,xinf,sinf,mu,sigma,alpha_x,alpha_z)
if (pl >= 1)
    fprintf(' ');
    fprintf('%4d    ', it);
    fprintf('%8.2e  ', xinf);
    fprintf('%8.2e  ', sinf);
    fprintf('%8.2e  ', mu);
end
if (pl >= 2)
    fprintf('  ');
    fprintf('%8.2e  ', sigma);
    fprintf('%8.2e  ', alpha_x);
    fprintf('%8.2e  ', alpha_z);
end
if (pl >= 1); fprintf('\n'); end
end

function [F] = Primal_Dual_Res(x,y,z,~,...
    g,A,b,H,~,...
    xk,yk,rho, delta)
F    = cell(2,1);
F{1} = H*x+g-A.'*y-z+rho*(x-xk); % xi_d
F{2} = A*x-b+delta.*(y-yk);      % xi_p
%F{3} = z(pos_vars).*x(posvars)-mu;
end

function [npres] = natural_prox_res(g,A,b,H,x,y,pos_var,xk,yk,rho,delta)
r     = cell(2,1);
r{1}  = H*x+g-A.'*y+rho.*(x-xk);
r{1}(pos_var) = x(pos_var)-max(0,x(pos_var)-r{1}(pos_var));
r{2}  = A*x-b+delta.*(y-yk);
npres = norm(r{1})+norm(r{2});
end

function NS = Newton_factorization(A,A_tr,Q,x,z,delta,rho,pos_vars,free_vars,pivot_thr,Struct)
% ==================================================================================================================== %
% Newton_factorization: Factorize the Newton matrix
% -------------------------------------------------------------------------------------------------------------------- %
% NS = Newton_factorization(A,A_tr,Q,x,z,delta,rho,pos_vars,free_vars) returns a MATLAB struct that holds the
%      factorization of the Newton matrix for solving the step equations in
%      the IPM, as well as relevant information concerning failure.
% Factorization Method
% --------------------
% 1: augmented system, LDL' factorization.
%
% Author: S.Cipolla, J. Gondzio.
% ==================================================================================================================== %
[m, n] = size(A);
NS = struct();
% ==================================================================================================================== %
% LDL' factorization of KKT matrix
% -------------------------------------------------------------------------------------------------------------------- %
%
% MATLAB uses MA57 when K is sparse, which is not available in OCTAVE.
% -------------------------------------------------------------------------------------------------------------------- %
NS.x = x;
NS.z = z;
Q_bar = zeros(n,1);
if (size(pos_vars,1) > 0)
    Q_bar(pos_vars) = z(pos_vars)./x(pos_vars) + rho;
    Q_bar(free_vars) = rho;
else
    Q_bar(:) = rho;
end

if strcmp(Struct.Fact,'ldl')
    K = [Q+spdiags(Q_bar,0,n,n), A_tr; A, -spdiags(delta.*ones(m,1),0,m,m)];
    %condest(K)
    [NS.L,NS.D,NS.pp] = ldl(K,pivot_thr,'vector'); %Small pivots allowed, to avoid 2x2 pivots.
    NS.nnz   =  nnz(NS.L);
else
    NS.A = A;
    NS.D =  Q+spdiags(Q_bar,0,n,n); %  1./Q_bar;
    K       = A*(NS.D\A_tr)+spdiags(delta.*ones(m,1),0,m,m);
    [NS.L,~,NS.pp] = chol(K, 'vector');
    NS.nnz   =  nnz(NS.L);
end

% ==================================================================================================================== %

% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
end


function [dx,dy,dz,instability] = Newton_backsolve(NS,res_p,res_d,res_mu,pos_vars,free_vars, Struct)
% ==================================================================================================================== %
% Newton_backsolve    Solve linear system with factorized matrix, by using backward substitution.
% -------------------------------------------------------------------------------------------------------------------- %
% OUTPUT:
%  [dx,dy,dz,instability] = newtonsolve(NS,res_p,res_d,res_mu,A,A_tr,pos_vars,free_vars)
%  i.e. the Newton direction and a boolean parameter indicating critical ill-conditioning.
%
% Author: S. Cipolla, J. Gondzio.
% ==================================================================================================================== %
m = size(res_p,1);
n = size(res_d,1);
instability = false;
dx = zeros(n,1);
dz = zeros(n,1);
dy = zeros(m,1);
if (size(pos_vars,1) > 0)
    temp_res = zeros(n,1);
end
% ==================================================================================================================== %
% Solve KKT system with LDL' factors.
% -------------------------------------------------------------------------------------------------------------------- %
if strcmp(Struct.Fact,'ldl')
    if (size(pos_vars,1) > 0)
        temp_res(pos_vars) =  res_mu(pos_vars)./(NS.x(pos_vars));
        rhs = [res_d+temp_res; res_p];
    else
        rhs = [res_d; res_p];
    end
    warn_stat = warning;
    warning('off','all');
    lhs = NS.L'\(NS.D\(NS.L\rhs(NS.pp)));
    if (nnz(isnan(lhs)) > 0 || nnz(isinf(lhs)) > 0)
        instability = true;
        return;
    end
    warning(warn_stat);
    lhs(NS.pp) = lhs;
    dx = lhs(1:n,1);
    dy = -lhs(n+1:n+m,1);
else % Chol
    if (size(pos_vars,1) > 0)
        temp_res(pos_vars) =  res_mu(pos_vars)./(NS.x(pos_vars));
        res_d = res_d+temp_res;
    end
    warn_stat = warning;
    %warning('off','all');
    rhs = res_p-NS.A*(NS.D\res_d);
    lhs  = NS.L\(NS.L.'\(rhs(NS.pp)));
    dy(NS.pp)  = lhs;
    if (nnz(isnan(dy)) > 0 || nnz(isinf(dy)) > 0)
        instability = true;
        return;
    end
    warning(warn_stat);
    dx = NS.D\(NS.A.'*dy+res_d);
end

if (size(pos_vars,1) > 0)
    dz(pos_vars) = (res_mu(pos_vars)-NS.z(pos_vars).*dx(pos_vars))./NS.x(pos_vars);
    dz(free_vars) = 0;
end

% ==================================================================================================================== %

% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
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
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %
