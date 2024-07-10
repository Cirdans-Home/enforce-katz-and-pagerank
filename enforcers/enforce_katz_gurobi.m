function [Delta,varargout] = enforce_katz_gurobi(A,alpha,muhat,P,beta)
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
%   NOTE: This version uses the Gurobi solver for the optimization problem:
%   https://www.gurobi.com, it is possible to get an Academic License for
%   it. If you don't want to use, the enforce_katz implementatio is written
%   in pure Matlab and doesn't require external tools.


if nargout >= 2
    varargout{1} = struct();
end

% General infos
N = size(A,1);

% Build the projector
proj         = pattern_projector(P);  % Projector Onto the Pattern of A
reduced_size = size(proj,1);

if beta == 1
    % Solve only with the Frobenius norm constraint
    % We build here the Gurobi model structure:
    % See: https://docs.gurobi.com/projects/optimizer/en/current/reference/matlab.html#matlabproblem
    % min x'Qx -2*c
    % s.t. A x = b
    model.modelname        = 'Katz-Frobenius';
    model.Q                = speye(reduced_size);
    c = proj*reshape(A,N*N,1);
    model.obj              = -2.0.*full(c);
    model.A                = kron(muhat.',speye(N))*proj.';
    model.rhs              = (1/alpha).*(muhat-1)-A*muhat+model.A*c;
    model.sense            = repmat('=',size(model.A,1),1);
    model.objcon           = 0.0;
    model.lb               = zeros(reduced_size,1);
    % Call GUROBI on the problem to solve
    result = gurobi(model);
    % Collect data if the user ask for them
    if nargout >= 2
        varargout{1}.time       = result.runtime;
        varargout{1}.iter       = result.itercount;
        varargout{1}.baritercount = result.baritercount;
    end
    % Recover Matrix and Desired Ranking
    [ival,jval,~] = find(P);
    xfinalvec    = result.x - c;
    Delta = sparse(ival,jval,xfinalvec);
    % Compute the optimizate Katz centrality
    if nargout == 3
        I = speye(N,N);
        e = ones(N,1);
        varargout{2} = (I - alpha*(A+Delta))\e;
    end

elseif beta < 1
    % Solve also with the sparsity constraint
    tau                    = (1-beta)/beta;
    model.modelname        = 'Katz-General';
    Q                      = speye(reduced_size);
    c                      = (proj*reshape(A,N*N,1));
    L                      = kron(muhat.',speye(N))*proj.';
    model.Q                = blkdiag(Q,sparse(reduced_size,reduced_size),sparse(reduced_size,reduced_size));
    model.obj              = [-2.*full(c); tau.*ones(reduced_size,1); tau.*ones(reduced_size,1)];
    model.A                = [L,                     sparse(N,reduced_size), sparse(N,reduced_size);...
                              - speye(reduced_size), speye(reduced_size),      -speye(reduced_size)];
    model.rhs              = full([(1/alpha).*(muhat-1)-A*muhat+L*c; -c]);
    model.sense            = repmat('=',size(model.A,1),1);
    model.objcon           = 0.0;
    model.lb               = zeros(3*reduced_size,1);
    % Call GUROBI on the problem to solve
    result = gurobi(model);
    % Collect data if the user ask for them
    if nargout >= 2
        varargout{1}.time       = result.runtime;
        varargout{1}.iter       = result.itercount;
        varargout{1}.baritercount = result.baritercount;
    end
    % Recover Matrix and Desired Ranking
    [ival,jval,~] = find(P);
    xfinalvec    = result.x(1:reduced_size) - c;
    Delta = sparse(ival,jval,xfinalvec);
    % Compute the optimizate Katz centrality
    if nargout == 3
        I = speye(N,N);
        e = ones(N,1);
        varargout{2} = (I - alpha*(A+Delta))\e;
    end
end



end

function [proj] = pattern_projector(P)
%PATTERN_PROJECTOR given a sparse matrix P this function returns the
%projector onto its sparsity pattern
k     = find(P);
N    = size(P,1);
proj = sparse(1:1:length(k),k,ones(length(k),1), length(k),N*N );
end

