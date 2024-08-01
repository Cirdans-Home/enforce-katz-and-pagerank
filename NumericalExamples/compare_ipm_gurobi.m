%% Comparison of IPM and GUROBI Optimization for PageRank
% This test compares the performances and results of the PageRank enforcers
% using either the IPM or the GUROBI optimization routines.

clear; clc; close all;
warning('off','all')
addpath('../enforcers/')
addpath('../testmatrices/')
addpath('/software/gurobi/gurobi1102/linux64/matlab');
gurobi_setup

% 2 6 13 15
% PGPgiantcompo ct2010 nh2010 vt2010
fid = fopen('compare_ipm_gurobi.txt','a+');
test_matrices = ["PGPgiantcompo","ct2010","nh2010","vt2010"];
betavals = [1,0.50,1/101,1/201];

fprintf(fid,'\n\n');
[info,machine] = system('hostname');
fprintf(fid,'Test %s machine %s',string(datetime("today")),machine);

ntests = length(test_matrices);
nbetas = length(betavals);

% Data collection array
TIME         = NaN(ntests,nbetas,2);
TAU          = NaN(ntests,nbetas,2);
RHAT         = NaN(ntests,nbetas,2);
RELDELTANORM = NaN(ntests,nbetas,2);
ABSDELTANORM = NaN(ntests,nbetas,2);
NNZ          = NaN(ntests,nbetas,2);
ITER         = NaN(ntests,nbetas,2);

fprintf(fid,['& & \\multicolumn{7}{c}{IPM} & \\multicolumn{7}{c}{GUROBI}' ...
    ' \\\\\n']);
fprintf(fid,['Matrix & $\\beta$ & Time & $\\tau$ & $\\hat{r}$ & Rel. & Abs. & nnz & ' ...
    'Iter. & $\\tau$ & $\\hat{r}$ & Rel. & Abs. & nnz & Iter. \\\\\n']);
for i = 1:length(test_matrices)
    % Load the matrix
    load(sprintf("%s.mat",test_matrices(i)));
    A = Problem.A;
    % Build the original PageRank
    n = size(Problem.A,1);
    I        = speye(n,n);
    e        = ones(n,1);
    alpha    = 0.9;
    v        = (1/n).*e;  % Teleportation
    deg      = Problem.A*e;
    mu       = (I - alpha*(spdiags(1./deg,0,n,n)*Problem.A).')\((1-alpha).*v) ;
    % Normalize it
    mu = mu./sum(mu);
    % Locate the top-10%
    [~,bestmu] = sort(mu,"descend");
    muhat                   = mu;
    nindex                  = round(0.1*n);
    muhat(bestmu(1:nindex)) = mean(mu(bestmu(1:nindex)));
    P                       = spones(A+I);
    [~,ind_muhat] = sort(muhat,"descend");
    normA = norm(A,"fro");
    fprintf(fid,'%s ',test_matrices(i));
    for j = 1:length(betavals)
        % Enforce with IPM
        tol  = 1e-9;
        beta = betavals(j);
        fprintf(fid,'& %1.4f & ',beta);
        [Delta_IPM,Info_IPM,picheck_IPM,rhat_IPM] = ...
            enforce_pagerank(A,alpha,muhat,v,P,beta,tol);
        [~,ind_ipm] = sort(picheck_IPM,"descend");
        % Enforce with GUROBI
        [Delta_GU,Info_GU,picheck_GU,rhat_GU] = ...
            enforce_pagerank_gurobi(A,alpha,muhat,v,P,beta);
        [~,ind_gu] = sort(picheck_GU,"descend");
        % Populate data
        TIME(i,j,1)         = Info_IPM.time;
        TAU(i,j,1)          = corr(picheck_IPM,muhat,"type","Kendall");
        RHAT(i,j,1)         = rhat_IPM;
        ABSDELTANORM(i,j,1) = norm(Delta_IPM,"fro");
        RELDELTANORM(i,j,1) = ABSDELTANORM(i,j,1)/normA;        
        NNZ(i,j,1)          = nnz(Delta_IPM);
        ITER(i,j,1)         = Info_IPM.IPMiter;
        
        TIME(i,j,2)         = Info_GU.time;
        TAU(i,j,2)          = corr(picheck_GU,muhat,"type","Kendall");
        RHAT(i,j,2)         = rhat_GU;
        ABSDELTANORM(i,j,2) = norm(Delta_GU,"fro");
        RELDELTANORM(i,j,2) = ABSDELTANORM(i,j,2)/normA;
        NNZ(i,j,2)          = nnz(Delta_GU);
        ITER(i,j,2)         = Info_GU.iter + Info_GU.baritercount;

        % Save data to table
        fprintf(fid,['%1.2e & %1.2f & %1.2f & %1.2e & %1.2e & %d & %d & ' ...
            '%1.2e & %1.2f & %1.2f & %1.2e & %1.2e & %d & %d'],...
            TIME(i,j,1),...
            TAU(i,j,1),...
            RHAT(i,j,1),...
            ABSDELTANORM(i,j,1),...
            RELDELTANORM(i,j,1),...        
            NNZ(i,j,1),...
            ITER(i,j,1),...
            TIME(i,j,2),...
            TAU(i,j,2),...
            RHAT(i,j,2),...
            ABSDELTANORM(i,j,2),...
            RELDELTANORM(i,j,2),...        
            NNZ(i,j,2),...
            ITER(i,j,2));
        fprintf(fid,'\\\\\n');
    end
end

fclose(fid);
save("compare_ipm_gurobi2.mat","ITER","NNZ","RELDELTANORM",...
    "ABSDELTANORM","RHAT","TAU","TIME");
