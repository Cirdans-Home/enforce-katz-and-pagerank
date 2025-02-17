%% Test with higher precision

clear; clc; close all;

addpath('/opt/matlabauxiliaries/advanpix-4.8.0/')

mp.Digits(64); % Set quadruple precision
addpath("../enforcers/"); % Add the folder containing the enforcers
load('../matrices/karate.mat');

tolvec = mp([1e-9,1e-12,1e-15,1e-18,1e-21,1e-24,1e-27,1e-30,1e-33,1e-36,1e-39,1e-42,1e-45]);
nnzvec = [];
normvec = [];
relnormvec = [];

nnzvec_std = [];
normvec_std = [];
relnormvec_std = [];

A = mp(abs(Problem.A));
N = size(A,1);
I = mp(speye(N,N));
e = mp(ones(N,1));
%% Katz
% Computation of the "original" Katz centrality
rhoA = max(eigs(A));
rhoA = abs(rhoA);
alpha = mp('0.5')/rhoA;
mu = (I - alpha*A)\e;

%% Target problem
% Selection of the target centrality
muhat = mu;
muhat(33) = mp(0.5)*mu(34);
muhat(5) = mp(1.5)*mu(5);

%% Selection of the pattern
% Pattern to be used
P = mp(spones(A));

beta = mp(0.5);
for tol = tolvec
    fprintf("Testing tolerance %e.\n",tol)

    if tol >= mp(1e-15)
        fprintf("\n\nRunning classical:\n\n")
        [Delta_origin,stat_origin,mucheck_origin] = ....
            enforce_katz(double(A),double(alpha),double(muhat),...
            double(P),double(beta),double(tol));
        nnzvec_std = [nnzvec_std,nnz(Delta_origin(abs(Delta_origin) > 1e-10))];
        normvec_std = [normvec_std,norm(double(muhat)-mucheck_origin)];
        relnormvec_std = [relnormvec_std,norm(double(muhat)-mucheck_origin)/norm(double(muhat))];
    end

    fprintf("\n\nRunning mp\n\n")
    [Delta,stat,mucheck] = enforce_katz_hp(A,alpha,muhat,P,beta,tol);

    nnzvec  = [nnzvec,nnz(Delta(abs(Delta) > 1e-10))];
    normvec = [normvec,norm(muhat-mucheck)];
    relnormvec = [relnormvec,norm(muhat-mucheck)/norm(muhat)];
end

%% Visualize results
figure(1)
semilogy(1:length(normvec),normvec,...
    1:length(relnormvec),relnormvec,...
    1:length(normvec_std),normvec_std,...
    1:length(normvec),tolvec(1:length(normvec)),'--','LineWidth',2)
ylabel('Error on Katz centrality vector')
legend({'Augmented precision (abs)','Augmented precision (rel)',...
    'Double precision','Tolerance'},'Location','westoutside')
xticks(1:length(relnormvec))
xticklabels(split(sprintf("%1.1e\n",tolvec(1:length(relnormvec)))))