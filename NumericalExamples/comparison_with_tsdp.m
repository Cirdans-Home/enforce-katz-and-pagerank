%% Comparison of our method for PR enforcing with the TSDP methods

clear; clc; close all;

addpath('../TSDP/');        % Gillis and Van Dooren Code
addpath('../enforcers/')    % Our PageRank Enforcer

 set(0, 'DefaultAxesFontSize', 14);

%% Load Test Problem
load('../matrices/karate.mat');
A = Problem.A;
n = size(A,1);
e = ones(n,1);
I = speye(n,1);
v = ones(n,1);
v = v/norm(v,1);

%% Compute the true PageRank vector
alpha = 0.85;
D = spdiags(sum(A,2),0,n,n);
P = D\A;
G = alpha*P + (1-alpha)*e*v';
[pi,lam] = eigs(G',1,'largestabs');
pi = pi/sum(pi);
figure(1)
plot(1:n,pi,'b--','LineWidth',2)
xlim([1 n])
title('PageRank')

%% What is the PageRank that we want?
pihat = pi;
pihat(1) = 0.5*pi(34);
pihat(34) = 0.5*pi(34);
pihat = pihat/sum(pihat);
figure(1)
plot(1:n,pi,'b--',1:n,pihat,'LineWidth',2)
xlim([1 n])
title('PageRank')
legend({'$\pi$','$\hat{\pi}$'},'Interpreter','latex')

%% Use the IPM to obtain it 
P = spones(A+I);
beta = 0.1;
tol = 1e-12;
[Delta,info,picheck,rhat,Phat] = enforce_pagerank(A,alpha,pihat,v,P,beta,tol);
Delta(abs(Delta)<1e-10) = 0;
figure(1)
subplot(3,1,1)
plot(1:n,pi,'b--',1:n,pihat,1:n,picheck,'ko','LineWidth',2)
xlim([1 n])
title('PageRank IPM')
legend({'$\pi$','$\hat{\pi}$','Optimized $\hat{\pi}$ IPM'},'Interpreter','latex','Location','eastoutside')

%% Use the GUROBI routine to obtain it
beta = 0.2;
[Delta_gu,info_gu,picheck_gu,rhat_gu,Phat_gu] = enforce_pagerank_gurobi(A,alpha,pihat,v,P,beta);
Delta_gu(abs(Delta_gu)<1e-10) = 0;
figure(1)
subplot(3,1,2)
plot(1:n,pi,'b--',1:n,pihat,1:n,picheck_gu,'ko','LineWidth',2)
xlim([1 n])
title('PageRank IPM')
legend({'$\pi$','$\hat{\pi}$','Optimized $\hat{\pi}$ GUROBI'},'Interpreter','latex','Location','eastoutside')
subplot(3,1,3)
semilogy(1:n,abs(pihat-picheck),'-',1:n,abs(pihat-picheck_gu),'--','LineWidth',2)
title('Absolute Error')
xlim([1 n])
legend({'IPM','GUROBI'},'Location','eastoutside')

%% Use the TDSP method (No target)
[D_not,results_not,model_not,timings_not] = supportTSDP(G,pihat);
[pitdsp,~] = eigs((G+D_not)',1,'largestabs');
pitdsp = pitdsp/sum(pitdsp);

%% Use the TDSP methd (sparsity of A target)
options.support = spones(P+I);
[D_tar,results_tar,model_tar,timings_tar] = supportTSDP(G,pihat,options);
[pitdsptar,~] = eigs((G+D_tar)',1,'largestabs');
pitdsptar = pitdsptar/sum(pitdsptar);

%% Compare results
figure(1)
subplot(1,4,1)
spy(Delta)
title('IPM')
subplot(1,4,2)
spy(Delta_gu)
title('Gurobi')
subplot(1,4,3)
spy(D_not)
title('TDSP (Dense Matrix)')
subplot(1,4,4)
spy(D_tar)
title('TDSP Pattern of A + I')

figure(2)
semilogy(1:n,abs(pihat-picheck),'-', ...
    1:n,abs(pihat-picheck_gu),'--', ...
    1:n,abs(pihat-pitdsptar),'o-',...
    1:n,abs(pihat-pitdsp),'^--',...
    'LineWidth',2)
legend({'IPM','GUROBI','TDSP Pattern of A + I','TDSP (Dense Matrix)'},...
    'Location','eastoutside')

figure(3)
fronorms = [norm(Delta,'fro'),norm(Phat-G,'fro'),norm(Delta_gu,'fro')...
    norm(Phat_gu-G,'fro'),norm(D_not,'fro'),norm(D_tar,'fro')];
X = ["||\Delta||_F","||G-\hat{P}||_F","||\Delta_{GU}||_F",...
    "||G-\hat{P}_{GU}||_F","||D||_F","||D_2||_F"];
bar(X,fronorms)


%% Does the TDSP preserve the PageRank structure?
G_tsdp_not = G + D_not;
P_tsdp_not = (G_tsdp_not - (1-alpha)*e*v')/alpha;
figure(4)
subplot(1,2,1)
hold on
spy(P_tsdp_not > 0)
spy(P_tsdp_not < 0,'rx')
hold off

%% And in the other
G_tsdp_tar = G + D_tar;
P_tsdp_tar = (G_tsdp_tar - (1-alpha)*e*v')/alpha;
subplot(1,2,2)
hold on
spy(P_tsdp_tar > 0)
spy(P_tsdp_tar < 0,'rx')
hold off
