%% Enforcing Pagerank

clear; clc; close all;

load("../testmatrices/uk.mat");
A       = Problem.A;
alpha = 0.8;
n        = size(A,1);
I         = speye(n,n);
e        = ones(n,1);
v        = (1/n).*e;  % Teleportation
deg    = A*e;

% Check Reducibility of the associated Markov Chain 
%M       = spdiags(1./deg,0,n,n)*A;
%isreducible(dtmc(M))

mu      = (I - alpha*(spdiags(1./deg,0,n,n)*A).')\((1-alpha).*v) ;

muhat                    = mu;
[~, max_ind]           = maxk(mu,10);
[~, min_ind]            = mink(mu,10);
muhat(max_ind) = mu(min_ind);
muhat(min_ind)  = mu(max_ind);
   


beta    = 0.08;
P         = spones(A+spdiags(e,0,n,n));
[Delta,stat,mucheck,rhat] = enforce_pagerank(A,alpha,muhat,v,P,beta,1e-8);
fprintf('Objective: %e\n',beta*norm(Delta,"fro")^2 + (1-beta)*norm(Delta,1));

figure(1)
plot(1:n,mu,'--',1:n,muhat,'o',1:n,mucheck,'x','LineWidth',2);

figure(2)
semilogy(0:stat.IPMiter,stat.primalres,0:stat.IPMiter,stat.dualres, ...
    0:stat.IPMiter,stat.compl,'LineWidth',2);
xlabel('IPM Iteration')
legend('Primal Residual','Dual Residual','Complementarity Residual');

figure(3)
%Delta(abs(Delta)<1e-10) = 0;
%Delta = sparse(Delta);
subplot(1,3,1)
spy(A)
subplot(1,3,2)
spy(Delta > 1e-10,'b+')
subplot(1,3,3)
spy(Delta < -1e-10,'r-')
