%% Enforcing Katz

clear; clc; close all;

load("../testmatrices/uk.mat");
A       = Problem.A;
rhoA    = eigs(A,1,"largestabs");
alpha   = 1/(2*rhoA);
n       = size(A,1);
I       = speye(n,n);
e       = ones(n,1);

mu      = (I - alpha*A)\e;

[~,bestmu] = sort(mu,"descend");

perc = 0.1;
muhat   = mu;
muhat(3) = mu(1);
muhat(1) = mu(3);
% nindex  = round(perc*n);
% muhat(bestmu(1:nindex)) = mean(mu(bestmu(1:nindex)));


beta    = 0.1;
P       = spones(A);
[Delta,stat,mucheck] = enforce_katz(A,alpha,muhat,P,beta,1e-12);
fprintf('Objective: %e\n',beta*norm(Delta,"fro")^2 + (1-beta)*norm(Delta,1));

figure(1)
plot(1:n,mu,'--',1:n,muhat,'o',1:n,mucheck,'x','LineWidth',2);

figure(2)
semilogy(0:stat.iter,stat.primalres,0:stat.iter,stat.dualres, ...
    0:stat.iter,stat.compl,'LineWidth',2);
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
