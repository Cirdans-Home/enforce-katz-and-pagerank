function [Phat,Delta,fronorm] = enforce_pagerank_closedform(A,alpha,v,pihat)
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

% Compute the original PageRank Matrix
n = size(A,1);
e = ones(n,1);
I = speye(n,n);
D = spdiags(sum(A,2),0,n,n);
P = D\A;
G = alpha*P + (1-alpha)*e*v';


% Compute the steady state vector of P
[w,~] = eigs(P',1,'largestabs');
w = w/sum(w);

% Compute the \sigma_* as in Proposition 4.1
M = ( I - P + e*w')';
sig_star =  M\((1-alpha)*v - (I - alpha*P)'*pihat);
sig_star = sig_star./(alpha*pihat);

% Compute the sigma satisfying the constraint sigma > -1
if any(sig_star < -1)
    gam = max(max((-1 - sig_star).*pihat./w));
    sig_star = sig_star + (gam*w)./pihat;  
end
% This is the Delta perturbation from Proposition 4.1
Delta = spdiags(sig_star,0,n,n)*(A - D);

% Build the Phat matrix as in Proposition 4.2
rhat = max(1 - alpha*diag(D\(A+Delta)));
alphahat = 1 - (1-alpha)/rhat;
Phat = 1/(rhat-1+alpha)*( alpha*(D\(A+Delta)) + (rhat - 1)*I);
Ghat = alphahat*Phat + (1-alphahat)*e*v';

fronorm = norm(G-Ghat,"fro");


end