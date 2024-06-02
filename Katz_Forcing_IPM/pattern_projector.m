%% Projector on the sparsity pattern for reduced form
function [proj] = pattern_projector(P)
k     = find(P);
N    = size(P,1);
proj = sparse([1:1:length(k)],k,ones(length(k),1), length(k),N*N );
end