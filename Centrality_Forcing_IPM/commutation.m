% Commutation Matrix
function [K] = commutation(n)
I = reshape(1:n*n, [n, n]); % initialize a matrix of indices of size(A)
I = I';    % Transpose it
I = I(:); % vectorize the required indices
K = speye(n*n); % Initialize an identity matrix
K = K(I,:); % Re-arrange the rows of the identity matrix
end
