function [rbo_ext,varargout] = rbosimilarity(mu,muhat,p,varargin)
%RBOSIMILARITY Returns RBO indefinite rank similarity metric, as described 
% in: Webber, W., Moffat, A., & Zobel, J. (2010). A similarity measure for 
% indefinite rankings. ACM Transactions on Information Systems. 
% doi:10.1145/1852102.1852106.

% Input check
narginchk(3,4)
n = length(mu);
if n ~= length(muhat)
    error("RBO: scores of different length")
end
% Output check
nargoutchk(1,2);

scores = transpose(1:n);

% Compute the ranking
[~,rank_mu] = sort(mu,"descend"); 
[~,perm_mu] = sort(rank_mu,"ascend");
rank_mu = scores(perm_mu);
[~,rank_muhat] = sort(muhat,1,"descend");
[~,perm_muhat] = sort(rank_muhat,"ascend");
rank_muhat = scores(perm_muhat);

% Convert MATLAB arrays to Python lists
l1_py = py.list(rank_mu);
l2_py = py.list(rank_muhat);

% Call the Python function calc_rbo
rbo_ext_py = pyrunfile("calc_rbo.py","z",l1=l1_py,l2=l2_py,p=p);

% Convert the result from Python to MATLAB
rbo_ext = double(rbo_ext_py);

% If the second argument is asked for it computes a table depicting the
% results
if nargout == 2
    if nargin == 4
        name = varargin{4};
    else
        name = string(num2str(transpose(1:n)));
    end
    varargout{1} = table(mu,rank_mu,muhat,rank_muhat, ...
        'RowNames',name);
end

end

