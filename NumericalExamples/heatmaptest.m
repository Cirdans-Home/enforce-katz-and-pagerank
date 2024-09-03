%% Testing on the average fixing

clear; clc; close all;

QP_problems_path = "../testmatrices/";
d = dir(fullfile(QP_problems_path,'*.mat')); 
fprintf("Test matrices are: \n")
fprintf("\t%s\n",d.name)

perc = [0.1,0.2,0.3,0.4,0.5];
betav = [0.8,0.6,0.4,0.2,1];
ndigits = 6;

Pert = NaN(length(d),length(perc),length(betav)); % Value of the Objective
Time = NaN(length(d),length(perc),length(betav)); % Time To Solution
Iter = NaN(length(d),length(perc),length(betav)); % IPM Iter to Solution
NNZ  = NaN(length(d),length(perc),length(betav)); % Number of nonzeros in Delta
NNZp = NaN(length(d),length(perc),length(betav)); % Number of pos nonzeros in Delta
NNZn = NaN(length(d),length(perc),length(betav)); % Number of neg nonzeros in Delta
KTAU = NaN(length(d),length(perc),length(betav)); % Number of neg nonzeros in Delta

for i=1:length(d)
    load(fullfile(QP_problems_path,d(i).name));
    
    A = Problem.A;
    P = spones(A);

    % Compute the true mu
    n = size(A,1);
    I = speye(n,n);
    e = ones(n,1);
    rhoA = eigs(A,1,"largestabs");
    alpha = 1/(2*rhoA);
    mu      = (I - alpha*A)\e;
    [~,bestmu] = sort(mu,"descend");

    for k = 1:length(betav)
        for j = 1:length(perc)
            muhat = mu;
            % Compute the changed mu
            nindex  = round(perc(j)*n);
            muhat(bestmu(1:nindex)) = mean(mu(bestmu(1:nindex)));
            beta = betav(k);
            [Delta,stat,mucheck] = enforce_katz(A,alpha,muhat,P,beta,1e-12);

            % Collect Data
            % Value of the Objective
            Pert(i,j,k) = beta*norm(Delta,"fro")^2 + (1-beta)*norm(Delta,1);
            Time(i,j,k) = stat.time; % Time To Solution
            Iter(i,j,k) = stat.IPMiter; % IPM Iter to Solution
            NNZ(i,j,k)  = nnz(Delta); % Number of nonzeros in Delta
            NNZp(i,j,k) = nnz(Delta > 0); % Number of pos nonzeros in Delta
            NNZn(i,j,k) = nnz(Delta < 0); % Number of neg nonzeros in Delta
            KTAU(i,j,k) = corr(round(mucheck,ndigits),round(muhat,ndigits),'type','Kendall');
        end
    end

end

%% Save Results
save('heatmapresults.mat','Pert','Time','Iter','NNZ','NNZp','NNZn','KTAU')
