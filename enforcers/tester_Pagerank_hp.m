%% Enforcing Pagerank

clear; clc; close all;

addpath('/opt/matlabauxiliaries/advanpix-4.8.0/')

mp.Digits(71); % Set quadruple precision

load("../matrices/karate.mat");
A       = mp(Problem.A);
alpha = mp(0.8);
n        = size(A,1);
I         = mp(speye(n,n));
e        = mp(ones(n,1));
v        = (mp(1)/mp(n)).*e;  % Teleportation
deg    = A*e;

visualize_instance = false;

mu      = (I - alpha*(spdiags(mp(1)./deg,0,n,n)*A).')\((mp(1)-alpha).*v) ;
muhat                    = mu;
[~, max_ind]           = maxk(double(mu),10);
[~, min_ind]            = mink(double(mu),10);
muhat(max_ind) = mu(min_ind);
muhat(min_ind)  = mu(max_ind);

%% Repetition at different tolerances
tolvec = mp([1e-9,1e-12,1e-15,1e-18,1e-21,1e-24,1e-27,1e-30,1e-33,1e-36,1e-39,1e-42,1e-45]);
betavec = [1,0.5];
P         = mp(spones(A+mp(spdiags(e,0,n,n))));

for i = 1:2
    beta = mp(betavec(i));
    nnzvec = [];
    normvec = [];
    relnormvec = [];

    nnzvec_std = [];
    normvec_std = [];
    relnormvec_std = [];


    for tol = tolvec
        fprintf("Testing tolerance %e.\n",tol)

        if tol >= mp(1e-13)
            fprintf("\n\nRunning classical:\n\n")
            [Delta_origin,stat_origin,mucheck_origin] = ....
                enforce_pagerank(double(A),double(alpha),double(muhat),...
                double(v),double(P),double(beta),tol);
            nnzvec_std = [nnzvec_std,nnz(Delta_origin(abs(Delta_origin) > 1e-10))];
            normvec_std = [normvec_std,norm(double(muhat)-mucheck_origin)];
            relnormvec_std = [relnormvec_std,norm(double(muhat)-mucheck_origin)/norm(double(muhat))];
        end

        fprintf("\n\nRunning mp\n\n")
        [Delta,stat,mucheck,rhat] = enforce_pagerank_hp(A,alpha,muhat,v,P,beta,tol);
        fprintf('Objective: %e\n',beta*norm(Delta,'fro')^2 + (1-beta)*norm(Delta,1));

        if visualize_instance
            figure(1)
            plot(1:n,mu,'--',1:n,muhat,'o',1:n,mucheck,'x','LineWidth',2);

            figure(2)
            semilogy(1:length(stat.primalres),stat.primalres,...
                1:length(stat.dualres),stat.dualres, ...
                1:length(stat.compl),stat.compl,'LineWidth',2);
            xlabel('IPM Iteration')
            legend('Primal Residual','Dual Residual','Complementarity Residual');

            figure(3)
            %Delta(abs(Delta)<1e-10) = 0;
            %Delta = sparse(Delta);
            subplot(1,3,1)
            spy(double(A))
            subplot(1,3,2)
            spy(double(Delta) > 1e-10,'b+')
            subplot(1,3,3)
            spy(double(Delta) < -1e-10,'r-')
        end

        nnzvec  = [nnzvec,nnz(Delta(abs(Delta) > 1e-10))];
        normvec = [normvec,norm(muhat-mucheck)];
        relnormvec = [relnormvec,norm(muhat-mucheck)/norm(muhat)];
    end
    figure(1)
    subplot(1,2,i)
    semilogy(1:length(normvec),normvec,...
        1:length(relnormvec),relnormvec,...
        1:length(normvec_std),normvec_std,...
        1:length(normvec),tolvec(1:length(normvec)),'--','LineWidth',2)
    ylabel('Error on PageRank')
    if i == 1
        legend({'Augmented precision (abs)','Augmented precision (rel)',...
            'Double precision','Tolerance'},'Location','northeast')
    end
    xticks(1:length(relnormvec))
    xticklabels(split(sprintf("%1.1e\n",tolvec(1:length(relnormvec)))))
    axis tight
    title(sprintf('$\\beta = %1.2f$',betavec(i)),'Interpreter','latex')
end




