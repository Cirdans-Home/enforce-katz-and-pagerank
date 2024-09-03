%This script loads various Graphs and solves the Katz Ranking enforcing
%Problem Using Primal Dual Regularized IPM
clear all;
close all;
clc;
addpath('../Centrality_Forcing_IPM/');
%The path on which all the Graphs are:
QP_problems_path = '../testmatrices'; 

%Finds all the Netlib problems and stores their names in a struct
d = dir(fullfile(QP_problems_path,'*.mat')); 
ndigits = 6;


% addpath('/Users/stefanocipolla/Library/CloudStorage/OneDrive-UniversityofSouthampton/Southampton_Work/Force_Ranking/export_fig')
%Open the file to write the results
fileID = fopen('./Results_Figures/Dataset_Info.txt','a+');
fprintf(fileID,'       Problem      &  Type   &   Size   & NNz   & Connected Comp \\\\   \n'); 

% fileID1 = fopen('./Results_Figures/K_Beta1_mu1.txt','a+');
% fprintf(fileID1,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs     & nnz fact & Kendal & RBO \\\\    \n'); 
% 
% fileID2 = fopen('./Results_Figures/K_L1_mu1.txt','a+');
% fprintf(fileID2,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs     & nnz fact & Kendal & RBO \\\\    \n'); 
% 
% fileID3 = fopen('./Results_Figures/K_Beta1_mu2.txt','a+');
% fprintf(fileID3,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs      & nnz fact & Kendal & RBO \\\\    \n'); 
% 
% fileID4 = fopen('./Results_Figures/K_L1_m2.txt','a+');
% fprintf(fileID4,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs       & nnz fact & Kendal & RBO \\\\    \n'); 
fileID1 = fopen('./Results_Figures/K_Beta1_mu1.txt','a+');
fprintf(fileID1,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs     & nnz fact & Kendal \\\\    \n'); 

fileID2 = fopen('./Results_Figures/K_L1_mu1.txt','a+');
fprintf(fileID2,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs     & nnz fact & Kendal  \\\\    \n'); 

fileID3 = fopen('./Results_Figures/K_Beta1_mu2.txt','a+');
fprintf(fileID3,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs      & nnz fact & Kendal  \\\\    \n'); 

fileID4 = fopen('./Results_Figures/K_L1_m2.txt','a+');
fprintf(fileID4,'Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs       & nnz fact & Kendal  \\\\    \n'); 



seed = 10;
rng(seed)
total_iters                       = 0;
total_time                       = 0;
total_IPM_iters              = 0;
scaling                           = 1;
scaling_option               = 1; % Do not Change Scaling Options
scaling_direction           = 'l'; % Do not Change Scaling Options
tol                                  = 1e-9;
pc_mode                       = 2;
print_mode                    = 3;
problems_converged    = 0;
plot_fig                         = 0;
rho                                = 1e-12;
Struct                           = struct();
Struct.Fact                   = 'chol';
delta                             = rho; 
for k =  5%1:length(d)
   model       = struct();
   load(fullfile(QP_problems_path,d(k).name));
   model.name = d(k).name
   n = size(Problem.A,1);
   %% Computation of the "original" Katz centrality
   I        = speye(n,n);
   e       = ones(n,1);
   rhoA  = eigs(Problem.A,1,"largestabs");
   rhoA  = abs(rhoA);
   alpha = 0.5/rhoA;
   mu     = (I - alpha*Problem.A)\e;
   %% Generating Target Katz Centralities:
   %muhat_1                 =  mu(randperm(n).');
   %muhat_2                 = mu;
   %[~, max_ind]            = maxk(mu,10);
   %[~, min_ind]            = mink(mu,10);
   %muhat_2(max_ind)        = mu(min_ind);
   %muhat_2(min_ind)        = mu(max_ind);
   [~,bestmu] = sort(mu,"descend");   
   muhat_1                   = mu;
   nindex                      = round(0.1*n);
   muhat_1(bestmu(1:nindex)) = mean(mu(bestmu(1:nindex)));

   muhat_2                   = mu;
   muhat_2(1:nindex)   = mu(nindex:-1:1);
   %[~, max_ind]              = maxk(mu,100);
   %swap                      = muhat_2(max_ind(1:50));
   %muhat_2(max_ind(1:50))    = muhat_2(max_ind(51:100));
   %muhat_2(max_ind(51:100))  = swap;



   %% Pattern of the Perturbation Equals the Pattern of A
   P                             = spones(Problem.A);
   proj                          = pattern_projector(P);  % Projector Onto the Pattern of A
   reduced_size           = size(proj,1);
   
 
   fprintf(fileID,'  %d &    %s      &  %s     &   %d   & %d   & 1 \\\\   \n', k, model.name, Problem.kind, n, nnz(P) ); 


   %%   ----> muhat_1 <----  Solution -- Without -- Sparsity Constraints --
   
   model.H         = 2*speye(reduced_size);
   model.g          = (proj*reshape(Problem.A,n*n,1));
   model.L          = kron(muhat_1.',speye(n))*proj.';
   model.b          = (1/alpha).*(muhat_1-1)-Problem.A*muhat_1+model.L*model.g;
   
  if (scaling == 1)
        DD = Scale_the_problem(model.L,scaling_option,scaling_direction);
        model.L = spdiags(DD,0,size(model.L,1),size(model.L,1))*model.L;  % Apply the left scaling.
        model.b = model.b.*DD;
    end


   % Running the solver
    free_variables = [];
    time                 = 0; 
    tic;
    [xfinalvec,y,z,Info] = PPM_IPM(-2*model.g,model.L,model.b,model.H,free_variables,tol,200,...
                                         pc_mode,print_mode,Struct,rho,delta); 
    time                 = time                 + toc;
    total_time        = total_time         + time;
    opt                   = Info.opt;
    iter                   = Info.ExIt;
    IPMiter             = Info.IPM_It;
    total_IPM_iters = total_IPM_iters+IPMiter;
    total_iters          = total_iters        + iter; % PPM Iters
     % Recover Matrix and Desired Ranking
     [ival,jval,~] = find(P);
     xfinalvec(abs(xfinalvec)<1e-13)=0;
     xfinalvec    = xfinalvec -model.g;
     xfinal = sparse(ival,jval,xfinalvec,n,n);
     % Compute the optimizate Katz centrality
     mufinal = (I - alpha*(Problem.A+xfinal))\e;
     % Compute Correlations
      [K_1,~]    = corr(round(mufinal,ndigits),round(muhat_1,ndigits),'type','Kendall');
      % [rbo_1,~] = rbosimilarity(mufinal,muhat_1,0.1);
     % Print Details
     if (opt == 1)
       % Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs & nnz fact & Kendal & RBO 
       problems_converged = problems_converged + 1;
       fprintf(fileID1,' %d ]  & %.3e                                       & %d          & %1.2f    & %.3e        & %d        & %d          & %d  &   %.2f \\\\   \n',...
                           k,      condest(model.L*model.L.'),  IPMiter,    time,       norm(xfinal,"fro")/norm(Problem.A,"fro"),...
                           full(sum(xfinal(:) > 1e-10)),...
                           full(sum(xfinal(:) < -1e-10)), Info.nnz, K_1); 
    else
      fprintf(fileID1,' %d ]  & %.3e                                         & %d            & %1.2f   & %.3e         & %d        & %d          & %d  &   %.2f -- Non Opt \\\\  \n',...
                           k,        condest(model.L*model.L.'),  IPMiter,    time,       norm(xfinal,"fro")/norm(Problem.A,"fro"),...
                           full(sum(xfinal(:) > 1e-10)),...
                           full(sum(xfinal(:) < -1e-10)), Info.nnz, K_1); 
    end
   
    %%   ----> muhat_2 <----  Solution -- Without -- Sparsity Constraints --
   
   model.H         = 2*speye(reduced_size);
   model.g          = (proj*reshape(Problem.A,n*n,1));
   model.L          = kron(muhat_2.',speye(n))*proj.';
   model.b          = (1/alpha).*(muhat_2-1)-Problem.A*muhat_2+model.L*model.g;
   

   if (scaling == 1)
        DD = Scale_the_problem(model.L,scaling_option,scaling_direction);
        model.L = spdiags(DD,0,size(model.L,1),size(model.L,1))*model.L;  % Apply the left scaling.
        model.b = model.b.*DD;
    end


   % Running the solver
    free_variables = [];
    time                 = 0; 
    tic;
    [xfinalvec,y,z,Info] = PPM_IPM(-2*model.g,model.L,model.b,model.H,free_variables,tol,200,...
                                         pc_mode,print_mode,Struct,rho,delta); 
    time                  = time                 + toc;
    total_time         = total_time        + time;
    opt                    = Info.opt;
    iter                    = Info.ExIt;
    IPMiter              = Info.IPM_It;
    total_IPM_iters = total_IPM_iters+IPMiter;
    total_iters          = total_iters        + iter; % PPM Iters
     % Recover Matrix and Desired Ranking
     [ival,jval,~] = find(P);
     xfinalvec(abs(xfinalvec)<1e-13)=0;
     xfinalvec    = xfinalvec -model.g;
     xfinal = sparse(ival,jval,xfinalvec,n,n);
     % Compute the optimizate Katz centrality
     mufinal = (I - alpha*(Problem.A+xfinal))\e;
     % Compute Correlations
      [K_1,~]    = corr(round(mufinal,ndigits),round(muhat_2,ndigits),'type','Kendall');
     % [rbo_1,~] = rbosimilarity(mufinal,muhat_2,0.1);
     % Print Details
     if (opt == 1)
       % Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs & nnz fact & Kendal & RBO 
       problems_converged = problems_converged + 1;
       fprintf(fileID3,' %d ]  & %.3e                                           & %d          & %1.2f    & %.3e      & %d        & %d          & %d  &   %.2f  \\\\   \n',...
                           k,        condest(model.L*model.L.'),  IPMiter,    time,       norm(xfinal,"fro")/norm(Problem.A,"fro"),...
                           full(sum(xfinal(:) > 1e-10)),...
                           full(sum(xfinal(:) < -1e-10)), Info.nnz, K_1); 
    else
      fprintf(fileID3,' %d ]  & %.3e                                         & %d       & %1.2f    & %.3e             & %d        & %d          & %d  &   %.2f  -- Non Opt \\\\ \n',...
                           k,        condest(model.L*model.L.'),  IPMiter,    time,        norm(xfinal,"fro")/norm(Problem.A,"fro"),...
                           full(sum(xfinal(:) > 1e-10)),...
                           full(sum(xfinal(:) < -1e-10)), Info.nnz, K_1); 
    end





    %% ----> muhat_1 <----  Solution -- With -- Sparsity Constraints
     tau                        = 100;
     Q                          = 2*speye(reduced_size);
     c                           = (proj*reshape(Problem.A,n*n,1));
     L                           = kron(muhat_1.',speye(n))*proj.';
     b                           = (1/alpha).*(muhat_1-1)-Problem.A*muhat_1 +L*c;
     model.H                = blkdiag(Q,sparse(reduced_size,reduced_size),sparse(reduced_size,reduced_size));
     model.g                 = [-2.*c; tau.*ones(reduced_size,1); tau.*ones(reduced_size,1)];
     model.L                 = [L,                                    sparse(n,reduced_size), sparse(n,reduced_size);...
                                          - speye(reduced_size), speye(reduced_size),      -speye(reduced_size)]; 
     model.b                 = [(1/alpha).*(muhat_1-1)-Problem.A*muhat_1+L*c; -c];
   
    
    if (scaling == 1)
        DD = Scale_the_problem(model.L,scaling_option,scaling_direction);
        model.L = spdiags(DD,0,size(model.L,1),size(model.L,1))*model.L;  % Apply the left scaling.
        model.b = model.b.*DD;
    end

     
    % Running the solver
    free_variables = [];
    time = 0; 
    tic;
    [xfinalvec_L1_long,y_L1,z_L1,Info_L1] = PPM_IPM(model.g,model.L,model.b,model.H,free_variables,tol,200,...
                                         pc_mode,print_mode,Struct,rho,delta); 
    time = time + toc;
    total_time = total_time + time;
    opt     = Info_L1.opt;
    iter    = Info_L1.ExIt;
    IPMiter = Info_L1.IPM_It;
    total_IPM_iters = total_IPM_iters+IPMiter;
    total_iters     = total_iters + iter; % PPM Iters
    % Recover Matrix
    xfinalvec_L1 = xfinalvec_L1_long(1:reduced_size);
    xfinalvec_L1(abs(xfinalvec_L1)<1e-13)=0;
    xfinalvec_L1    = xfinalvec_L1  -c;
    xfinal_L1          = sparse(ival,jval,xfinalvec_L1,n,n);
    % Compute the optimizate Katz centrality
    mufinal_L1      = (I - alpha*(Problem.A+xfinal_L1))\e;
     % Compute Correlations
      [K_1_L1,~]    = corr(round(mufinal_L1,ndigits),round(muhat_1,ndigits),'type','Kendall');
      %[rbo_1_L1,~] = rbosimilarity(mufinal_L1,muhat_1,0.1);
     % Print Details
     if (opt == 1)
       % Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs & nnz fact & Kendal & RBO 
       problems_converged = problems_converged + 1;
       fprintf(fileID2,' %d ]  & %.3e                                        & %d        & %1.2f    & %.3e         & %d        & %d          & %d   &  %.2f \\\\  \n',...
                           k,       condest(model.L*model.L.'),  IPMiter,    time,       norm(xfinal,"fro")/norm(Problem.A,"fro"), ...
                           full(sum(xfinal_L1(:) > 1e-10)),...
                           full(sum(xfinal_L1(:) < -1e-10)),  Info_L1.nnz, K_1_L1); 
    else
      fprintf(fileID2,' %d ]  & %.3e                                         & %d       & %1.2e    & %.3e      & %d        & %d          & %d  &   %.2f -- Non Opt \\\\ \n',...
                           k,      condest(model.L*model.L.'),  IPMiter,    time,        norm(xfinal_L1,"fro")/norm(Problem.A,"fro"),...
                           full(sum(xfinal_L1(:) > 1e-10)),...
                           full(sum(xfinal_L1(:) < -1e-10)), Info_L1.nnz, K_1); 
    end
    
   %% ----> muhat_2 <----  Solution -- With -- Sparsity Constraints
     tau                        = 100;
     Q                          = 2*speye(reduced_size);
     c                           = (proj*reshape(Problem.A,n*n,1));
     L                           = kron(muhat_2.',speye(n))*proj.';
     b                           = (1/alpha).*(muhat_2-1)-Problem.A*muhat_2 +L*c;
     model.H                = blkdiag(Q,sparse(reduced_size,reduced_size),sparse(reduced_size,reduced_size));
     model.g                 = [-2.*c; tau.*ones(reduced_size,1); tau.*ones(reduced_size,1)];
     model.L                 = [L,                                    sparse(n,reduced_size), sparse(n,reduced_size);...
                                          - speye(reduced_size), speye(reduced_size),      -speye(reduced_size)]; 
     model.b                 = [(1/alpha).*(muhat_2-1)-Problem.A*muhat_2+L*c; -c];
   
    if (scaling == 1)
        DD = Scale_the_problem(model.L,scaling_option,scaling_direction);
        model.L = spdiags(DD,0,size(model.L,1),size(model.L,1))*model.L;  % Apply the left scaling.
        model.b = model.b.*DD;
    end
 
     
    % Running the solver
    free_variables = [];
    time = 0; 
    tic;
    [xfinalvec_L1_long,y_L1,z_L1,Info_L1] = PPM_IPM(model.g,model.L,model.b,model.H,free_variables,tol,200,...
                                         pc_mode,print_mode,Struct,rho,delta); 
    time = time + toc;
    total_time = total_time + time;
    opt     = Info_L1.opt;
    iter    = Info_L1.ExIt;
    IPMiter = Info_L1.IPM_It;
    total_IPM_iters = total_IPM_iters+IPMiter;
    total_iters     = total_iters + iter; % PPM Iters
    % Recover Matrix
    xfinalvec_L1 = xfinalvec_L1_long(1:reduced_size);
    xfinalvec_L1(abs(xfinalvec_L1)<1e-13)=0;
    xfinalvec_L1    = xfinalvec_L1  -c;
    xfinal_L1          = sparse(ival,jval,xfinalvec_L1,n,n);
    % Compute the optimizate Katz centrality
    mufinal_L1      = (I - alpha*(Problem.A+xfinal_L1))\e;
     % Compute Correlations
      [K_1_L1,~]    = corr(round(mufinal_L1,ndigits),round(muhat_2,ndigits),'type','Kendall');
      %[rbo_1_L1,~] = rbosimilarity(mufinal_L1,muhat_2,0.1);
     % Print Details
     if (opt == 1)
       % Prob. & Dim. & IPM Iter & Time & Norm Sol & Pos Arcs & Neg Arcs & nnz fact & Kendal & RBO 
       problems_converged = problems_converged + 1;
       fprintf(fileID4,' %d ]  & %.3e                                         & %d       & %1.2f    & %.3e      & %d        & %d          & %d  &   %.2f \\\\   \n',...
                           k,       condest(model.L*model.L.'),  IPMiter,    time,       norm(xfinal_L1,"fro"),...
                           full(sum(xfinal_L1(:) > 1e-10)),...
                           full(sum(xfinal_L1(:) < -1e-10)), Info_L1.nnz, K_1_L1); 
    else
      fprintf(fileID4,' %d ]  & %.3e                                         & %d       & %1.2e   & %.3e      & %d        & %d          & %d  &   %.2f  -- Non Opt \\\\ \n',...
                           k,       condest(model.L*model.L.'),  IPMiter,    time,        norm(xfinal_L1,"fro")/norm(Problem.A,"fro"),...
                           full(sum(xfinal_L1(:) > 1e-10)),...
                           full(sum(xfinal_L1(:) < -1e-10)),        Info_L1.nnz, K_1); 
    end

    % Print Final Statistics

    %fprintf(fileID,'  Total time = %.1e, Total IPM Iter =  & %d  \n',...
                              %total_time, total_IPM_iters  );

   if plot_fig == 1
      figure(1)
      spy(Problem.A,25);
      title('Adjacency','FontSize',30)
      figure(2)
      title('Correction sign','FontSize',30)
      hold on
      spy(xfinal > 1e-10,'b+',25)
      spy(xfinal < -1e-10,'r_',25)
      hold off
      xlabel(sprintf('(Blue) Positive %d - (Red) Negative %d',full(sum(xfinal(:) > 1e-10)), ...
                full(sum(xfinal(:) < -1e-10))),'FontSize',30)
      figure(3)
      title('Correction sign L1','FontSize',30)
      hold on
      spy(xfinal_L1 > 1e-10,'b+',25)
      spy(xfinal_L1 < -1e-10,'r_',25)
      hold off
      xlabel(sprintf('(Blue) Positive %d - (Red) Negative %d',full(sum(xfinal_L1(:) > 1e-10)), ...
                full(sum(xfinal_L1(:) < -1e-10))),'FontSize',30)
      figure(4)
      yyaxis left
      plot(1:n,mu,'bv',1:n,muhat_1,'rx',1:n,mufinal,'ro',1:n,mufinal_L1,'r^','MarkerSize',15)
      yyaxis right
      semilogy(1:n,abs(muhat_1-mufinal));
      legend({'Katz','Desired Katz','Obtained Katz','Obtained Katz L1','Error'},'Location','northoutside', ...
                   'FontSize',13,'NumColumns',3)
      xlim([1 n])
      % h5 = figure(5);
      % hm = heatmap(log10(abs(xfinal)),'GridVisible','off','FontSize',30);
      % s=struct(hm);
      % s.XAxis.Visible='off';
      % s.YAxis.Visible='off';
      % h6 = figure(6);
      % hm1 = heatmap(log10(abs(xfinal_L1)),'GridVisible','off','FontSize',30);
      % s1=struct(hm1);
      % s1.XAxis.Visible='off';
      % s1.YAxis.Visible='off';
     figure(7)
     spy(Normal_equation,25);
     title('Normal Equation Pattern','FontSize',30)
     figure(8)
     spy(Normal_equation_L1,25);
     title('Normal Equation L1 Pattern','FontSize',30)

       try
        figure(1)
        set(gcf,'Color','white')
        export_fig(['./Results_Figures/',model.name,'_adjacency.pdf'])
        figure(2)
        set(gcf,'Color','white')
        export_fig(['./Results_Figures/',model.name,'_correction_sign.pdf'])
        figure(3)
        set(gcf,'Color','white')
        export_fig(['./Results_Figures/',model.name,'_correction_sign_l1.pdf'])
        figure(4)
        set(gcf,'Color','white')
        export_fig(['./Results_Figures/',model.name,'_error_l1.pdf'])
        % figure(5)
        % set(gcf,'Color','white')
        % export_fig(['./Results_Figures/',model.name,'_heatmap.pdf'])
        % figure(6)
        % set(gcf,'Color','white')
        % export_fig(['./Results_Figures/',model.name,'_heatmap_l1.pdf'])
        figure(7)
        set(gcf,'Color','white')
        export_fig(['./Results_Figures/',model.name,'_normal_equation.pdf'])
        figure(8)
        set(gcf,'Color','white')
        export_fig(['./Results_Figures/',model.name,'_normal_equation_l1.pdf'])
        catch
        fprintf("Needs matlab2tikz and export_fig to save figures.\n")
       end
   end
end
fprintf(fileID,' The total PPM iterates were: %d, \n The total IPPM iterates were: %d, \n The total time was: %d and %d problems converged.\n',...
        total_iters,total_IPM_iters,total_time,problems_converged);
fprintf(fileID,' The Stopping tol was: %.2e \n ',tol);
fclose(fileID);
