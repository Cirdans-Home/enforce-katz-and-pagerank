%% Small example
% This is a small descriptive test case showing the general idea of the
% algorithm and using MATLAB quadprog and fmincon routines. These are 
% feasible only for small problems.

clear; clc; close all;

matrix = input(['Select the matrix:\n' ...
    '1) Karate\n' ...
    '2) Ragusa\n' ...
    '3) Grenoble 1107\n' ...
    '4) Sioux Falls\n' ...
    'matrix= ']);
if matrix ~= 1 && matrix ~= 2 && matrix ~= 3 && matrix ~= 4
    matrix = 4;
end

solverstring = sprintf(['Select solver:\n' ...
    '1) fmincon (interior-point)\n' ...
    '2) quadprog\n'...
    '3) quadprog for L1\n'...
    'solver= ']);

solver = input(solverstring);
switch solver
    case 1
        solvername = 'fmincon (interior-point)';
    case 2
        solvername = 'quadprog';
    case 3
        solvername = 'quadprog';
end
doiprofile = input('Do you want a profile (y,n)? ','s');
switch upper(doiprofile)
    case 'Y'
        doiprofile = true;
    case 'N'
        doiprofile = false;
    otherwise
        doiprofile = false;
end


%% Load a test matrix
switch matrix
    case 1
        load('../matrices/karate.mat');
    case 2
        load('../matrices/Ragusa16.mat');
    case 3
        load('../matrices/gre_1107.mat');
    case 4
        load('../matrices/siouxfalls.mat');
end
A = abs(Problem.A);%./norm(Problem.A(:),"inf");
N = size(A,1);
I = speye(N,N);
e = ones(N,1);
%% Katz
% Computation of the "original" Katz centrality
rhoA = eigs(A,1,"largestabs");
rhoA = abs(rhoA);
alpha = 0.5/rhoA;
beta = 0.2;

mu = (I - alpha*A)\e;

%% Plot the graph?
% Do you want to plot the graph? It is usually a messy plot for general
% graphs. Thus the default is false.
doiplot = input("Plot the Matlab graph object? (y/n): ",'s');
switch upper(doiplot)
    case 'Y'
        doiplot = true;
    case 'N'
        doiplot = false;
    otherwise
        doiplot = false;
end

if doiplot
    figure(1)
    if issymmetric(A)
        G = graph(A);
    else
        G = digraph(A);
    end
    subplot(1,4,1);
    if isfield(Problem,'aux')
        if isfield(Problem.aux,'coord')
            x = Problem.aux.coord(:,1);
            y = Problem.aux.coord(:,2);
            if size(Problem.aux.coord,2) == 3
                z = Problem.aux.coord(:,3);
            end
        end
    end
    if exist("x","var") && exist("y","var") && exist("z","var")
        plot(G,"MarkerSize",mu,"XData",x,"YData",y,"ZData",z, ...
            "LineWidth",G.Edges.Weight)
    elseif  exist("x","var") && exist("y","var")
        plot(G,"MarkerSize",mu,"XData",x,"YData",y,...
            "LineWidth",G.Edges.Weight)
    else
        hp = plot(G,"MarkerSize",mu,"LineWidth",G.Edges.Weight);
        x = hp.XData;
        y = hp.YData;
    end
end

%% Target state:
% we modify here the vector obtained with the standard Katz centrality
% as the one we would like to have.
muhat = mu;
switch matrix
    case 1
        % For Karate
        muhat(33) = 0.5*mu(34);
        muhat(5) = 1.5*mu(5);
    case 2
        % For Ragusa
        muhat(22) = 1;
        muhat(5) = 0.5*mu(5);
    case 3
        muhat(10) = 6;
    case 4
        muhat(10) = mu(11);
        muhat(13) = mu(12);
end

%% Pattern selection

patternselectionstring = sprintf(['Select a pattern:\n' ...
    '1) Original pattern of A (default)\n' ...
    'Pattern = ']);
whatpattern = input(patternselectionstring);
switch  whatpattern
    case 1
        % The pattern is the same of the matrix A.
        P = spones(A);
end
X0 = P;
%% Optimization algorithms
% Here the optimization is executed
if doiprofile
    profile on
end

fprintf('Starting %a Solver:\n',solvername);
switch solver
    case 1
        options = optimoptions('fmincon','Algorithm','interior-point',...
            'SpecifyObjectiveGradient',true, ...
            'SpecifyConstraintGradient',true,...
            'Display','iter-detailed','HessianFcn',@hessianfcn,...
            'OptimalityTolerance',1e-13,'HonorBounds',true,...
            'EnableFeasibilityMode',true);
        [ival,jval,~] = find(P);
        lb = -full(reshape(A(sub2ind([N,N],ival,jval)),length(ival),1));
        X0vec = full(reshape(X0(sub2ind([N,N],ival,jval)),length(ival),1));

        tic;
        [xfinalvec,~,~,output] = fmincon(@costFun_fmincon,X0vec,[],[],[], ...
            [],lb,[],@(x) fidelity(x,ival,jval,A,I,alpha,muhat),options);
        timetoc = toc;
        xfinal = sparse(ival,jval,xfinalvec);
    case 2
        options          = optimoptions('quadprog','Algorithm', ...
            'interior-point-convex',...
            'Display','iter-detailed',...
            'OptimalityTolerance',1e-13);
        proj                = pattern_projector(P);
        % Preparing Quantities for quadprog
        reduced_size = size(proj,1);
        Q                     = 2*speye(reduced_size);
        c                      = (proj*reshape(A,N*N,1));
        L                      = kron(muhat.',speye(N))*proj.';
        b                      = (1/alpha).*(muhat-1)-A*muhat+L*c;
        tic;
        [xfinalvec,~,~,output] = ...
            quadprog(Q,-2.*c,[],[],L,b,zeros(reduced_size,1),[],[],options);
        timetoc = toc;
        [ival,jval,~] = find(P);
        xfinalvec    = xfinalvec -c;
        xfinal = sparse(ival,jval,xfinalvec);
    case 3   
        options = optimoptions('quadprog', ...
            'Algorithm','interior-point-convex',...
            'Display','iter-detailed', 'OptimalityTolerance',1e-13);
        proj                 = pattern_projector(P);
        % Preparing Quantities for quadprog
        reduced_size = size(proj,1);
        % L1 Penalty Parameter
        tau                        = (1-beta)/beta;
        Q                          = 2*speye(reduced_size);
        c                           = (proj*reshape(A,N*N,1));
        L                           = kron(muhat.',speye(N))*proj.';
        b                           = (1/alpha).*(muhat-1)-A*muhat +L*c;
        hatQ                     = blkdiag(Q,sparse(reduced_size,reduced_size),sparse(reduced_size,reduced_size));
        hatc                      = [-2.*c; tau.*ones(reduced_size,1); tau.*ones(reduced_size,1)];
        hatL                      = [L,                                    sparse(N,reduced_size), sparse(N,reduced_size);...
            - speye(reduced_size), speye(reduced_size),      -speye(reduced_size)];
        hatb                      = [(1/alpha).*(muhat-1)-A*muhat+L*c; -c];
        tic;
        [xfinalvec_long,~,~,output] = quadprog(hatQ,hatc,[],[],hatL,hatb,zeros(3*reduced_size,1),[],[],options);
        timetoc = toc;
        [ival,jval,~] = find(P);
        xfinalvec    = xfinalvec_long(1:reduced_size) -c;
        xfinal = sparse(ival,jval,xfinalvec);
end
if doiprofile
    profile off
end

%% Compute the optimizate Katz centrality
mufinal = (I - alpha*(A+xfinal))\e;

%% Visualize
figure(1)
subplot(1,4,3)
spy(A);
title('Adjacency')
figure(1)
subplot(1,4,4)
title('Correction sign')
hold on
spy(xfinal > 1e-10,'b+',5)
spy(xfinal < -1e-10,'r_',5)
hold off
xlabel(sprintf('(Blue) Positive %d - (Red) Negative %d',full(sum(xfinal(:) > 1e-10)), ...
    full(sum(xfinal(:) < -1e-10))))
h4 = figure(4);
hm = heatmap(log10(abs(xfinal)),'GridVisible','off','FontSize',30);
s=struct(hm);
s.XAxis.Visible='off';
s.YAxis.Visible='off';



%% Plot graph
if doiplot
    figure(1)
    G2 = digraph(A+xfinal);
    subplot(1,4,2);
    if exist("x","var") && exist("y","var") && exist("z","var")
        hp = plot(G2,"MarkerSize",mu,"XData",x, ...
            "YData",y,"ZData",z);
    else
        hp = plot(G2,"MarkerSize",mu,"XData",x,"YData",y);
    end

    for l=1:G2.numedges
        i = G2.Edges.EndNodes(l,1);
        j = G2.Edges.EndNodes(l,2);
        if xfinal(i,j) < 0
            EdgeColor = 'red';
            LineStyle = '--';
            if A(i,j) + xfinal(i,j) > 0
                highlight(hp,i,j,'EdgeColor',EdgeColor,'LineStyle',LineStyle, ...
                    'LineWidth',A(i,j)+xfinal(i,j))
            end
        elseif xfinal(i,j) > 0
            EdgeColor = 'blue';
            LineStyle = '-';
            if A(i,j) + xfinal(i,j) > 0
                highlight(hp,i,j,'EdgeColor',EdgeColor,'LineStyle',LineStyle, ...
                    'LineWidth',abs(A(i,j)+xfinal(i,j)))
            end
        end

    end
end
%% Print statistics
fprintf("The norm of the perturbation is: %1.4f\n",norm(xfinal,"fro"));
[~,~,values] = find(xfinal);
minabs = full(min(abs(values),[],"all"));
minval = full(min(xfinal,[],"all"));
maxval = full(max(xfinal,[],"all"));
fprintf("Max and min entries of the perturbation are: %e and %e\n", ...
    minval,maxval);
fprintf("Minimal absolute value of the perturbation: %e\n",minabs);

if solver == 1 || solver == 2 || solver == 3
    [~,violation] = fidelity(xfinalvec,ival,jval,A,I,alpha,muhat);
    maxviolation = max(abs(violation),[],'all');
    meanviolation = mean(abs(violation));
end

fprintf("Max violation is: %e\n",full(maxviolation));
fprintf("Mean violation is %e\n",full(meanviolation));

figure(5)
yyaxis left
plot(1:N,mu,'b--',1:N,muhat,'rx',1:N,mufinal,'ro')
yyaxis right
semilogy(1:N,abs(muhat-mufinal));
legend({'Katz','Desired Katz','Obtained Katz','Error'},'Location','northoutside', ...
    'FontSize',13,'NumColumns',4)
xlim([1 N])

%% If profile is active plot the results:
if doiprofile
    profile viewer
end

%% Export figures with matlab2tikz
try
    figure(2)
    matlab2tikz('ragusa_example_adjacency_qpl1.tex','width','0.3\columnwidth')
    figure(3)
    matlab2tikz('ragusa_example_correction_sign_qpl1.tex','width','0.3\columnwidth')
    figure(4)
    set(gcf,'Color','white')
    export_fig('ragusa_example_correction_magnitude_qpl1.pdf')
    figure(5)
    matlab2tikz('ragusa_example_error_qpl1.tex','width','0.6\columnwidth')
catch
    fprintf("Needs matlab2tikz and export_fig to save figures.\n")
end


%% Functionals for fmincon
% 
function [f,g] = costFun_fmincon(y)
% Implementation of both the cost functional and gradient evaluation. The
% gradient is the unconstrained one.
f = norm(y,2)^2;
if nargout > 1
    g = 2*y;
end

end

function Hout = hessianfcn(x,~)
% Implementation of the Hessian, the constraints are linear and therefore
% they do not appear here.
Hout = 2*speye(size(x,1));
end

%% Reduced fmincon
% Objective function for the fmincon function.
function [c,ceq,gc,gceq] = fidelity(c,ival,jval,A,I,alpha,muhat)

Delta = sparse(ival,jval,c);
c = [];
ceq = 1 - (I - alpha*(A+Delta))*muhat;
gc = [];

gceq = zeros(length(ival),size(A,2));
for i = 1:length(ival)
    gceq(i,ival(i)) = alpha*muhat(jval(i));
end

end

%% Projector on the sparsity pattern for reduced form

function [proj] = pattern_projector(P)
%PATTERN_PROJECTOR given a sparse matrix P this function returns the
%projector onto its sparsity pattern
k     = find(P);
N    = size(P,1);
proj = sparse(1:1:length(k),k,ones(length(k),1), length(k),N*N );
end