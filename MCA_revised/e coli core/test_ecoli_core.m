% Protein allocation formulation (without translation sector)
%
% --------------------- Pedro Saa UC 2023 ----------------------------------
clc,clearvars,close all

% Show results
showResults = 3; % 1,2 or 3

% Let us consider the following LP problem
% max. vi
% s.t.
%      S*v          = 0
%      Kin*v  -  E  = 0
%      v           <= ub
%     -v           <= -lb
%      E           <= E_max
%     -E           <= -E_min
%  wUE*v + sum(Ei) <= phi_P0 - phi_U0 - phi_T0
% Variables
%          v_i, E_i
%

%% Load data
load model_data.mat

%% Initialize enzymatic parameters (from excel)
phi_P0 = .14;                   % total enzyme pool
phi_U0 = 0;                     % unused enzyme pool
phi_T0 = 0.049920;              % translation pool
phi_E0 = phi_P0-phi_U0-phi_T0;  % net enzyme pool
enz    = sum(~isnan(keff));     % mapped enzymes
rxnID  = rxns;

% Newtork parameters
[m,n] = size(S);
Kinv  = diag(keff.^-1);
Kinv(isnan(keff),:) = [];

% Index of important reactions
idBio = 25;
idGlc = 50;

% Parameters of the translation sector
wT = zeros(1,n);
wT(idBio) = 0.002944;    % g*h/gDCW

% Parameters of the unused enzyme sector
wUE = zeros(1,n);
wUE(idGlc) = 0;

w = wUE + wT;

% Definition of enzyme bounds
Emin    = zeros(enz,1);
Emax    = phi_E0*ones(enz,1);

% Set up optimization problem
model.obj= zeros(n+enz,1);    % Null objective
model.A  = sparse([S,zeros(m,enz);...           % Coefficients matrix
    Kinv,-eye(enz);...
    eye(n),zeros(n,enz);...
    -eye(n),zeros(n,enz);...
    zeros(enz,n),eye(enz);...
    zeros(enz,n),-eye(enz);...
    w,ones(1,enz)]);
b  = [zeros(m,1);zeros(enz,1);ub;-lb;Emax;-Emin;phi_E0];     % Right-hand side definition
LB = -1e6*ones(n+enz,1);      % Bounds definition (unconstrained)
UB = 1e6*ones(n+enz,1);

% Assign values to the gurobi variable structure
model.rhs = b;              % Right-hand side
model.lb  = LB;             % Bounds
model.ub  = UB;

% Define constraint and model sense
for ix = 1:size(model.A,1)
    if ix <= m+enz
        model.sense(ix) = '=';
    else
        model.sense(ix) = '<';
    end
end
model.modelsense = 'max';   % Model sense
model.vtype      = 'C';     % Variable type

% Set optimization parameters
params.FeasibilityTol = 1e-9;
params.OptimalityTol  = 1e-9;
params.OutputFlag     = 0;          % Gurobi parameter


%% Main loop for Capacity Allocation Coeficients calculations
Crhs = [];
bdual = b(m+enz+1:end);
for jx = 1:n

    % Objective vector
    model.obj(jx) = 1;

    % Solve model
    sol = gurobi(model,params);

    % Extract solution
    vz = sol.objval;            % optimal objective value
    mu = sol.pi(m+enz+1:end);   % shadow prices of inequality constraints

    % Calculate control coefficients
    Ctemp = (bdual.*mu)'/vz;
    Crhs  = [Crhs;Ctemp];

    % Show results for the biomass optimization
    if jx==idBio
        disp('Max. biomass growth')
        opt_mu = sol.objval
        nonZeros = (abs(Ctemp)>1e-8);
        disp('Non-zero Capacity Control Coefficients')
        Ctemp(nonZeros)
        disp('Capacity Control Coefficients sum')
        sum(Ctemp)
    end

    % Clear objetive vector
    model.obj(jx) = 0;
end

% Check summation result
% assert(all(max(abs(sum(Crhs,2)-1))<1e-6))
if showResults==1
    disp('% Deviation from theoretical result')
    control_sum = sum(Crhs,2);
    dev_from_theory = 100*abs(sum(Crhs,2)-1);
    table(rxns,control_sum,dev_from_theory)
end

% Plot heatmap
figure(1)
subplot(3,1,1)
h = heatmap(Crhs);
colormap jet
h.XDisplayLabels=repmat({''},numel(h.XDisplayLabels),1);
h.YDisplayLabels=repmat({''},numel(h.YDisplayLabels),1);
h.Title = 'Capacity Control Coefficients';
xlabel('Right-hand side (Capacity parameter)')
ylabel('Reaction')

%% Main loop for Flux Control Coefficients
Cvars = [];
for jx = 1:n

    % Objective vector
    model.obj(jx) = 1;

    % Solve model
    sol = gurobi(model,params);

    % Extract solution
    vz = sol.objval;                        % optimal objective value
    v  = sol.x(1:n);                        % optimal primal vector (fluxes)
    e  = sol.x(n+1:end);                    % optimal primal vector (enzymes)
    mu_e  = sol.pi(m+1:m+enz);              % shadow prices of the E-constraints
    mu_ub = sol.pi(m+enz+1:m+enz+n);        % shadow prices of the mu_UB
    mu_lb = sol.pi(m+enz+n+1:m+enz+2*n);    % shadow prices of the mu_LB
    mu_E  = sol.pi(end);                    % shadow price of the total enzyme pool

    % Calculate control coefficients
    Ctemp = [(mu_ub-mu_lb+mu_E*wUE').*v;mu_e.*e]/vz;
    Cvars = [Cvars;Ctemp'];

    % Show results for the biomass optimization
    if jx==idBio
        disp('Max. biomass growth')
        opt_mu = sol.objval
        nonZeros = (abs(Ctemp)>1e-8);
        disp('Non-zero Flux Control Coefficients')
        Ctemp(nonZeros)
        disp('Flux Control Coefficients sum')
        sum(Ctemp)
    end

    % Clear objetive vector
    model.obj(jx) = 0;
end

% Check summation result
if showResults==2
    disp('% Deviation from theoretical results')
    control_sum = sum(Cvars,2);
    dev_from_theory = 100*abs(sum(Cvars,2)-1);
    table(rxnID,control_sum,dev_from_theory);
end

% Plot heatmap
figure(1)
subplot(3,1,2)
h = heatmap(Cvars);
colormap jet
h.XDisplayLabels=repmat({''},numel(h.XDisplayLabels),1);
h.YDisplayLabels=repmat({''},numel(h.YDisplayLabels),1);
h.Title = 'Flux Control Coefficients';
xlabel('Decision variable (flux, enz)')
ylabel('Reaction')

%% Main loop Connectivity Relationships for enzyme-related dual variables
Cvars = [];
for jx = 1:n

    % Objective vector
    model.obj(jx) = 1;

    % Solve model
    sol = gurobi(model,params);

    % Extract solution
    vz = sol.objval;                              % optimal objective value
    v  = sol.x(1:n);                              % optimal primal vector (fluxes)
    e  = sol.x(n+1:end);                          % optimal primal vector (enzymes)
    alpha = sum(e)/phi_E0;                        % utilized enzyme fraction
    mu_e  = sol.pi(m+1:m+enz);                    % shadow prices of the E-constraints
    mu_ub = sol.pi(m+enz+1:m+enz+n);              % shadow prices of the mu_UB
    mu_lb = sol.pi(m+enz+n+1:m+enz+2*n);          % shadow prices of the mu_LB
    e_max = sol.pi(m+2*n+enz+1:m+2*n+2*enz);      % shadow prices of the epsilon_max
    e_min = sol.pi(m+2*n+2*enz+1:m+2*n+3*enz);    % shadow prices of the epsilon_min
    lambda = sol.pi(1:m);
    pi     = sol.pi(end);

    % Calculate control coefficients
    Ctemp = [e.*(-mu_e + e_max - e_min);alpha*pi]/vz;
    Cvars = [Cvars;Ctemp'];

    % Show results for the biomass optimization
    if jx==idBio
        disp('Max. biomass growth')
        opt_mu = sol.objval
        nonZeros = (abs(Ctemp)>1e-8);
        disp('Non-zero Enzyme Allocation Coefficients')
        Ctemp(nonZeros)
        disp('Enzyme Allocation Coefficients sum')
        sum(Ctemp)
    end

    % Clear objetive vector
    model.obj(jx) = 0;
end

% Check summation result
if showResults==3
    disp('Deviation from theoretical results')
    control_sum = sum(Cvars,2);
    dev_from_theory = sum(Cvars,2);
    table(rxnID,control_sum,dev_from_theory)
end

% Plot heatmap
figure(1)
subplot(3,1,3)
h = heatmap(Cvars);
colormap jet
h.XDisplayLabels=repmat({''},numel(h.XDisplayLabels),1);
h.YDisplayLabels=repmat({''},numel(h.YDisplayLabels),1);
h.Title = 'Enzyme Allocation Coefficients';
xlabel('Enzymes + Total enzyme pool')
ylabel('Reaction')