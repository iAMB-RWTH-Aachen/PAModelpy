% Protein allocation formulation (without translation and unused sectors)
%
% --------------------- Pedro Saa UC 2023 ----------------------------------
clc,clearvars,close all

% Let us consider the following LP problem
% max. vi
% s.t.
%      S*v          = 0
%      Kinv*v  - E  = 0
%      v           <= ub
%     -v           <= -lb
%      E           <= E_max
%     -E           <= E_min
%      sum(Ei)     <= phi_P0
% Variables
%          v_i, E_i (positive)
%

%% Model set up
% Stoichiometric matrix definition
%      G  C  A  P  B
R1 = [ 1, 0, 0, 0, 0];
R2 = [-1, 1, 1, 0, 0];
R3 = [-1, 0, 0, 1, 0];
R4 = [ 0, 0,-1, 0, 1];
R5 = [ 0, 0, 0,-1, 1];
R6 = [ 0, 0, 0, 0,-1];
R7 = [ 0,-1, 0, 0, 0];
S  = [R1;R2;R3;R4;R5;R6;R7]';

% Additional enzymatic parameters
phi_P0 = 0.45;                  % total enzyme pool
kcat = [1,0.5,0.25,1.5,0.5];    % catalytic turnovers
enz  = numel(kcat);

% Newtork parameters
[m,n] = size(S);
Kinv  = [diag(1./kcat),zeros(enz,n-enz)];

% Definition of reaction bounds
lb      = zeros(n,1);
ub      = 1e2*ones(n,1);
ub(1:2) = 10;               % Glucose uptake (positive direction for consistency)
Emin    = zeros(enz,1);
Emax    = phi_P0*ones(enz,1);

% Set up optimization problem
params.OutputFlag = 0;      % Gurobi parameter
model.obj= zeros(n+enz,1);  % Null objective
model.A  = sparse([S,zeros(m,enz);...           % Coefficients matrix
                    Kinv,-eye(enz);...
                    eye(n),zeros(n,enz);...
                    -eye(n),zeros(n,enz);...
                    zeros(enz,n),eye(enz);...
                    zeros(enz,n),-eye(enz);...
                    zeros(1,n),ones(1,enz)]);
b  = [zeros(m,1);zeros(enz,1);...
    ub;-lb;Emax;Emin;phi_P0];     % Right-hand side definition
LB = zeros(n+enz,1);              % Bounds definition (unconstrained)
UB = 1e6*ones(n+enz,1);
model.rhs = b;              % Right-hand side
model.lb  = LB;             % Bounds
model.ub  = UB;

% Constraint and model sense
for ix = 1:size(model.A,1)
    if ix <= m+enz
        model.sense(ix) = '=';
    else
        model.sense(ix) = '<';
    end
end
model.modelsense = 'max';   % Model sense
model.vtype      = 'C';     % Variable type

%% Main loop for Control Coefficients calculation for rhs parameters
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
    %     mu = reshape(mu,n,2);
    %     mu = [mu(:,1);mu(:,2)];    % Extract dual variable pairs

    % Calculate control coefficients
    Crhs = [Crhs;(bdual.*mu)'/vz];

    % Clear objetive vector
    model.obj(jx) = 0;
end

% Check summation result
assert(all(max(abs(sum(Crhs,2)-1))<1e-6))

% Plot heatmap
figure(1)
subplot(2,1,1)
heatmap(Crhs)
colormap jet
xlabel('Right-hand side parameter')
ylabel('Reaction')

%% Main loop for Control Coefficients calculation for decision variables
Cvars = [];
for jx = 1:n

    % Objective vector
    model.obj(jx) = 1;

    % Solve model
    sol = gurobi(model,params);
    
    % Extract solution
    vz = sol.objval;            % optimal objective value
    v  = sol.x(1:n);            % optimal primal vector (fluxes)
    e  = sol.x(n+1:end);        % optimal primal vector (enzymes)
    mu_e  = sol.pi(m+1:m+enz);              % shadow prices of the E-constraints
    mu_ub = sol.pi(m+enz+1:m+enz+n);        % shadow prices of the mu_UB
    mu_lb = sol.pi(m+enz+n+1:m+enz+2*n);    % shadow prices of the mu_LB

    % Calculate control coefficients
    C  = [(mu_ub-mu_lb).*v;mu_e.*e]/vz;
    Cvars = [Cvars;C'];

    % Clear objetive vector
    model.obj(jx) = 0;
end

% Check summation result
assert(all(max(abs(sum(Cvars,2)-1))<1e-6))

% Plot heatmap
figure(1)
subplot(2,1,2)
heatmap(Cvars)
colormap jet
xlabel('Decision variable (flux and enz)')
ylabel('Reaction')
