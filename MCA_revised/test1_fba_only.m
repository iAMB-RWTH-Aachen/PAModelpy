% Conventional FBA test
%
% --------------------- Pedro Saa UC 2023 ----------------------------------
clc,clearvars,close all

% Let us consider the following LP problem
% max. vi
% s.t.
%      S*v          = 0
%      v           <= ub
%     -v           <= -lb
% Bounds
%     v_i positive (unbounded)
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

% Newtork parameters
[m,n] = size(S);

% Definition of reaction bounds
lb      = zeros(n,1);
ub      = 1e2*ones(n,1);
ub(1:2) = 10;               % Glucose uptake (positive direction for consistency)

% Set up optimization problem
params.OutputFlag = 0;      % Gurobi parameter
model.obj= zeros(n,1);      % Null objective
model.A  = sparse([S;...    % Coefficients matrix
                eye(n);...
                -eye(n)]);
b  = [zeros(m,1);ub;-lb];   % Right-hand side definition
LB = zeros(n,1);            % Bounds definition (positive unbounded)
UB = 1e6*ones(n,1);
model.rhs = b;              % Right-hand side
model.lb  = LB;             % Bounds
model.ub  = UB;

% Constraint and model sense
for ix = 1:size(model.A,1)
    if ix <= m
        model.sense(ix) = '=';
    else
        model.sense(ix) = '<';
    end
end
model.modelsense = 'max';   % Model sense
model.vtype      = 'C';     % Variable type

%% Main loop for Control Coefficients calculation for rhs parameters
Crhs = [];
bdual = b(m+1:end);
for jx = 1:n

    % Objective vector
    model.obj(jx) = 1;

    % Solve model
    sol = gurobi(model,params);

    % Extract solution
    vz = sol.objval;            % optimal objective value
    mu = sol.pi(m+1:end);       % shadow prices of inequality constraints    
    
    % Calculate control coefficients
    Crhs = [Crhs;(mu.*bdual)'/vz];

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
    v  = sol.x;                 % optimal primal vector
    mu_ub = sol.pi(m+1:m+n);        % shadow prices of the mu_UB
    mu_lb = sol.pi(m+n+1:m+2*n);    % shadow prices of the mu_LB

    % Calculate control coefficients
    C  = (mu_ub-mu_lb).*v/vz;
    Cvars = [Cvars;C(:)'];

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
xlabel('Decision variable (flux)')
ylabel('Reaction')
