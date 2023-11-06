% Protein allocation formulation (Toy Model)
%
% --------------------- Pedro Saa UC 2023 ----------------------------------
clc,clearvars,close all

% Let us consider the following LP problem
% max. vi
% s.t.
%      S*v          = 0
%      Kin*v  -  E  = 0
%      v           <= ub
%     -v           <= -lb
%      E           <= E_max
%     -E           <= -E_min
%  sum(Ei) + w'*v  <= phi_E0
% Variables
%          v_i, E_i
%

%% Model set up
% Stoichiometric matrix definition
%     sub int byp atp co2 pre bio 
R1 = [ 1,  0,  0,  0,  0,  0,  0];
R2 = [-1,  1,  0,  0,  1,  0,  0];
R3 = [ 0, -1,  1,  1,  0,  0,  0];
R3r= -R3;
R4 = [ 0, -1,  0,  2,  1,  0,  0];
R5 = [ 0, -1,  0,  0,  0,  1,  0];
R6 = [ 0,  0,  0, -1,  0, -1,  1];
R7 = [ 0,  0,  0,  0,  0,  0, -1];
R8 = [ 0,  0,  0,  0, -1,  0,  0];
R9 = [ 0,  0, -1,  0,  0,  0,  0];
S  = [R1;R2;R3;R3r;R4;R5;R6;R7;R8;R9]';

% Additional enzymatic parameters
phi_P0 = 0.60;                  % total enzyme pool
phi_U0 = 0.10;                  % unused enzyme pool
phi_T0 = 0.01;                  % translation sector
phi_E0 = phi_P0-phi_U0-phi_T0;  % net enzyme pool
kcat   = [1, 0.5, 1, 1, 0.5, 0.45, 1.5];  % catalytic turnovers
enz    = numel(kcat);         

% Newtork parameters
[m,n] = size(S);
Kinv  = [diag(kcat.^-1),zeros(enz,n-enz)];
w     = [-0.01, 0, 0, 0, 0, 0, 0, 0.01, 0, 0];     % flux weights allocation

% Definition of reaction bounds
lb      = zeros(n,1);
ub      = 1e2*ones(n,1);
Emin    = zeros(enz,1);
Emax    = phi_E0*ones(enz,1);

% Set up optimization problem
params.OutputFlag = 0;                           % Gurobi parameter
model.obj = zeros(n+enz,1);                      % Null objective
model.A   = sparse([S,zeros(m,enz);...           % Coefficients matrix
                    Kinv,-eye(enz);...
                    eye(n),zeros(n,enz);...
                    -eye(n),zeros(n,enz);...
                    zeros(enz,n),eye(enz);...
                    zeros(enz,n),-eye(enz);...
                    w,ones(1,enz)]);
b  = [zeros(m+enz,1);ub;-lb;Emax;-Emin;phi_E0];     % Right-hand side definition
LB = zeros(n+enz,1);        % Bounds definition (unconstrained)
UB = 1e3*ones(n+enz,1);
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

%% Main loop for CACs and FACs calculation
GRate = []; 
ES    = [];
FS    = [];
PS    = [];
EFS   = [];
ixBio = 8;
vmax  = 1e-3:1e-3:1e-1;
model.obj(ixBio) = 1;       % Objective vector
for jx = 1:numel(vmax)

    % Set new bound for vs_max
    model.rhs(m+enz+1) = vmax(jx);
    sol = gurobi(model,params);

    % Extract solution
    vz    = sol.objval;            % optimal objective value
    GRate = [GRate;vz];
    mu    = sol.pi(m+enz+1:end);   % shadow prices of inequality constraints    
    bdual = model.rhs(m+enz+1:end);

    % Extract solution
    v       = sol.x(1:n);                           % optimal primal vector (fluxes)
    e       = sol.x(n+1:end);                       % optimal primal vector (enzymes)
    mu_e    = sol.pi(m+1:m+enz);                    % shadow prices of e-constraints
    mu_v_ub = sol.pi(m+enz+1:m+enz+n);              % shadow prices of v_max
    mu_v_lb = sol.pi(m+enz+n+1:m+enz+2*n);          % shadow prices of v_min
    mu_e_ub = sol.pi(m+enz+2*n+1:m+2*enz+2*n);      % shadow prices of e_max
    mu_e_lb = sol.pi(m+2*enz+2*n+1:m+3*enz+2*n);    % shadow prices of e_min
    mu_pi   = sol.pi(end);                          % shadow price of total available enzyme pool
    alpha   = sum(e)/phi_E0;

    % Calculate Capacity Sensitivity Coefficients
    FCS = (mu_v_ub - mu_v_lb).*(v/vz);              % Flux Capacity Sensitivities
    ECS = (mu_e_ub - mu_e_lb).*(e/vz);              % Enzyme Capacity Sensitivities
    PCS = mu_pi*(phi_E0/vz);                        % Proteome Capacity Sensitivities
    
    % Check theorem (Strong Duality)    
    assert(abs(sum(FCS)+sum(ECS)+PCS-1)<1e-8)
    
    % Calculate Flux-Enzyme Sensitivity
    FES = mu_e.*(e/vz);

    % Check realtionships
    [sum(FES) + sum(FCS) + (1-alpha)*PCS,-sum(FES) + sum(ECS) + alpha*PCS]
   
    % Save results
    FS  = [FS;FCS(1)];
    ES  = [ES;ECS(1)];
    EFS = [EFS;FES(1)];
    PS  = [PS;PCS];
end

% Plot results for the substrate consumption and growth
subplot(2,1,1)
plot(vmax,GRate,'-k','LineWidth',1)
ylabel('\mu')

subplot(2,1,2)
plot(vmax,FS,'-r','LineWidth',1)
hold on
plot(vmax,EFS,'-b','LineWidth',1)
plot(vmax,PS,'-g','LineWidth',1)
ylabel('Sensitivity coefficient value (-)')
xlabel('v_{s,max}')
legend({'FCS_1','FES_1','PCS'},'Location','northwest')
% 
% subplot(3,1,3)
% plot(vmax,FACs(:,1),'-r','LineWidth',1)
% hold on
% plot(vmax,FACs(:,n+1),'-b','LineWidth',1)
% ylabel('FACs')
% legend({'v_s','e_s'},'Location','northwest')