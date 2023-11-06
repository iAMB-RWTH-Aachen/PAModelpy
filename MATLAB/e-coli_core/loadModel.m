clc, clearvars

load e_coli_core.mat

% Reduce model manually
nRxns   = numel(model.rxns);
lb      = zeros(nRxns,1);
ub      = lb;
model.c = ub;
vTol    = 1e-8;                 %  flux tolerance
for ix = 1:nRxns
    model.c(ix) = 1;
    solMin = optimizeCbModel(model,'min');
    solMax = optimizeCbModel(model,'max');
    if abs(solMin.f)>vTol
        lb(ix) = solMin.f;
    end
    if abs(solMax.f)>vTol
        ub(ix) = solMax.f;
    end
    model.c(ix) = 0;
end

% Find zero reactions
ixZeroRxns = (abs(lb)==0)&(abs(ub)==0);
rxns = model.rxns;
S    = model.S;
rxns(ixZeroRxns) = [];
lb(ixZeroRxns)   = [];
ub(ixZeroRxns)   = [];
S(:,ixZeroRxns)  = [];
orphanMets       = all(S==0,2);
S(orphanMets,:)  = [];

% Figure out irreversible rxns
rev = sign(lb).*sign(ub);
rev(rev~=-1) = 0;
rev(rev==-1) = 1;
revRxns      = rxns(rev==1);

% Load enzymatic data
load enzymatic_data.mat
MW   = MW/1e3;                % g E/mmol E
kcat = kcat*3600;             % mmol S/mmol E/h
keff = kcat./MW;              % mmol S/g E/h

% Pre-process enzyme Ids
hits = 0;
keff_temp  = nan(numel(rxns),1);
for jx = 1:numel(rxnIds)

    testID = regexp(rxnIds{jx},'_f','split');
    testID = regexp(testID{1},'_b','split');

    % Map data
    hit = strcmp(rxns,testID{1});
    if any(hit)       
%         [testID{1},rxns{hit}]
        keff_temp(hit) = keff(jx);
        hits = hits+1;
    end
end
hits

% Make model irreversible manually
for ix = 1:numel(rev)
    if rev(ix)==1
        keff_temp(end+1) = keff_temp(ix);
        S(:,end+1) = -S(:,ix);
        ub(end+1)  = -lb(ix);       % make positve upper bound
        lb(end+1)  = 0;             % make new rxn irreversible
        lb(ix)     = 0;             % make original rxn irrversible
        rxns{end+1} = [rxns{ix},'_neg'];
    else
        if lb(ix)<0
            lb_temp = lb(ix);
            lb(ix)  = -ub(ix);
            ub(ix)  = -lb_temp;
            S(:,ix) = -S(:,ix);
        end
    end
end
keff = keff_temp;
table(rxns,keff_temp,lb,ub)

% % Test problem
% problem.A=sparse(S);
% problem.lb=lb;
% problem.ub=ub;
% problem.obj=zeros(numel(lb),1);
% problem.sense='=';
% problem.modelsense='max';
% problem.obj(25)=1;
% sol=gurobi(problem)

% Extract fields for later use
clearvars -except S lb ub keff rxns
save model_data