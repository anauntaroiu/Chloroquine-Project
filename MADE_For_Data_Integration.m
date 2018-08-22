

%% Initialize COBRA
initCobraToolbox % this only needs to be run once
changeCobraSolver('gurobi6');

%% Initialize TIGER toolbox
start_tiger
set_solver('gurobi');
set_solver_option('MaxTime',60*60); % max time a TIGER simulation is allowed to run
set_solver_option('IntFeasTol',1e-8); % cutoff number for interpreting a value as 0

%% import expression data 
gene_expression_4hr = readtable('Expression_Data_For_MADE_CQ_4hr.csv','Delimiter',',');
gene_expression_24hr = readtable('Expression_Data_For_MADE_CQ_24hr.csv','Delimiter',',');

%% import model (in both cobra and tiger format for speed)
load pf_model_cobra.mat;
pf_cobra = model;
load pf_model_tiger.mat;

%% Test each model's ability to produce biomass 
% Values should be the same for cobra and tiger models 
optimizeCbModel(pf_cobra) % Solve a flux balance analysis problem
fba(pf_tiger)

%% Specify a minimum growth rate (biomass flux) for each model
pf_obj_flux_result = fba(pf_tiger);
pf_obj_flux = pf_obj_flux_result.val;

% Want both models to grow
obj_value_desired = 1 / 24 / 7;
pf_obj_frac = obj_value_desired / pf_obj_flux;
required_growth = pf_obj_frac;

% set p value cutoff, this means that only genes with a p value less than
% this will be interpreted as changing
% note, a significance threshold of 0.05 is still used
fdr_threshold = 0.5;

% ensure genes are in the same format as in model
gene_in_pf_tiger_4hr = ismember(gene_expression_4hr.ORF_old, pf_tiger.genes);
gene_expression_made_4hr = gene_expression_4hr(gene_in_pf_tiger_4hr,:);

gene_in_pf_tiger_24hr = ismember(gene_expression_24hr.ORF_old, pf_tiger.genes);
gene_expression_made_24hr = gene_expression_24hr(gene_in_pf_tiger_24hr,:);

    
    %% RUN MADE with growth threshold 80%

    tic
    pf_made_4hr_80 = made(...
        pf_tiger,... % model
        gene_expression_made_4hr.logFC,... % fold change (log)
        gene_expression_made_4hr.P_Value,... % p values
        'gene_names',gene_expression_made_4hr.ORF_old,...
        'obj_frac', 0.8, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_24hr_80 = made(...
        pf_tiger,...
        gene_expression_made_24hr.logFC,...
        gene_expression_made_24hr.P_Value,...
        'gene_names',gene_expression_made_24hr.ORF_old,...
        'obj_frac', 0.8, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% RUN MADE with growth threshold 70%

    tic
    pf_made_4hr_70 = made(...
        pf_tiger,... % model
        gene_expression_made_4hr.logFC,... % fold change (log)
        gene_expression_made_4hr.P_Value,... % p values
        'gene_names',gene_expression_made_4hr.ORF_old,...
        'obj_frac', 0.7, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_24hr_70 = made(...
        pf_tiger,...
        gene_expression_made_24hr.logFC,...
        gene_expression_made_24hr.P_Value,...
        'gene_names',gene_expression_made_24hr.ORF_old,...
        'obj_frac', 0.7, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% RUN MADE with growth threshold 60%

    tic
    pf_made_4hr_60 = made(...
        pf_tiger,... % model
        gene_expression_made_4hr.logFC,... % fold change (log)
        gene_expression_made_4hr.P_Value,... % p values
        'gene_names',gene_expression_made_4hr.ORF_old,...
        'obj_frac', 0.6, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_24hr_60 = made(...
        pf_tiger,...
        gene_expression_made_24hr.logFC,...
        gene_expression_made_24hr.P_Value,...
        'gene_names',gene_expression_made_24hr.ORF_old,...
        'obj_frac', 0.6, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
        
%% RUN MADE with growth threshold 50%

    tic
    pf_made_4hr_50 = made(...
        pf_tiger,... % model
        gene_expression_made_4hr.logFC,... % fold change (log)
        gene_expression_made_4hr.P_Value,... % p values
        'gene_names',gene_expression_made_4hr.ORF_old,...
        'obj_frac', 0.5, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_24hr_50 = made(...
        pf_tiger,...
        gene_expression_made_24hr.logFC,...
        gene_expression_made_24hr.P_Value,...
        'gene_names',gene_expression_made_24hr.ORF_old,...
        'obj_frac', 0.5, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    
%% RUN MADE with growth threshold 40%

    tic
    pf_made_4hr_40 = made(...
        pf_tiger,... % model
        gene_expression_made_4hr.logFC,... % fold change (log)
        gene_expression_made_4hr.P_Value,... % p values
        'gene_names',gene_expression_made_4hr.ORF_old,...
        'obj_frac', 0.4, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_24hr_40 = made(...
        pf_tiger,...
        gene_expression_made_24hr.logFC,...
        gene_expression_made_24hr.P_Value,...
        'gene_names',gene_expression_made_24hr.ORF_old,...
        'obj_frac', 0.4, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc    

    
%% RUN MADE with growth threshold 30%

    tic
    pf_made_4hr_30 = made(...
        pf_tiger,... % model
        gene_expression_made_4hr.logFC,... % fold change (log)
        gene_expression_made_4hr.P_Value,... % p values
        'gene_names',gene_expression_made_4hr.ORF_old,...
        'obj_frac', 0.3, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc
    % all other parameters at default
    
    tic
    pf_made_24hr_30 = made(...
        pf_tiger,...
        gene_expression_made_24hr.logFC,...
        gene_expression_made_24hr.P_Value,...
        'gene_names',gene_expression_made_24hr.ORF_old,...
        'obj_frac', 0.3, ... % default = 0.3
        'p_thresh',fdr_threshold,... 
        'set_IntFeasTol',1e-20,...
        'log_fold_change', true);
    toc

  

%% confirm growth %% NOTE Models are sometimes infeasible- if so, next double check they are feasible in next code block

% Save models for each condition
four_hr_files = {pf_made_4hr_30;pf_made_4hr_40;pf_made_4hr_50;...
    pf_made_4hr_60;pf_made_4hr_70;pf_made_4hr_80};
two_four_hr_files = {pf_made_24hr_30;pf_made_24hr_40;pf_made_24hr_50;...
    pf_made_24hr_60;pf_made_24hr_70;pf_made_24hr_80};

% Generate cell to hold results
tiger_fba_results = cell(5,7,1); %preallocate results
tiger_fba_results{1,1} = 'Model';
tiger_fba_results{1,2} = '30'; tiger_fba_results{1,3} = '40';
tiger_fba_results{1,4} = '50'; tiger_fba_results{1,5} = '60';
tiger_fba_results{1,6} = '70'; tiger_fba_results{1,7} = '80';
tiger_fba_results{2,1} = '4hr_CQ'; tiger_fba_results{3,1} = '4hr_No_CQ';
tiger_fba_results{4,1} = '24hr_CQ'; tiger_fba_results{5,1} = '24hr_No_CQ';

for i = 1:6 % Loop through each model
    four_hr_models = four_hr_files{i};
    four_hr_CQ = fba(four_hr_models.models{1,1});
    four_hr_No_CQ = fba(four_hr_models.models{1,2});
    two_four_hr_models = two_four_hr_files{i};
    two_four_hr_CQ = fba(two_four_hr_models.models{1,1});
    two_four_hr_No_CQ = fba(two_four_hr_models.models{1,2});
    tiger_fba_results{2,(i+1)} = four_hr_CQ.val;
    tiger_fba_results{3,(i+1)} = four_hr_No_CQ.val;
    tiger_fba_results{4,(i+1)} = two_four_hr_CQ.val;
    tiger_fba_results{5,(i+1)} = two_four_hr_No_CQ.val;
end

clearvars four_hr_models two_four_hr_models two_four_hr_CQ four_hr_No_CQ i two_four_hr_CQ two_four_hr_No_CQ four_hr_files two_four_hr_files 

%% Save MADE results for COBRA
four_hr_files = {pf_made_4hr_30;pf_made_4hr_40;pf_made_4hr_50;...
    pf_made_4hr_60;pf_made_4hr_70;pf_made_4hr_80};
two_four_hr_files = {pf_made_24hr_30;pf_made_24hr_40;pf_made_24hr_50;...
    pf_made_24hr_60;pf_made_24hr_70;pf_made_24hr_80};

gene_states = cell(2,6);
for i = 1:6
    four_hr_models = four_hr_files{i};
    two_four_hr_models = two_four_hr_files{i};
    % row 1 = 4hr CQ treatment
    four_hr_table = table(four_hr_models.genes, four_hr_models.gene_states(:,1), four_hr_models.gene_states(:,2),'VariableNames',...
        {'genes','CQ_state','No_CQ_state'});
    gene_states{1,i} = [four_hr_table.Properties.VariableNames; table2cell(four_hr_table)];
    % row 2 = 24hr CQ treatment
    two_four_hr_table = table(two_four_hr_models.genes, two_four_hr_models.gene_states(:,1),two_four_hr_models.gene_states(:,2),'VariableNames',...
        {'genes','CQ_state','No_CQ_state'});
    gene_states{2,i} = [two_four_hr_table.Properties.VariableNames; table2cell(two_four_hr_table)];
    % columns = growth thresholds
end
clearvars i

cobra_flux = cell(4,6); % each column = growth threshold
% rows = four_hr_CQ, four_hr_No_CQ, two_four_hr_CQ, two_four_hr_No_CQ
for i = 1:6 % for each threshold
    % identify genes to delete (if gene state == 0)
    CQ_index = (cell2mat(gene_states{1,i}(2:end,2)) == 0);
    CQ_index = vertcat(logical(0), CQ_index); % shift as first element is header
    four_hr_deleted_CQ = gene_states{1,i}(CQ_index,1);
    CQ_index = (cell2mat(gene_states{1,i}(2:end,3)) == 0);
    CQ_index = vertcat(logical(0), CQ_index);
    four_hr_deleted_No_CQ = gene_states{1,i}(CQ_index,1);
    CQ_index = (cell2mat(gene_states{2,i}(2:end,2)) == 0);
    CQ_index = vertcat(logical(0), CQ_index);
    two_four_hr_deleted_CQ = gene_states{2,i}(CQ_index,1);
    CQ_index = (cell2mat(gene_states{2,i}(2:end,3)) == 0);
    CQ_index = vertcat(logical(0), CQ_index);
    two_four_hr_deleted_No_CQ = gene_states{2,i}(CQ_index,1);
    
    % delete genes from model
    four_hr_CQ = deleteModelGenes(pf_cobra,four_hr_deleted_CQ);
    four_hr_No_CQ = deleteModelGenes(pf_cobra,four_hr_deleted_No_CQ);
    two_four_hr_CQ = deleteModelGenes(pf_cobra,two_four_hr_deleted_CQ);
    two_four_hr_No_CQ = deleteModelGenes(pf_cobra,two_four_hr_deleted_No_CQ);
    
    %predict flux once genes are deleted
    four_hr_CQ_flux = optimizeCbModel(four_hr_CQ);
    four_hr_No_CQ_flux = optimizeCbModel(four_hr_No_CQ);
    two_four_hr_CQ_flux = optimizeCbModel(two_four_hr_CQ);
    two_four_hr_No_CQ_flux = optimizeCbModel(two_four_hr_No_CQ);

    % store data
    cobra_flux{1,i} = {four_hr_CQ four_hr_CQ_flux.f};
    cobra_flux{2,i} = {four_hr_No_CQ four_hr_No_CQ_flux.f};
    cobra_flux{3,i} = {two_four_hr_CQ two_four_hr_CQ_flux.f};
    cobra_flux{4,i} = {two_four_hr_No_CQ two_four_hr_No_CQ_flux.f};
end

for i = 1:6
    for j = 1:4
    flux = cobra_flux{j,i}{2};
        if (flux < .1)
            disp(i)
            disp(j)
            warning('Cobra flux == 0')
        end
        if (isempty(flux))
            disp(i)
            disp(j)
            warning('Infeasible solution')
        end
    end
end
clearvars four_hr_deleted_CQ four_hr_deleted_No_CQ two_four_hr_deleted_CQ two_four_hr_deleted_No_CQ four_hr_CQ
clearvars two_four_hr_CQ two_four_hr_No_CQ four_hr_No_CQ four_hr_files four_hr_table i four_hr_CQ_flux four_hr_models four_hr_No_CQ_flux
clearvars four_hr_CQ_flux_result CQ_index two_four_hr_files two_four_hr_models two_four_hr_table
clearvars two_four_hr_No_CQ_flux two_four_hr_CQ_flux
clearvars pf_made_*

%% Rxn KO studies
rxn_KO = cell(4,6); % each column = growth threshold
rxn_infeas = cell(4,6);

% Find Essential Rxns for 80% growth
for i = 6
    disp(i)
    for j = 1:4
        disp(j)
        model = cobra_flux{j,i}{1};
        [grRatio,~,~,~,~,~] = singleRxnDeletion(model,'FBA',model.rxns);
        lethalKO = model.rxns(grRatio < 0.1);
        infeas = model.rxns(isnan(grRatio));
        rxn_KO{j,i} = lethalKO;
        rxn_infeas{j,i} = infeas;
    end
end
rxn_KO_80 = rxn_KO(:,6);
rxn_Infeas_80 = rxn_infeas(:,6);
[m1, ~] = size(rxn_Infeas_80{1});[m2, ~] = size(rxn_Infeas_80{2});
[m3, ~] = size(rxn_Infeas_80{3});[m4, ~] = size(rxn_Infeas_80{4});
if (m1+m2+m3+m4)>0
    warning('infeasible KOs')
end
clearvars m1 m2 m3 m4 rxn_infeas

% 80 consensus CQ essential
consensus_CQ_80 = intersect(rxn_KO_80{1},rxn_KO_80{3});
% 80 consensus No CQ essential
consensus_No_CQ_80 = intersect(rxn_KO_80{2},rxn_KO_80{4});
% 80 unique 4hr
Unique_4hr_80 = setdiff(setdiff(rxn_KO_80{1}, rxn_KO_80{3}), consensus_No_CQ_80);
writetable(cell2table(Unique_4hr_80),'Unique_4hr_EssentialRxns.xls')
% 80 unique 24hr
Unique_24hr_80 = setdiff(setdiff(rxn_KO_80{3}, rxn_KO_80{1}), consensus_No_CQ_80);
writetable(cell2table(Unique_24hr_80),'Unique_24hr_EssentialRxns.xls')
% 80 unique CQ essential
unique_CQ_80 = setdiff(consensus_CQ_80,consensus_No_CQ_80);
% 80 unique No CQ essential
unique_No_CQ_80 = setdiff(consensus_No_CQ_80,consensus_CQ_80);
%essential to all
all_80 = intersect(consensus_No_CQ_80,consensus_CQ_80);
