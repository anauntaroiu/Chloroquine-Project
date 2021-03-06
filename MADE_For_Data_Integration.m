

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

% We want both models to grow
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

% Columns represent different growth thresholds
tiger_fba_results{1,2} = '30'; tiger_fba_results{1,3} = '40';
tiger_fba_results{1,4} = '50'; tiger_fba_results{1,5} = '60';
tiger_fba_results{1,6} = '70'; tiger_fba_results{1,7} = '80';

% Rows represent different condition models
tiger_fba_results{2,1} = '4hr_CQ'; tiger_fba_results{3,1} = '4hr_No_CQ';
tiger_fba_results{4,1} = '24hr_CQ'; tiger_fba_results{5,1} = '24hr_No_CQ';

for i = 1:6 % Loop through each model generated from different growth thres.
    % 4hr Condition
    four_hr_models = four_hr_files{i}; % Call a set of models
    four_hr_CQ = fba(four_hr_models.models{1,1}); % Solve for biomass for CQ model
    four_hr_No_CQ = fba(four_hr_models.models{1,2}); % Solve for biomass for No CQ model
    
    % 24hr Condition
    two_four_hr_models = two_four_hr_files{i}; % Call a set of models
    two_four_hr_CQ = fba(two_four_hr_models.models{1,1}); % Solve for biomass for CQ model
    two_four_hr_No_CQ = fba(two_four_hr_models.models{1,2}); % Solve for biomass for No CQ model
    
    % Fill in cell table with biomass values for each model
    tiger_fba_results{2,(i+1)} = four_hr_CQ.val;
    tiger_fba_results{3,(i+1)} = four_hr_No_CQ.val;
    tiger_fba_results{4,(i+1)} = two_four_hr_CQ.val;
    tiger_fba_results{5,(i+1)} = two_four_hr_No_CQ.val;
end

clearvars four_hr_models two_four_hr_models two_four_hr_CQ four_hr_No_CQ i two_four_hr_CQ two_four_hr_No_CQ four_hr_files two_four_hr_files 

%% Save MADE results for COBRA

% Save models for each condition
four_hr_files = {pf_made_4hr_30;pf_made_4hr_40;pf_made_4hr_50;...
    pf_made_4hr_60;pf_made_4hr_70;pf_made_4hr_80};
two_four_hr_files = {pf_made_24hr_30;pf_made_24hr_40;pf_made_24hr_50;...
    pf_made_24hr_60;pf_made_24hr_70;pf_made_24hr_80};

% Generate object to hold information on gene states
% Columns = Different growth thrsholds (30 - 80%)
% Rows = 4hr CQ or 24hr CQ
gene_states = cell(2,6);

for i = 1:6 % Loop through models generated at different growth thresholds 
    
    four_hr_models = four_hr_files{i}; % Call 4hr models
    two_four_hr_models = two_four_hr_files{i}; % Call 24hr models
    
    % Generate a table of gene states for 4hr CQ and 4hr No CQ models
    four_hr_table = table(four_hr_models.genes, four_hr_models.gene_states(:,1), four_hr_models.gene_states(:,2),'VariableNames',...
        {'genes','CQ_state','No_CQ_state'});
    
    % Save results in gene state table
    gene_states{1,i} = [four_hr_table.Properties.VariableNames; table2cell(four_hr_table)];
    
    % Generate a table of gene states for 24hr CQ and 24hr No CQ models
    two_four_hr_table = table(two_four_hr_models.genes, two_four_hr_models.gene_states(:,1),two_four_hr_models.gene_states(:,2),'VariableNames',...
        {'genes','CQ_state','No_CQ_state'});
    
    % Save results in gene state table
    gene_states{2,i} = [two_four_hr_table.Properties.VariableNames; table2cell(two_four_hr_table)];
    
end
clearvars i

% Generate an objec to hold flux information
% Columns = Different growth thresholds
% Rows = four_hr_CQ, four_hr_No_CQ, two_four_hr_CQ, two_four_hr_No_CQ
cobra_flux = cell(4,6);

for i = 1:6 % Loop through each threshold
    
    % identify genes to delete (if gene state == 0)
    CQ_index = (cell2mat(gene_states{1,i}(2:end,2)) == 0); % Find gene states that equal zero
    CQ_index = vertcat(logical(0), CQ_index); % shift as first element is header
    four_hr_deleted_CQ = gene_states{1,i}(CQ_index,1); % Store gene state to be deleted
    
    CQ_index = (cell2mat(gene_states{1,i}(2:end,3)) == 0); % Find gene states that equal zero
    CQ_index = vertcat(logical(0), CQ_index); % shift as first element is header
    four_hr_deleted_No_CQ = gene_states{1,i}(CQ_index,1); % Store gene state to be deleted
    
    CQ_index = (cell2mat(gene_states{2,i}(2:end,2)) == 0); % Find gene states that equal zero
    CQ_index = vertcat(logical(0), CQ_index); % shift as first element is header
    two_four_hr_deleted_CQ = gene_states{2,i}(CQ_index,1); % Store gene state to be deleted
    
    CQ_index = (cell2mat(gene_states{2,i}(2:end,3)) == 0); % Find gene states that equal zero
    CQ_index = vertcat(logical(0), CQ_index); % shift as first element is header
    two_four_hr_deleted_No_CQ = gene_states{2,i}(CQ_index,1); % Store gene state to be deleted
    
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

% Check models to make sure they are working proberly

for i = 1:6 % Loop through each growth threshold (30 - 80%)
    
    for j = 1:4 % Loop through each model type (4hr CQ, 4 hr No CQ, 24hr CQ, 24 hr No CQ)
        
    flux = cobra_flux{j,i}{2}; % Record model's flux value
    
        if (flux < .1) % If flux value is tiny
            disp(i)
            disp(j)
            warning('Cobra flux == 0')
        end
        
        if (isempty(flux)) % If there is no flux value
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

% Generate object to hold essential reactions
rxn_KO = cell(4,6); % each column = growth threshold
rxn_infeas = cell(4,6);
% rows = 4hr_CQ, 4hr_No_CQ, 24hr_CQ, 24hr_No_CQ

for i = 6 % Select 80% growth threshold for Rxn KO study
    disp(i)
    for j = 1:4 % Loop through each model type (4hr_CQ, 4hr_No_CQ, 24hr_CQ, 24hr_No_CQ)
        disp(j)
        model = cobra_flux{j,i}{1}; % Select model
        
        % Run single reaction deletion on all reactions
        [grRatio,~,~,~,~,~] = singleRxnDeletion(model,'FBA',model.rxns); 
        
        % Save results with growth threshold less than 0.1
        lethalKO = model.rxns(grRatio < 0.1);
        
        % Save results with NA growth threshold
        infeas = model.rxns(isnan(grRatio));
        
        % Save results in table
        rxn_KO{j,i} = lethalKO;
        rxn_infeas{j,i} = infeas;
    end
end

% Select 80% growth threshold model
rxn_KO_80 = rxn_KO(:,6);
rxn_Infeas_80 = rxn_infeas(:,6);

% Find number of infeasible KOs in all model types
[m1, ~] = size(rxn_Infeas_80{1});[m2, ~] = size(rxn_Infeas_80{2});
[m3, ~] = size(rxn_Infeas_80{3});[m4, ~] = size(rxn_Infeas_80{4});

if (m1+m2+m3+m4)>0 % If infeasible KOs were found
    warning('infeasible KOs')
end

clearvars m1 m2 m3 m4 rxn_infeas

% 80 consensus CQ essential
consensus_CQ_80 = intersect(rxn_KO_80{1},rxn_KO_80{3});
writetable(cell2table(consensus_CQ_80),'Consensus_CQ_EssentialRxns.xls')

% 80 consensus No CQ essential
consensus_No_CQ_80 = intersect(rxn_KO_80{2},rxn_KO_80{4});
writetable(cell2table(consensus_No_CQ_80),'Consensus_No_CQ_EssentialRxns.xls')

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

