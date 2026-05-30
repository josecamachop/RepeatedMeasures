function varargout = RM_LiMM_PCA(data, options)

% Function for doing repeated measures-LiMM_PCA.
% Requires the Statistics and Machine Learning toolbox.
% version X
% [M_A, M_B, M_C, ZU, E, Ymodel, fig1, fig2, fig3] = RM_LiMM_PCA(data, options)
% The output elements M_A, M_B, and M_C describe the time effect, treatment + time*treatment effect, 
% and time + treatment + time*treatment effect, respectively. Each is a structure with the following
% fields:
% 
% M: The effect matrix for the factor(s)
% loadings: The back-transformed loadings from PCA on the the effect matrix.
% scores: The scores from PCA on the effect matrices.
% eigen: Eigenvalues for the principal components for the effect matrices.
% pval: P-value for the effects calculated from the GLLR-test
% GLLR_b: Vector containing the bootstrapped GLLR-values
% GLLR_obs: The observed GLLR-value. The p-value is calculated as sum(GLLR_b >
% GLLR_obs)/options.iterations
% scores_boot: A cell array containing the bootstrapped scores.
% loadings_boot: A cell array containing the bootstrapped loadings.
% 
% The outputs ZU and E are the random intercept and residual matrix, respectively. 
% Ymodel is a structure containing the PCA-model describing Y. 
% [fig1, fig2, fig3] are figure handles for the figures for M_A, M_B, and M_C, respectively. 
%
% The function requires a data table (data), which in addition to the response variables also includes the following variables: 
% data.ID: Subject ID vector
% data.timepoint: Vector indicating timepoint
% data.treatment: Vector indicating treatment assignment
% 
% options is a structure with the following required fields:
% options.baseline: Set to either 'ucLDA' or 'cLDA'.
% options.Y_vars: String array with response variable names
% options.ncomp: Number of components from pca(Y) to be included 
% 
% Optional fields:
% options.GLLR: Set to 'yes' to get p-values with the GLLR-test
% options.coding: Determines coding for the treatment variable. Default is sum coded treatment variable
% Set to 'PRC' if treatment should be reference coded. If the treatment variable is a string, the group appearing first in the
% alphabet is the reference group. If it is a number, then lowest number is
% the reference group. If it is logical, then 'false' is the reference category.
% options.iterations: Number of iterations to use when bootstrapping (1000 by default).
% options.CI: Confidence intervals are always calculated by default. To
% avoid this, set options.CI = [];
% options.plot: By default, plots are always produced. To avoid this, set
% options.plot = [];
% options.show_varnames: Set to 'yes' if variable names should be displayed
% on x-axis for loadings.
% options.color: P x 3 matrix where P is the number of response variables,
% and every row is a color vector (e.g. [0 1 0]) indicating the color of
% the variable in the loading plot

%% Generate design- and response matrix
% Make handle for a function to create the design matrix
designmat = @(T, G) [ones(size(T,1),1), T, G, T(:,repmat([1:size(T,2)], 1, size(G,2))).*G(:,kron(1:size(G,2),ones(1,size(T,2))))];

% Set baseline as the reference timepoint ==> Pepe, 20/08/2025, time was
% PRC and I transform it to sum coding
%dummy_time = dummyvar(data.timepoint);
%dummy_time(:,1) =  [];
timepoints = sort(unique(data.timepoint));
for i = 1:length(timepoints)-1
    dummy_time(:,i) = (data.timepoint == timepoints(i)) - (data.timepoint == timepoints(length(timepoints)));
end

% Set sum coding as default for treatment factor
if isfield(options, 'coding') == 0
    options.coding = 'ASCA'; % 'ASCA' is set to default
end

% Create fixed effect design matrix
if isequal(options.coding,'PRC')
    dummy_treatment = dummyvar(nominal(data.treatment));
    dummy_treatment(:,1) = [];
elseif isequal(options.coding,'ASCA')
    treatments = sort(unique(data.treatment));
    for i = 1:length(treatments)-1
        dummy_treatment(:,i) = (data.treatment == treatments(i)) - (data.treatment == treatments(length(treatments)));
    end
end

% Create fixed effect design matrix and specify column number for the effects
% Pepe 20/08/2025: Modify the effects of interest
if isequal(options.baseline,'cLDA')
    X = designmat(dummy_time, dummy_treatment);
    X(:,(1+length(unique(data.timepoint))):(length(unique(data.timepoint))+size(dummy_treatment,2))) = [];
    options.GLLR_effects = {[2:length(unique(data.timepoint))], [1+length(unique(data.timepoint)):length(unique(data.treatment))+length(unique(data.timepoint))-1], [length(unique(data.timepoint))+length(unique(data.treatment)):size(X,2)]};
elseif isequal(options.baseline,'ucLDA')
    X = designmat(dummy_time, dummy_treatment);
    options.GLLR_effects = {[2:length(unique(data.timepoint))], [1+length(unique(data.timepoint)):length(unique(data.treatment))+length(unique(data.timepoint))-1], [length(unique(data.timepoint))+length(unique(data.treatment)):size(X,2)]};
end

Y_vars = options.Y_vars;
Y = data{:,Y_vars};
data_X = data(:, setdiff(data.Properties.VariableNames, options.Y_vars));
G = data.ID;
Z = ones(size(data,1),1);
effects = options.GLLR_effects;

% List all timepoints and treatments
timepoints = sort(unique(data.timepoint), 'ascend');
treatments = sort(unique(data.treatment), 'ascend');

% Pepe 20/08/2025: Scale outside this function
% % Scaling response variables to their baseline standard deviation
% scalingfactor = std(Y(data.timepoint==timepoints(1),:)); % Scale to baseline standard deviation
% for d = 1:size(Y,2)
%     Y(:,d) = Y(:,d)./scalingfactor(d);
% end

%% Fit model
% Pretransform data by PCA
[loadings_Y, scores_Y, eigen_Y] = pca(Y, 'NumComponents', options.ncomp);

% Run LMMs and collect coefficients, residuals and log likelihoods
for i = 1:size(scores_Y,2)
    lme = fitlmematrix(X, scores_Y(:,i), Z, G, 'FitMethod', 'REML');
    b(:,i) = lme.Coefficients.Estimate;
    E(:,i) = residuals(lme);
    LL_obs_h1(i) = lme.ModelCriterion.LogLikelihood;
    U(:,i) = randomEffects(lme);
    [vca,mse(i)] = covarianceParameters(lme); % Pepe 25/08/2025: compute the MS from the variance components
    msca(i) = mse(i) + length(timepoints)*vca{1};
end

% Make effect matrices for each of the effect combinations specified in 'effects'-input
for k = 1:length(effects)
    M_f = X(:,effects{k})*b(effects{k},:);
    [loadings_f, scores_f, eigen_f] = pca(M_f, 'NumComponents', rank(M_f), 'Centered', options.center);
    varargout{k} = struct('M', M_f, 'loadings', loadings_Y*loadings_f, 'scores', scores_f, 'eigen', eigen_f);
end

% Get residual matrix E and random effect matrix ZU and add to output
varargout{length(varargout)+1} = struct('M',full(designMatrix(lme, 'random'))*U,'MSS',msca);  % Pepe 25/08/2025: compute the MS from the variance components
varargout{length(varargout)+1} = struct('M',E,'MSS',mse);  % Pepe 25/08/2025: compute the MS from the variance components

% Get scores, loadings, and eigenvalues for PCA on original Y-data and add to output
varargout{length(varargout)+1} = struct('scores', scores_Y, 'loadings', loadings_Y, 'eigen', eigen_Y);

%% GLLR-test
if isfield(options, 'iterations') == 0
    options.iterations = 1000; % set default nummber to 1000
end

if isfield(options, 'GLLR') == 0
    options.GLLR = ''; % set no GLLR-test as default
end

if isequal(options.GLLR,'yes')
    for e = 2
        
        for i = 1:size(scores_Y,2)
            lme = fitlmematrix(X(:,setdiff(1:size(X,2), options.GLLR_effects{e})) , scores_Y(:,i), Z, G, 'FitMethod', 'REML'); % Mixed model without the tested effects
            b_h0(:,i) = lme.Coefficients.Estimate;
            rEffects_h0 = randomEffects(lme);
            [randomeffectvariance_h0(e,i), residualvariance(e,i)] = covarianceParameters(lme);
            LL_obs_h0(i) = lme.ModelCriterion.LogLikelihood;
        end
   
        b_h0_2{e} = b_h0; % Get the coefficients for the null model for each set of tested effects
        gllr_obs(e) = 2*(sum(LL_obs_h1 - LL_obs_h0)); % Get the GLLRs for each set of tested effects
        clear b_h0
    end
    
    dummy_z = full(designMatrix(lme,'random'));
    
    % Estimate distribution of the GLLRs under null hypotheses by bootstrapping
    for e = 2
        for i = 1:options.iterations
            
            % Generate bootstrapped random effects and residuals
            for f = 1:size(scores_Y,2)
                U_b(:,f) = normrnd(0, sqrt(randomeffectvariance_h0{e,f}), size(rEffects_h0, 1), 1);
                E_b(:,f) = normrnd(0, sqrt(residualvariance(e,f)), size(scores_Y, 1), 1);
            end
            
            % Generate bootstrapped Y-matrix
            Y_b = X(:, setdiff(1:size(X,2), options.GLLR_effects{e}))*b_h0_2{1,e} + dummy_z*U_b + E_b; 
            
            % Estimate null model for bootstrapped Y
            for t = 1:size(scores_Y,2)
                lme = fitlmematrix(X(:, setdiff(1:size(X,2), options.GLLR_effects{e})), Y_b(:,t), Z, G);
                LL_h0_b(t) = lme.ModelCriterion.LogLikelihood;
            end
            
            % Estimate full model for bootstrapped Y
            for t = 1:size(scores_Y,2)
                lme = fitlmematrix(X, Y_b(:,t), Z, G);
                LL_h1_b(t) = lme.ModelCriterion.LogLikelihood;
            end
            
            teststat(i) = 2*(sum(LL_h1_b - LL_h0_b)); % Compute test statistic
        end
        bootstrapped_gll{e} = teststat;
    end
    
    for r = 2
        varargout{r}.pval = (sum(bootstrapped_gll{r} >=gllr_obs(r)) + 1)/(options.iterations + 1);
        varargout{r}.GLLR_b = bootstrapped_gll{r};
        varargout{r}.GLLR_obs = gllr_obs(r);
    end
end


% %% Permutation test: Pepe's choice 20/08/2025
% 
if options.permute == 'yes'

    for ff = 1:4  % Pepe 20/08/2025: Add the individual effect  % Pepe 25/08/2025: compute the MS from the variance components
        if ff~=4
            teststat_perm_obs(ff) = sum(sum((varargout{ff}.M).^2));
            teststat_perm_obs2(ff) = (sum(sum((varargout{ff}.M).^2)))/sum(varargout{5}.MSS);
        else
            teststat_perm_obs(4) = sum(varargout{4}.MSS);
            teststat_perm_obs2(4) = sum(varargout{4}.MSS)/sum(varargout{5}.MSS); 
        end
        if ff~=2
            teststat_perm_obs3(ff) = teststat_perm_obs2(ff);
        else
            teststat_perm_obs3(2) = sum(sum((varargout{2}.M).^2))/sum(varargout{4}.MSS); 
        end 
    end

    for i = 1:options.iterations

        perms = randperm(size(scores_Y,1)); % Pepe 20/08/2025: Unconstrained permutation on raw data, equivalent to the GLM-based ASCA
        Y_perm = scores_Y(perms,:);

        for u = 1:size(scores_Y,2)
            lme_perm = fitlmematrix(X, Y_perm(:,u), Z, G, 'FitMethod', 'REML');
            b_perm(:,u) = lme_perm.Coefficients.Estimate;
            E_perm(:,u) = residuals(lme_perm);
            U_perm(:,u) = randomEffects(lme_perm);
            [vca,mse(u)] = covarianceParameters(lme_perm); % Pepe 25/08/2025: compute the MS from the variance components
            msca(u) = mse(u) + length(timepoints)*vca{1};
        end

        for ff = 1:3
            M_treatment_perm = X(:, options.GLLR_effects{ff})*b_perm(options.GLLR_effects{ff}, :);
            % Pepe 25/08/2025: compute the MS from the variance components
            teststat_perm(ff,i) = sum(sum((M_treatment_perm).^2));
            teststat_perm2(ff,i) = sum(sum((M_treatment_perm).^2))/sum(mse); 
            if ff~=2
                teststat_perm3(ff,i) = teststat_perm2(ff,i);
            else
                teststat_perm3(2,i) = sum(sum((M_treatment_perm).^2))/sum(msca); 
            end
        end
        teststat_perm(4,i) = sum(msca);
        teststat_perm2(4,i) = sum(msca)/sum(mse);
        teststat_perm3(4,i) = teststat_perm2(4,i);
    end

    for ff = 1:4 
        varargout{ff}.pval_perm1p = (sum(teststat_perm(ff,:) > teststat_perm_obs(ff))+1)/(options.iterations+1);
        varargout{ff}.pval_perm2p = (sum(teststat_perm2(ff,:) > teststat_perm_obs2(ff))+1)/(options.iterations+1);
        varargout{ff}.pval_perm3p = (sum(teststat_perm3(ff,:) > teststat_perm_obs3(ff))+1)/(options.iterations+1);    
    end

end


%% Nonparametric bootstrapping to calculate confidence intervals
if isfield(options, 'CI') == 0
    options.CI = 'yes'; % do bootstrapping-CI by default
end

if isequal(options.CI,'yes')
    allpatients = unique(G, 'stable'); % Create vector containing all unique patient IDs

    for i = 1:size(allpatients,1)
        groups_allpatients(i) = unique(data.treatment(find(G==allpatients(i)))); % Find group membership for each patient
    end

    for i = 1:options.iterations;
        
        % Make bootstrap sample
        X_b = [];
        Y_b = [];
        G_b = 0;  

        data_X_b = table();

        if isequal(options.newID,'false')
            for t = 1:length(treatments)
                for k = 1:length(allpatients(groups_allpatients==treatments(t)))
                    idx = find(groups_allpatients == treatments(t));
                    idx = idx(randi(length(idx)));

                    X_add = X(data.ID == allpatients(idx,:),:);
                    X_b = [X_b; X_add];

                    data_X_add = data_X(data.ID == allpatients(idx,:),:);
                    data_X_b = [data_X_b; data_X_add];

                    Y_add = Y(data.ID == allpatients(idx,:),:);
                    Y_b = [Y_b; Y_add];
                end
            end
        elseif isequal(options.newID,'true')
            for t = 1:length(treatments)
                for k = 1:length(allpatients(groups_allpatients==treatments(t)))
                    idx = find(groups_allpatients == treatments(t));
                    idx = idx(randi(length(idx)));

                    X_add = X(data.ID == allpatients(idx,:),:);
                    X_b = [X_b; X_add];

                    add = data(data.ID == allpatients(idx,:),:);
                    if size(data_X_b,1) < 1
                        data_X_b = add;
                        data_X_b.ID = ones(size(add,1),1);
                    else
                        add.ID = ones(size(add,1),1) * (max(unique(data_X_b.ID))+1);
                        data_X_b = [data_X_b; add];
                    end
                end
            end
        end
        G_b = data_X_b.ID;

        clear Y_add X_add G_add
        Z_b = ones(size(X_b,1),1);
        
        % Y_b = data_X_b{:, Y_vars};
        % Scale to baseline SD
        % Pepe 20/08/2025: Scale outside this function
        % scalingfactor = std(Y_b(data_X_b.timepoint==timepoints(1),:)); % Scale to baseline standard deviation
        % for d = 1:size(Y_b,2)
        %    Y_b(:,d) = Y_b(:,d)./scalingfactor(d);
        % end
        
        % PCA-reduction on bootstrapped dataset
        [loadings_Y_b, scores_Y_b, eigen_Y_b] = pca(Y_b, 'NumComponents',options.ncomp);
        
       % Run LMMs and collect coefficients and residuals
        for d = 1:size(scores_Y_b,2)
            lme = fitlmematrix(X_b, scores_Y_b(:,d), Z_b, G_b);
            b(:,d) = lme.Coefficients.Estimate;
        end
           
        % Create variable for storing each of the sub-models
        model_b = cell(1, length(effects));
        model_b2 = cell(1, length(effects));

        % Make effect matrices for each of the effect combinations specified in 'effects'-input
        for q = 1:length(effects)
            model_b{q}.M = X_b(:,effects{q})*b(effects{q},:);
            [model_b{q}.loadings, model_b{q}.scores, model_b{q}.eigen] = pca(model_b{q}.M, 'NumComponents', rank(model_b{q}.M), 'Centered', options.center);
            model_b{q}.loadings = loadings_Y_b*model_b{q}.loadings;
               
            % Get the original scores from the non-bootstrapped data
            clear compare_original compare_rot compare_rot2
            for g = 1:length(treatments)
                for t = 1:length(timepoints)
                    compare_original{g}{t,:} = varargout{q}.scores(find(data.timepoint == timepoints(t) & data.treatment == treatments(g), 1, 'first'),:);
                end
            end
            
            % Store the original scores in a matrix
            compare_original = [cell2mat(compare_original{1}); cell2mat(compare_original{2})];
            
            % Orient all bootstrap-scores so that they match optimally with the sample
            % First we make a matrix (C) containing all combinations of signs
            % to be tested
            n = rank(model_b{q}.scores);
            C = cell(1,n);
            [C{:}] = ndgrid([-1,1]);
            C = cellfun(@(a)a(:),C,'Uni',0); 
            signMatrix = [C{:}]; clear C
            
            V = zeros(1, size(signMatrix,1)); % vector for storing comparison results

            % This loop applies each reflection to the bootstrapped data,
            % compares it to the reference, and the stores the difference
            % in V

            for s = 1:size(signMatrix,1)
                
                R = diag(signMatrix(s,:));
                model_b2{q} = struct('M', model_b{q}.M, 'loadings', model_b{q}.loadings*R, 'scores', model_b{q}.scores*R, 'eigen', model_b{q}.eigen);
                
                [model_b2{q}.loadings, rotationmatrix] = rotatefactors(model_b2{q}.loadings, 'Method', 'procrustes', 'type', 'orthogonal', 'Target', varargout{q}.loadings); 
                model_b2{q}.scores = model_b2{q}.scores*inv(rotationmatrix);
            
                for g = 1:length(treatments)
                    for t = 1:length(timepoints)
                        compare_rot{g}{t,:} = model_b2{q}.scores(find(data_X_b.timepoint == timepoints(t) & data_X_b.treatment == treatments(g), 1, 'first'),:);
                    end
                end
                
                compare_rot2 = [cell2mat(compare_rot{1}); cell2mat(compare_rot{2})];

                V(s) = (norm(compare_original - compare_rot2, 'fro'))^2;
            end
               
            % Choose the reflection yielding the lowest value for V and
            % apply it to the model
            R = diag(signMatrix(find((V == min(V)) == 1, 1, 'first'),:));
            model_b{q} = struct('M', model_b{q}.M, 'loadings', model_b{q}.loadings*R, 'scores', model_b{q}.scores*R, 'eigen', model_b{q}.eigen);
            [model_b{q}.loadings, rotationmatrix] = rotatefactors(model_b{q}.loadings, 'Method', 'procrustes', 'type', 'orthogonal', 'Target', varargout{q}.loadings); 
            model_b{q}.scores = model_b{q}.scores*inv(rotationmatrix);

            % Finding scores to plot
            for r = 1:size(varargout{q}.scores, 2)
                if q == 1
                    for f = 1:length(unique(timepoints))
                        varargout{q}.scores_boot{r,1}(i,f) = model_b{q}.scores(find(data_X_b.timepoint == timepoints(f),1, 'first'), r);
                    end
                    varargout{q}.loadings_boot{r,1}(i,:) = model_b{q}.loadings(:,r)';
                elseif q > 1
                    for t = 1:length(treatments)
                        for f = 1:length(timepoints)
                            varargout{q}.scores_boot{r,t}(i,f) = model_b{q}.scores(find(data_X_b.timepoint == timepoints(f) & data_X_b.treatment == treatments(t),1, 'first'), r);
                        end
                    end
                    varargout{q}.loadings_boot{r,1}(i,:) = model_b{q}.loadings(:,r)';
                end
            end
        end
    end
end

end