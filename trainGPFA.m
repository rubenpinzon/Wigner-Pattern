function M = trainGPFA(D, zDim, showpred, folds)
%TRAINGPFA trains an cross validates a gpfa model with the data in D, using 3 folds
%           fields required in D: y: spike trains
%ruben pinzon 2015

% folds           = 3;

mask            = false(1,length(D)); % for cross validation
cv_trials       = randperm(length(D));
fold_indx       = floor(linspace(1,length(D)+1, folds+1));

mse             = zeros(1,folds);
like            = zeros(1,folds);
like_tr         = cell(1,folds);
paramsGPFA      = cell(1, folds);
test_trials     = cell(1, folds);
train_trials    = cell(1, folds);

for ifold = 1 : folds  % two-fold cross-validation        
    % prepare masks:
    % test_mask isolates a single fold, train_mask takes the rest
    
    test_mask       = mask;
    test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
    train_mask = ~test_mask;
    train_data = D(train_mask);
    test_data  = D(test_mask);
    
    test_trials{ifold}  = [train_data.trialId];
    train_trials{ifold} = [test_data.trialId];
    
    fprintf('training with trials %s\n',sprintf('%d, ',find(train_mask==1)))
    %training of the GPFA
    [params, gpfa_traj, ll_tr] = gpfa_mod(train_data,zDim);

    %Posterior of test data given the trained model
    [traj, ll_te] = exactInferenceWithLL(test_data, params,'getLL',1);
    % orthogonalize the trajectories4
    [Xorth, Corth] = orthogonalize([traj.xsm], params.C);
    traj = segmentByTrial(traj, Xorth, 'data');
% 
    %Validation with LNO
    cv_gpfa_cell = struct2cell(cosmoother_gpfa_viaOrth_fast...
                              (test_data,params,zDim));

    true_data      = [test_data.y];
    T              = [0 cumsum([test_data.T])];
    cvdata         = zeros(size(true_data));
    for i = 1 : length(test_data)
       cvdata(:, T(i)+1:T(i+1)) = cell2mat(cv_gpfa_cell(end,:,i));
    end
    mse_fold        = sum(sum((cvdata-true_data).^2));

    if showpred
       plot_firing(cvdata, true_data, T)            
    end

    mse(ifold)          = mse_fold;
    like(ifold)         = ll_te;
    like_tr{ifold}      = ll_tr;
    paramsGPFA{ifold}   = params;
    
    fprintf('Trained/validated fold %d\n',ifold)
    clear train_data test_data cvdata cv_gpfa* params
end

M.params      = paramsGPFA;
M.mse         = mse;
M.like_test   = like;
M.like_train  = like_tr;
M.testTrials  = test_trials;
M.trainTrials = train_trials;



clear result params* mse like

