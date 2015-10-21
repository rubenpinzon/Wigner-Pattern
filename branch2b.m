% BRANCH 2b - GPFA analysis of population data for a Linear track with
% simulated cells
% git:
% evernote : 
% The synthetic cells are store in txt files with two colums, time spike
% and cell to which it belongs.


clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/synthetic/db6/';
pattern         = 'spikes_*.txt';
Fs              = 1000;
T_real          = 1.610 * Fs; 

%======Extraxt the spikes for each cell for each lap======================%

D               = dir([basepath pattern]);
n_laps          = length(D);

for ilap = 1 : n_laps 
   spikes   = load([basepath D(ilap).name]); 
   n_cells  = max(spikes(:,2)) + 1;
   n_cells  = 99;
   T(ilap)  = 0; % for saving the lap duration
   
   %data structure for the GPFA
   D(ilap).data       = zeros(n_cells,T_real);
   D(ilap).trialId    = ilap;
   D(ilap).condition  = 'linear Track';
   
   for icell = 1 : n_cells
       s = spikes(spikes(:,2)==icell-1,1);
       spk{icell, ilap} = s;
       if ~isempty(s)
           T(ilap)          = max(T(ilap), max(spk{icell, ilap}));
           D(ilap).data(icell,floor(s*Fs)+2) = 1;
       end
   end   
end

bin_size        = 0.02;  %20 ms
zDim            = 10;    % Target latent dimensions
min_firing      = 0.5;
[D,keep_cell]   = segment(D, bin_size, Fs, min_firing);
showpred        = true; %show the predicted and real firing rates

%% ==================Using DataHigh Library==============================%%
folds           = 3;
mask            = false(1,length(D)); % for cross validation
cv_trials       = randperm(length(D));
fold_indx       = floor(linspace(1,length(D)+1, folds+1));
saveplot        = true;

%DataHigh(D,'DimReduce')

for ifold = 1 : 1%folds  % two-fold cross-validation        
    % prepare masks:
    % test_mask isolates a single fold, train_mask takes the rest
    test_mask       = mask;
    test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
    train_mask = ~test_mask;
    train_data = D(train_mask);
    test_data  = D(test_mask);
    %training of the GPFA
    [params, gpfa_traj, ll_tr] = gpfa_mod(train_data,zDim);

    %Posterior of test data given the trained model
    [traj, ll_te] = exactInferenceWithLL(test_data, params,'getLL',1);
    % orthogonalize the trajectories4
    [Xorth, Corth] = orthogonalize([traj.xsm], params.C);
    traj = segmentByTrial(traj, Xorth, 'data');
% 
    %Validation with LNO
    cv_gpfa_cell = struct2cell(cosmoother_gpfa_viaOrth_fast(test_data,params,zDim));

    true_data      = [test_data.y];
    T              = [0 cumsum([test_data.T])];
    cvdata         = zeros(size(true_data));
    for i = 1 : length(test_data)
       cvdata(:, T(i)+1:T(i+1)) = cv_gpfa_cell{end,i};
    end
    mse_fold        = sum(sum((cvdata-true_data).^2));

    if showpred
        plot_firing(cvdata, true_data, T)    
    
        figure(20)
        set(gcf, 'position', [1983,1,1424,973], 'color', 'w')
        color = jet(length(traj));
        start_traj = []; end_traj = [];
        for ilap = 1 : length(traj)
           lap_t = T(ilap)+1:T(ilap+1);

           plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),Xorth(3,lap_t),[1 2 4 5 7 8],{'X_1','X_2','X_3'},color(ilap,:))           
           plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),[],3,{'X_1','X_2'},color(ilap,:))
           plot_xorth(Xorth(2,lap_t),Xorth(3,lap_t),[],6,{'X_2','X_3'},color(ilap,:))
           plot_xorth(Xorth(1,lap_t),Xorth(3,lap_t),[],9,{'X_1','X_3'},color(ilap,:))          

           start_traj(ilap, :) = Xorth(1:3,lap_t(1));
           end_traj(ilap, :)   = Xorth(1:3,lap_t(end));
        end
        ellipse_eig(end_traj(:,1:2), 3, [1, 0, 0])
        ellipse_eig(end_traj(:,2:3), 6,[1, 0, 0])
        ellipse_eig(end_traj(:,[1,3]), 9,[1, 0, 0])
        ellipse_eig(start_traj(:,1:2), 3, [0, 0, 1])
        ellipse_eig(start_traj(:,2:3), 6,[0, 0, 1])
        ellipse_eig(start_traj(:,[1,3]), 9,[0, 0, 1])
        subplot(3,3,3)
        text(20, -0.2, 'start','color','b')
        text(-0.3, -0.5, 'end','color','r')
        if saveplot
            print(gcf,[basepath sprintf('x_orth_cond(fold%d).png',ifold)],'-dpng')
            title_span(gcf,sprintf('Neural Space (SVD ort1ho)(fold %d)', ifold)); 
        end
    end
    
    mse(ifold)  = mse_fold;
    like(ifold) = ll_te;
    paramsGPFA{ifold} = params;
    fprintf('Trained/validated fold %d\n',ifold)
    clear train_data test_data cvdata cv_gpfa* params
end
result.params = paramsGPFA;
result.mse = mse;
result.like = like;
result.cv_trials = cv_trials;
result.foldidx = fold_indx;

%% ======================Using Aecker Library=============================

folds           = 3;
mask            = false(1,length(D)); % for cross validation
cv_trials       = randperm(length(D));
fold_indx       = floor(linspace(1,length(D)+1, folds+1));
saveplot        = true;


for ifold = 1 : 1%folds  % two-fold cross-validation        
    % prepare masks:
    % test_mask isolates a single fold, train_mask takes the rest
    test_mask       = mask;
    test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
    train_mask = ~test_mask;
    train_data = D(train_mask);
    test_data  = D(test_mask);
    
    Y = reshape([train_data.y],[48, 80, numel(train_data)]);
    model = GPFA('Tolerance', 1e-8, 'Verbose',true);
    model = model.fit(Y, 10, 'hist');

end