% BRANCH 2b - GPFA analysis of The synthetic database created with artificial_ca1.py
%
%            DESCRIPTION:
%            (1) Sets working variables.
%            (2) Extracts the spikes for each lap in the files produced by the script artificial_ca1.py
%            . It also wraps up the data into a datahigh struct for easing the processing with the GPFA
%            library. (3) Trains and validates the GPFA model using cross validation. (4) If spws were
%            generated in the synthetic database, here the P(spw|model) is computed.
%
%            USAGE:
%            set the variable basepath to the folder where the files produced by the script
%            artificial_ca1.py reside. Set the number of folds, latent dimension (zDim) for
%            minimum firing rate (min-firing) for the cross validation during
%            the training of the GPFA and then you are good to go.
%
%            see also artificial_ca1.py, branch2.m, branch2_cleaned.m
%
%Version 1.0 Ruben Pinzon@2015

clc, close all; clear all;

%=============(1) Variables  ===============================

basepath        = '/media/bigdata/synthetic/db11/';
pattern         = 'spwspikes_*.txt';
description     = 'Training the gpfa on spw';
Fs              = 1000;
T_real          = 1.1 * Fs; 
pattern_rates   = 'spw_rates_*.txt';
bin_size        = 0.008;  %80 ms
zDim            = 10;    % Target latent dimensions
min_firing      = 0.1;
pattern         = 'spwspikes_*.txt';
T_real          = 0.1 * Fs;
bin_size_spw    = 0.002;  %2 ms
min_firing_spw  = 0;
folds           = 3;
saveplot        = true;

%======(2) Extract the spikes for each cell for each lap ======================%

D               = dir([basepath pattern]);
n_laps          = length(D);

for ilap = 1 : n_laps 
   spikes   = load([basepath D(ilap).name]); 
%    rates    = load([basepath strrep(D(ilap).name,'spikes','rates')])';
   
   %n_cells  = max(spikes(:,2)) + 1;
   n_cells   = 90;
   T(ilap)   = 0; % for saving the lap duration
   
   %data structure for the GPFA
   D(ilap).data       = zeros(n_cells,T_real);
   D(ilap).trialId    = ilap;
   D(ilap).condition  = 'linear Track';
%    D(ilap).rates      = rates;
   for icell = 1 : n_cells
       s = spikes(spikes(:,2)==icell-1,1);
       spk{icell, ilap} = s;
       if ~isempty(s)
           T(ilap)          = max(T(ilap), max(spk{icell, ilap}));
           D(ilap).data(icell,floor(s*Fs)+2) = 1;
       end
   end
   F(ilap).data = D(ilap).data;
   F(ilap).condition = 'l_track';
end


[D,keep_cell]   = segment(D, bin_size, Fs, min_firing);
showpred        = true; %show the predicted and real firing rates
%DataHigh(F,'DimReduce'); %Use the DataHigh library to train the GPFA

%% ====== (3) Train GPFA Using DataHigh Library ==============================%%

mask            = false(1,length(D)); % for cross validation
cv_trials       = randperm(length(D));
fold_indx       = floor(linspace(1,length(D)+1, folds+1));

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
        plot_firing(cvdata, true_data, T);    
    
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
            %title_span(gcf,description); 
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

%========        (4) LogLike p(spw|model)                =================%

S               = dir([basepath pattern]);
n_laps          = length(S);

for ilap = 1 : n_laps 
   spikes   = load([basepath S(ilap).name]); 
   n_cells  = sum(keep_cell);
   %n_cells  = 90;
  
   %data structure for the GPFA
   S(ilap).data       = zeros(n_cells,T_real);
   S(ilap).trialId    = ilap;
   S(ilap).condition  = 'linear Track';
   
   for icell = 1 : length(keep_cell)
       if keep_cell(icell)
           s = spikes(spikes(:,2)==icell-1,1);
           spk{icell, ilap} = s;
           if ~isempty(s)
               T(ilap)          = max(T(ilap), max(spk{icell, ilap}));
               S(ilap).data(icell,floor(s*Fs)+2) = 1;
           end
       end
   end
   G(ilap).data = S(ilap).data;
   G(ilap).condition = 'l_track';
end

S               = segment(S, bin_size_spw, Fs, min_firing);

%log like
[traj_spw, ll_te_spw] = exactInferenceWithLL(S, paramsGPFA{1},'getLL',1);
[Xorth_spw, Corth_spw]  = orthogonalize([traj_spw.xsm], paramsGPFA{1}.C);
T              = [0 cumsum([S.T])];

figure(20)
for ilap = 1 : length(traj_spw)
   lap_t = T(ilap)+1:T(ilap+1);
   plot_xorth(Xorth_spw(1,lap_t),Xorth_spw(2,lap_t),Xorth_spw(3,lap_t),...
       [1 2 4 5 7 8],{'X_1','X_2','X_3'},0.5*color(ilap,:))  
end

