% BRANCH 2 - During maze experiments, using a model trained in the theta 
% regime, predict activity in SPW regime keeping all the parameters of the 
% model unchanged but the time constant in the GPs.
% 
% Addendum Gergo: the purpose of this exercise is to be able to extend 
% beyond the conclusions of the linear maze by checking the translation 
% of the theta-regime subspace to the SPW-regime subspace in cases when 
% there are different behavioral contexts. Proposed analyses:
%
% (1), train a GPFA on left-arm trajectories and right-arm trajectories 
% separately and also jointly, then analysing SPWS following different 
% choices made, 
%
% (2), as a control a similar analysis of SPWs preceding specific choices ;
%
% (3), compare a two-D FA
%
% (4), a similar analysis for SPW following wheel running: using either 
%
% GPFA trained on the exploration or during wheel running. Comparing the 
% SPWs following exploration and those following wheel running 


clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);


%========================Variables of Interest===========================
animal          = 6;
data            = load(files{animal});
clusters        = data.Spike.totclu;
laps            = data.Laps.StartLaps(data.Laps.StartLaps~=0); %@1250 Hz
laps(end+1)     = data.Par.SyncOff;
mazesect        = data.Laps.MazeSection;
events          = data.Par.MazeSectEnterLeft;
Fs              = data.Par.SamplingFrequency;
X               = data.Track.X;
Y               = data.Track.Y;
eeg             = data.Track.eeg;
speed           = data.Track.speed;
isIntern        = data.Clu.isIntern;
numLaps         = numel(laps)-1;
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
n_cells         = size(spk_lap,2);
color           = jet(numLaps+1);
conditions      = {'_left', '_right', ''};

%% =======================================================================%
%==============   Extract Running Sections        ========================%
%=========================================================================%
debug           = false; %to show diganostic plots

% Extract spks when the mouse is running and in the wheel to calculate
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    
    idx_run                 = [events{lap}(2,1), sum(events{lap}(5:6,2))];
    int_at_maze(lap, :)     = idx_run;
    run_len(lap)            = (idx_run(2)-idx_run(1))/Fs;
    X_at_maze               = X(idx_run(1):idx_run(2));
    Y_at_maze               = Y(idx_run(1):idx_run(2));
    speed_lap               = speed(idx_run(1):idx_run(2));

    %sect 1:enter, 6:exit
    for neu=1:n_cells
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        %aligned to the start of the section
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx) - idx_run(1) + 1;
    end
    if debug
        figure(2)
        plot(X_at_maze, Y_at_maze, 'color', color(lap,:),...
            'displayname',sprintf('Lap %d',lap))
        hold on       
    end
    %Type of trial
    trial{lap}                = typetrial{data.Laps.TrialType(laps(lap))};
    color_trial(lap,:)        = color(data.Laps.TrialType(laps(lap)),:);
end
if debug
   figure(2), title('Position of animal per Lap in section Run')

   figure(3)
   plot(1:numLaps,run_len,'-s'), hold on
   xlabel('Lap Num.'), ylabel('time (s)')
   title('Duration of sections per lap')
end

% Convert to DataHigh format without segmenting, that is, using the whole
% time that the animal spent in the runing section. This implies laps with
% different lenghts

MaxTimeE            = floor(Fs * run_len);
onlyCorrectTrial    = true;
%Data processed for datahigh without interneuorns
SpkRun_DH           = get_high(SpkRun_lap(:,isIntern==0), MaxTimeE,...
                     trial, color, 'run', onlyCorrectTrial);

%% =======================================================================%
%======== (1) Train on Left/Right arms separately        =================%
%=========================================================================%

bin_size        = 0.04;  %20 ms
zDim            = 10;    % Target latent dimensions
results(1).bin  = bin_size;
min_firing      = 1;
[D,keep_cell]   = segment(SpkRun_DH, bin_size, Fs, min_firing);
[D_left, D_right] = split_trails(D);
showpred        = true; %show the predicted and real firing rates
folds           = 3;
try
    load([roots{animal} '_branch2_results40ms.mat'])
catch
    % train separately left/right/all 
    for s = 1 : 1%length(conditions)

        Data = eval(sprintf('D%s',conditions{s}));
        %preallocating cross validation variables, two folds    

        mask            = false(1,length(Data)); % for cross validation
        cv_trials       = randperm(length(Data));
        fold_indx       = floor(linspace(1,length(Data)+1, folds+1));

        for ifold = 1 : folds  % two-fold cross-validation        
            % prepare masks:
            % test_mask isolates a single fold, train_mask takes the rest
            test_mask       = mask;
            test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
            train_mask = ~test_mask;
            train_data = Data(train_mask);
            test_data  = Data(test_mask);
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
               cvdata(:, T(i)+1:T(i+1)) = cell2mat(cv_gpfa_cell(7,:,i));
            end
            mse_fold        = sum(sum((cvdata-true_data).^2));

            if showpred
               plot_firing(cvdata, true_data, T)            
            end

            mse(ifold)  = mse_fold;
            like(ifold) = ll_tr;
            paramsGPFA{ifold} = params;
            fprintf('Trained/validated fold %d\n',ifold)
            clear train_data test_data cvdata cv_gpfa* params
        end
        result.params = paramsGPFA;
        result.mse = mse;
        result.like = like;
        result.cv_trials = cv_trials;
        result.foldidx = fold_indx;

        eval(sprintf('result_D%s = result',conditions{s}))

        clear result params* mse like

    end
    save([roots{animal} '_branch2_results40ms.mat'],...
        'result_D', 'result_D_left', 'result_D_right', 'D', 'keep_cell')
end

%% %%%%%%%%show orthogonalized latent variables:%%%%%%%%%%%%%%%%%%%%

saveplot = false;

for s = 1 : length(conditions)
    
    Data            = eval(sprintf('D%s',conditions{s}));
    Result          = eval(sprintf('result_D%s;',conditions{s}));
    mask            = false(1,length(Data)); % for cross validation
    cv_trials       = Result.cv_trials;
    fold_indx       = Result.foldidx;
    
    
    for ifold = 1 : folds  % two-fold cross-validation        
        % prepare masks:
        % test_mask isolates a single fold, train_mask takes the rest
        test_mask       = mask;
        test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
        test_data  = Data(test_mask);
        %Posterior of test data given the trained model
        [traj, ll_te] = exactInferenceWithLL(test_data, Result.params{ifold},'getLL',1);
        % orthogonalize the trajectories4
        [Xorth, Corth] = orthogonalize([traj.xsm], Result.params{ifold}.C);
        T              = [0 cumsum([test_data.T])];
        
        figure(10)
        set(gcf, 'position', [1983,1,1424,973], 'color', 'w')

        color = jet(length(traj));
        start_traj = []; end_traj = [];
        for ilap = 1 : length(traj)
           lap_t = T(ilap)+1:T(ilap+1);
           
           plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),Xorth(3,lap_t),[1 2 4 5 7 8],{'X_1','X_2','X_3'},[0.8 0.8 0.8])           
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
        text(-0.8, -0.2, 'start','color','b')
        text(-0.3, -0.5, 'end','color','r')
        if saveplot
            print(gcf,[roots{animal} sprintf('x_orth_cond%s(fold%d).png',conditions{s},ifold)],'-dpng')
            title_span(gcf,sprintf('Neural Space (SVD ort1ho) Condition %s (fold %d)',conditions{s}(2:end), ifold)); 
            
        end
        close gcf
    end    
    title_span(gcf,sprintf('Condition %s (Two folds)',conditions{s}(2:end)));
    
    d_diff = sqrt(sum((Result.params{1}.d - Result.params{2}.d).^2));
    C_diff = sqrt(sum((Result.params{1}.Corth - Result.params{2}.Corth).^2));
    R_diff = sqrt(sum((diag(Result.params{1}.R) - diag(Result.params{2}.R)).^2));    
        
end


%% =======================================================================%
%======== (1) Get SPWs in the reward or wheel area       =================%
%=========================================================================%

markers = 1.25*load([roots{animal} '_spws.txt']); %from ms to samples
updateGP = true;

figure(20)
set(gcf, 'position', [0 1 1000 300], 'color', 'w')
plot(eeg./max(eeg),'color',[0.8, 0.8, 0.8]), hold on
plot(0.08*mazesect-0.5,'b')
ylim([-0.6 0.6])
plot(repmat(markers(:,1),1,2),ylim,'r')
plot(repmat(markers(:,2),1,2),ylim,'c')
plot(0.1*data.Laps.TrialType,'r')

spw_tag = {'after_left', 'after_right', 'after_left_error', 'after_right_error'};
color   = jet(4);
for sp = 1 : length(markers)   
    spw_type{sp}            = spw_tag{data.Laps.TrialType(ceil(markers(sp,1)))+1};
    color_spw(sp,:)         = color(data.Laps.TrialType(ceil(markers(sp,1)))+1,:);
    %spikes during the SPW
    
    idx_run                 = ceil(markers(sp,1):markers(sp,2));
    spw_len(sp)            = (idx_run(end)-idx_run(1))/Fs;
    cnt = 1;
    %sect 1:enter, 6:exit
    for neu=1:n_cells
        if isIntern(neu)==0
            idx = spk{neu}>=idx_run(1) & spk{neu}<=idx_run(end);
            %aligned to the start of the section
            SpkSPW{sp,cnt} = spk{neu}(idx) - idx_run(1) + 1;
            cnt = cnt +  1;
            plot(repmat(spk{neu}(idx),1,2)', repmat(0.01*[neu neu+0.9],length(spk{neu}(idx)),1)', 'linewidth',2)
        end
    end
end
ylim([0 1.2])
%Parameters to used from trained models
params          = result_D_left.params{ifold};

SpkSPW_DH       = get_high(SpkSPW(:,keep_cell==1), ceil(spw_len*Fs),...
                     spw_type, color_spw, 'spw', 0);

D               = segment(SpkSPW_DH, 0.004, Fs, 0.01); %2ms bin size
D               = filter_condition(D, spw_tag{1}, 2);
%D               = D(20:30);
[traj, ll_te]   = exactInferenceWithLL(D, params,'getLL',1);
[Xorth, Corth] = orthogonalize([traj.xsm], params.C);

T              = [0 cumsum([D.T])];

%retrained the GPs
if updateGP
    [paramsUP,seq,ll]      = update_gps(params, D, 150);
    [traj, ll_te]  = exactInferenceWithLL(D, paramsUP,'getLL',1);
    [Xorth, Corth] = orthogonalize([traj.xsm], paramsUP.C);
end

%firing rates
if showpred
%Validation with LNO
    cv_gpfa_cell = struct2cell(cosmoother_gpfa_viaOrth_fast...
                              (D,paramsUP,zDim));

    true_SPW       = [D.y];
    T              = [0 cumsum([D.T])];
    cvspw          = zeros(size(true_SPW));
    for i = 1 : length(D)
       cvspw(:, T(i)+1:T(i+1)) = cell2mat(cv_gpfa_cell(7,:,i));
    end

   plot_firing(cvspw, true_SPW, T)            
end

        

figure(1)
set(gcf, 'position', [1983,1,1424,973], 'color', 'w')
color = jet(length(traj));
for ilap = 1 : length(traj)
   lap_t = T(ilap)+1:T(ilap+1);
   plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),Xorth(3,lap_t),[1 2 4 5 7 8],{'X_1','X_2','X_3'},color(ilap,:))           

%    plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),[],3,{'X_1','X_2'},color(ilap,:))
%    plot_xorth(Xorth(2,lap_t),Xorth(3,lap_t),[],6,{'X_2','X_3'},color(ilap,:))
%    plot_xorth(Xorth(1,lap_t),Xorth(3,lap_t),[],9,{'X_1','X_3'},color(ilap,:))          

end

title_span(gcf,sprintf('Neural Space (SVD ort1ho) Condition %s (SPW) 2 ms bin','after_left'));        
print(gcf,[roots{animal} sprintf('x_orth_cond%s(SPW)2ms.png','after_left')],'-dpng')