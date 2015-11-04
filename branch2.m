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
color           = jet(55);
conditions      = {'_left', '_right', ''};
middle_arm      = true;
%% =======================================================================%
%==============   Extract Running Sections        ========================%
%=========================================================================%
debug           = true; %to show diganostic plots
%this is to remove/add the section in the middle arm of the maze
sect            = [3, 4];   
if middle_arm
   sect = 2; 
end
% Extract spks when the mouse is running and in the wheel to calculate
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    
    idx_run                 = [sum(events{lap}(sect,1)), sum(events{lap}(5:6,2))];
    int_at_maze(lap, :)     = idx_run;
    run_len(lap)            = (idx_run(2)-idx_run(1))/Fs;
    X_lap{lap}          = X(idx_run(1):idx_run(2));
    Y_lap{lap}          = Y(idx_run(1):idx_run(2));
    speed_lap               = speed(idx_run(1):idx_run(2));

    %sect 1:enter, 6:exit
    for neu=1:n_cells
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        %aligned to the start of the section
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx) - idx_run(1) + 1;
    end
    if debug
        figure(2)
        plot(X_lap{lap}, Y_lap{lap}, 'color', color(lap,:),...
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
                     trial, color_trial, 'run', onlyCorrectTrial);

%% =======================================================================%
%======== (1) Train on Left/Right arms separately        =================%
%=========================================================================%
name = '_branch2_noMidArm.mat';
if middle_arm
    name = '_branch2_results40ms.mat';
end

bin_size        = 0.04;  %20 ms
zDim            = 10;    % Target latent dimensions
results(1).bin  = bin_size;
min_firing      = 1.0;
[D,keep_cell]   = segment(SpkRun_DH, bin_size, Fs, min_firing);
[D_left, D_right] = split_trails(D);
showpred        = true; %show the predicted and real firing rates
folds           = 3;
try
    load([roots{animal} name])
catch
    % train separately left/right/all 
    for s = 1 : length(conditions)

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

        eval(sprintf('result_D%s = result',conditions{s}))

        clear result params* mse like

    end
    save([roots{animal} name],...
        'result_D', 'result_D_left', 'result_D_right', 'D', 'keep_cell', 'X_lap','Y_lap')
end

%% %%%%%%%%show orthogonalized latent variables:%%%%%%%%%%%%%%%%%%%%

saveplot = true;

for s = 1 : length(conditions)
    
    Data            = eval(sprintf('D%s',conditions{s}));
    Result          = eval(sprintf('result_D%s;',conditions{s}));
    mask            = false(1,length(Data)); % for cross validation
    cv_trials       = Result.cv_trials;
    fold_indx       = Result.foldidx;
    fprintf('Condition %s loaded\n',conditions{s})
    
    for ifold = 1 : folds  % two-fold cross-validation        
        % prepare masks:
        % test_mask isolates a single fold, train_mask takes the rest
        test_mask       = mask;
        test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
        test_data  = Data(test_mask);
        %Posterior of test data given the trained model
        [traj, ll_te] = exactInferenceWithLL(test_data, Result.params{ifold},'getLL',1);
        % orthogonalize the trajectories4
        [Xorth, Corth, TT, EE] = orthogonalize([traj.xsm], Result.params{ifold}.C);
        T              = [0 cumsum([test_data.T])];
        
        figure(10)
        set(gcf, 'position', [1,1,1424,973], 'color', 'w')
        figure(11)
         set(gcf, 'position', [1,1,1424,973], 'color', 'w')

        color = jet(length(traj));
        start_traj = []; end_traj = [];
        for ilap = 1 : length(traj)
           lap_t = T(ilap)+1:T(ilap+1);
           c = color(ilap,:); 
           if s == 3
               c = [0 0 1];
               if strcmp(test_data(ilap).condition,'right')
                    c = [1 0 0];
               end
           end            

           figure(10)
           plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),Xorth(3,lap_t),[1 2 4 5 7 8],{'X_1','X_2','X_3'},c,num2str(test_data(ilap).trialId))           
           plot_xorth(Xorth(1,lap_t),Xorth(2,lap_t),[],3,{'X_1','X_2'},c)
           plot_xorth(Xorth(2,lap_t),Xorth(3,lap_t),[],6,{'X_2','X_3'},c)
           plot_xorth(Xorth(1,lap_t),Xorth(3,lap_t),[],9,{'X_1','X_3'},c)          
           
           start_traj(ilap, :) = Xorth(1:3,lap_t(1));
           end_traj(ilap, :)   = Xorth(1:3,lap_t(end));
           figure(11)
           plot_latent(Xorth(:,lap_t), c)
        end
        figure(10)
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
            figure(10)
            title_span(gcf,sprintf('Neural Space (SVD ort1ho) Condition %s (fold %d) noMiddleArm',conditions{s}(2:end), ifold)); 
            print(gcf,[roots{animal} sprintf('x_orth_cond%s(fold%d)_noMidArm.png',conditions{s},ifold)],'-dpng')

            figure(11)
            title_span(gcf,sprintf('Latents Condition %s (fold %d) noMiddleArm',conditions{s}(2:end), ifold));
            print(gcf,[roots{animal} sprintf('Latents_%s(fold%d)_noMidArm.png',conditions{s},ifold)],'-dpng')

            
        end
                
        close gcf
    end    
    title_span(gcf,sprintf('Condition %s (Two folds)',conditions{s}(2:end)));
    
    d_diff = sqrt(sum((Result.params{1}.d - Result.params{2}.d).^2));
    C_diff = sqrt(sum((Result.params{1}.Corth - Result.params{2}.Corth).^2));
    R_diff = sqrt(sum((diag(Result.params{1}.R) - diag(Result.params{2}.R)).^2));    
        
end

%% =====Effects of changing the kernel lenght on the latent variables=====

s_gamma = [0.1 0.5 1 2 10];
c           = lines(length(s_gamma));
figure(11)
set(gcf, 'position', [1,1,1424,400], 'color', 'w')
for i = 1 :length(s_gamma)
    s               = 1;
    Data            = eval(sprintf('D%s',conditions{s}));
    Result          = eval(sprintf('result_D%s;',conditions{s}));
    Result.params{ifold}.gamma = s_gamma(i)*Result.params{ifold}.gamma;
    mask            = false(1,length(Data)); % for cross validation
    cv_trials       = Result.cv_trials;
    fold_indx       = Result.foldidx;
    ifold           = 1;
    ilap            = 1;

    test_mask       = mask;
    test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
    test_data  = Data(test_mask);
    %Posterior of test data given the trained model
    [traj, ll_te] = exactInferenceWithLL(test_data, Result.params{ifold},'getLL',1);
    % orthogonalize the trajectories4
    [Xorth, Corth, TT, EE] = orthogonalize([traj.xsm], Result.params{ifold}.C);
    T              = [0 cumsum([test_data.T])];
    

    color = jet(length(traj));
    start_traj = []; end_traj = [];
    lap_t = T(ilap)+1:T(ilap+1); 
    plot_latent(Xorth(:,lap_t), c(i,:), sprintf('mult. length %2.2f',s_gamma(i)))
    
end


%% =======================================================================%
%======== (2) Get SPWs in the reward or wheel area       =================%
%=========================================================================%
D           = SpkRun_DH;
markers     = 1.25*load([roots{animal} '_spws.txt']); %from ms to samples
updateGP    = true;
removeInh   = true;

figure(20), hold on
set(gcf, 'position', [0 1 1000 300], 'color', 'w')
%plot(eeg./max(eeg),'color',[0.8, 0.8, 0.8]), hold on
plot(0.08*mazesect-0.5,'-.k')
ylim([-0.6 0.6])
plot(repmat(markers(:,1),1,2),ylim,'b')
plot(repmat(markers(:,2),1,2),ylim,'c')
plot(0.1*data.Laps.TrialType,'r','linewidth',2)

%selects those SPWs that are close to the runs and extracts spks
%saves struct S

run_end = int_at_maze(:,2);
max_d   = 4 * Fs; 
min_d   = 0 * Fs; 
cnt     = 1;
try
    clear S;
catch
   disp('Data cleaned') 
end

for sp = 1 : length(markers)
    d_spw2run   = run_end - markers(sp,1)*ones(length(run_end),1);
    lap_spw     = find(abs(d_spw2run(d_spw2run<0)) <= max_d & abs(d_spw2run(d_spw2run<0)) >= min_d);
    if ~isempty(lap_spw)
        S(cnt).marker = markers(sp,:);
        S(cnt).lap_It_Belongs = lap_spw;
        S(cnt).lap_type = data.Laps.TrialType(ceil(markers(sp,1)));
        S(cnt).name_lap_It_Belongs = typetrial{S(cnt).lap_type};
        S(cnt).markerId = sp;
        %extract spikes
        idx_run               = ceil(markers(sp,1):markers(sp,2));
        S(cnt).duration       = (idx_run(end)-idx_run(1));
        c_cnt                 = 1;
        S(cnt).spk_train      = zeros(n_cells, S(cnt).duration);
        plot(repmat(markers(sp,1),1,2),ylim,'m')
        
        for neu=1:n_cells
            idx = spk{neu}>=idx_run(1) & spk{neu}<=idx_run(end);
            %aligned to the start of the section
            spk_spw        = spk{neu}(idx) - idx_run(1) + 1;
            if ~isempty(spk_spw)
               S(cnt).spk_train(neu, spk_spw) = 1; 
            end
            c_cnt = c_cnt + 1;
        end
        cnt = cnt + 1;
    end 
    
end
n_spws = length(S);

if removeInh
    for sp = 1 : n_spws
        S(sp).spk_train(isIntern==1,:) = [];
    end 
end

%compare the firing rate in the run for each cell with the SPW
for sp = 1 : n_spws
    sct_spw(sp,:) = sum(S(sp).spk_train,2);
    sct_run(sp,:) = sum(SpkRun_DH(S(sp).lap_It_Belongs).data,2);
end

%Spikes from left alts. Since there are no SPWs after L alt. they have to
%be extracted directly from D
cnt_r = 1;
cnt_l = 1;
laps_left = [];
laps_right = [];
clear scn*
for lap = 1 : length(D)
   if strcmp(D(lap).condition, 'left')
       scn_left(:,cnt_l) = sum(D(lap).data,2);
       laps_left(end+1)  = D(lap).trialId;
       cnt_l = cnt_l + 1;
   elseif strcmp(D(lap).condition, 'right')
       scn_right(:,cnt_r) = sum(D(lap).data,2);
       laps_right(end+1)  = D(lap).trialId;

       cnt_r = cnt_r + 1;
   else
       D(lap).condition
   end    
end

%Feature scaling, min = zero
Th_left         = mean(scn_left,2);
Th_left_sc      = Th_left./max(Th_left);
Th_right        = mean(scn_right,2);
Th_right_sc     = Th_right./max(Th_right);
Th_spw_right    = mean(sct_run([S.lap_type]==1,:))';
Th_spw_right_sc = Th_spw_right./max(Th_spw_right);
Th_spw_left    = mean(sct_run([S.lap_type]==2,:))';
Th_spw_left_sc = Th_spw_left./max(Th_spw_left);

SPW_right       = mean(sct_spw([S.lap_type]==1,:))';
SPW_right_sc    = SPW_right./max(SPW_right);
SPW_left        = mean(sct_spw([S.lap_type]==2,:))';
SPW_left_sc     = SPW_left./max(SPW_left);

plot(Th_left_sc), hold on
plot(Th_right_sc)
plot(SPW_left_sc, '--')
plot(SPW_right_sc, '--')

%=========================================================================%
%====================      Correlations    ===============================%
%=========================================================================%
% Theta Right vs Theta Left
% Theta Right vs SPW right
% Theta Left vs SPW right
% SPW LeftE vs SPW right
%%
X = {'Th_right_sc', 'Th_left_sc', 'Th_right_sc', 'Th_spw_right_sc',...
    'Th_spw_right_sc', 'Th_spw_left_sc', 'SPW_right_sc',...
    'Th_spw_right_sc','Th_spw_left_sc',};

Y = {'Th_spw_right_sc', 'Th_spw_left_sc', 'Th_left_sc', 'Th_spw_left_sc',...
    'SPW_right_sc', 'SPW_left_sc', 'SPW_left_sc',...
    'SPW_left_sc','SPW_right_sc'};

%(Right theta vs Error Left, i.e., right)
for cv = 1 : length(X)
    x = eval(X{cv});
    y = eval(Y{cv});

    figure(cv), hold on
    set(gcf, 'position', [1000 100 483 300], 'color', 'w')
    color = parula(length(x));

    for k = 1 : length(x)
        plot(x(k),y(k), 'ok',...
            'Markerfacecolor',color(k,:))
    end
    %robtus fit
    brob = robustfit(x,y);
    y_reg = brob(1)+brob(2)*x;
    plot(x,y_reg,'r-')
    R = 1 - sum((y - y_reg).^2)/(sum((y - mean(y)).^2));
    title(sprintf('R^2 = %1.4f => %3.3f+%3.3fx ',R,brob(1),brob(2)))

    hold on, grid on
    axis([0 1 0 1])
    xlabel(strrep(X{cv},'_',' '))
    ylabel(strrep(Y{cv},'_',' '))
    set(gca,'FontSize',12)
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 6.25])
    
    handaxes2 = axes('Position', [0.6 0.3 0.3 0.3]);
    hold on
    for k = 1 : length(x)
        plot(x(k),y(k), 'ok',...
            'Markerfacecolor',color(k,:))
    end
    plot(x,y_reg,'r-')
    axis([0 0.1 0 0.1])
    set(handaxes2, 'Box', 'off','fontsize',9)
    
    %print(gcf,[roots{animal} sprintf('Corr_%s_%s.png',X{cv},Y{cv})],'-dpng')
    
    clear x y brob
end
%(Right theta vs right spw)

%% =======================================================================%
%========            (3) Test model on SPWs              =================%
%=========================================================================%

params          = result_D_left.params{ifold};
 
 SpkSPW_DH       = get_high(SpkSPW(:,keep_cell==1), ceil(spw_len*Fs),...
                      spw_type, color_spw, 'spw', 0);
 
 D               = segment(SpkSPW_DH, 0.004, Fs, 0.01); %2ms bin size
 D               = filter_condition(D, spw_tag{1}, 2);
 %D               = D(20:30);
 [traj, ll_te]   = exactInferenceWithLL(D, params,'getLL',1);
 [Xorth, Corth]  = orthogonalize([traj.xsm], params.C);
 
 T              = [0 cumsum([D.T])];
 
 %retrained the GPs
 if updateGP
     [paramsUP,seq,ll]      = update_gps(params, D, 150);
     [traj, ll_te]  = exactInferenceWithLL(D, paramsUP,'getLL',1);
     [Xorth, Corth] = orthogonalize([traj.xsm], paramsUP.C);
 end

