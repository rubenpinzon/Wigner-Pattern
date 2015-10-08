%% =======================================================================%
%=================      Preprocessing             ========================%
%=========================================================================%
clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);


%========================Variables of Interest===========================
animal          = 6;
obj             = load(files{animal});
clusters        = obj.Spike.totclu;
laps            = obj.Laps.StartLaps(obj.Laps.StartLaps~=0); %@1250 Hz
%Adding the end of the last lap because Laps only contains the start
laps(end+1)     = obj.Par.SyncOff;
mazesect        = obj.Laps.MazeSection;
wheelspeed      = obj.Laps.WhlSpeedCW;
X               = obj.Track.X;
Y               = obj.Track.Y;
events          = obj.Par.MazeSectEnterLeft;
Fs              = obj.Par.SamplingFrequency;
eeg             = obj.Track.eeg;
speed           = obj.Track.speed;
isIntern        = obj.Clu.isIntern;
numLaps         = numel(laps)-1;
time_seg        = linspace(0, length(eeg)/Fs, length(eeg));
[spk, spk_lap]  = get_spikes(clusters, obj.Spike.res,...
                                                laps);
N               = size(spk_lap,2);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
tcolor          = jet(numLaps+1);
%% =======================================================================%
%==============   Extract Running/wheel Section   ========================%
%=========================================================================%
debug           = 1; %to show diganostic plots

% Extract spks when the mouse is running and in the wheel to calculate
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    idx_lap = laps(lap):laps(lap+1);

    wheelNonZero        = find(wheelspeed(idx_lap)~=0) + laps(lap);
    wheelSpeed_lap{lap} = wheelspeed(wheelNonZero);
    if ~isempty(wheelNonZero)
        wheel_len(lap)   = numel(wheelNonZero)/Fs;

        if debug
            figure(1)
            x   = linspace(0, wheel_len(lap), numel(wheelSpeed_lap{lap})); 
            plot(x, wheelSpeed_lap{lap}, 'color', tcolor(lap,:)), hold on
            text(10, 100 + 20*lap, ['Lap ' num2str(lap)], 'color', tcolor(lap,:));
        end
        
        %last lap does not have wheel
        for neu=1:N
            idx = spk_lap{lap,neu}>=wheelNonZero(1)...
                          & spk_lap{lap,neu}<=wheelNonZero(end);
            %aligned to the start of the section
            SpkWheel_lap{lap,neu} = spk_lap{lap,neu}(idx) - wheelNonZero(1) + 1;
        end
    end
    
    %(b) Runing in the maze. Extracted based on the
    %EnterSection time stamp without considering left-right
    idx_run                 = [events{lap}(2,1), sum(events{lap}(5:6,2))];
    int_at_maze(lap, :)     = idx_run;
    run_len(lap)            = (idx_run(2)-idx_run(1))/Fs;
    X_at_maze               = X(idx_run(1):idx_run(2));
    Y_at_maze               = Y(idx_run(1):idx_run(2));
    speed_lap               = speed(idx_run(1):idx_run(2));

    %sect 1:enter, 6:exit
    for neu=1:N
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        %aligned to the start of the section
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx) - idx_run(1) + 1;
    end
    if debug
        figure(2)
        plot(X_at_maze, Y_at_maze, 'color', tcolor(lap,:)), hold on
        text(700, 400 + 30*lap, ['Lap ' num2str(lap)], 'color', tcolor(lap,:));
        figure(22)
        plot(speed_lap, 'color', tcolor(lap,:)), hold on
        text(700, 800 + 150*lap, ['Lap ' num2str(lap)], 'color', tcolor(lap,:));
    end
    %Type of trial
    trial{lap}          = typetrial{obj.Laps.TrialType(laps(lap))};
    color(lap,:)        = tcolor(obj.Laps.TrialType(laps(lap)),:);
end
figure(1), title('Wheel speed per Lap')
figure(2), title('Position of animal per Lap in section Run')
if debug
   figure(3)
   plot(1:numLaps,run_len,'-s'), hold on
   plot(1:numLaps-1,wheel_len, '-o')
   xlabel('Lap Num.'), ylabel('time (s)')
   title('Duration of sections per lap')
   legend('Run', 'Wheel')
end
figure(22)
title('Linear distance spanned in the wheel')
%% =======================================================================%
%==============   Extract Running/wheel Section   ========================%
%=========================================================================%
% Convert to DataHigh format without segmenting, that is, using the whole
% time that the animal spent in the runing section. This implies laps with
% different lenghts

MaxTimeE = floor(Fs * run_len);
onlyCorrectTrial = true;
%Data processed for datahigh without interneuorns
SpkRun_DH = get_high(SpkRun_lap(:,isIntern==0), MaxTimeE,...
                     trial, color, 'run', onlyCorrectTrial);

if debug
    figure(4)
    for ilap = 1: numLaps
        if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'left')
            idx_run                 = [events{ilap}(2,1), events{ilap}(2,1) + MaxTimeE(ilap)];
            X_at_maze               = X(idx_run(1):idx_run(2));
            Y_at_maze               = Y(idx_run(1):idx_run(2));
            plot(X_at_maze, Y_at_maze, 'color', tcolor(ilap,:)), hold on
            plot(X_at_maze(end), Y_at_maze(end), 's', 'color',...
                 tcolor(ilap,:), 'MarkerSize', 8, 'markerfacecolor',tcolor(ilap,:))
            text(700, 400 + 30*ilap, ['Lap ' num2str(ilap)], 'color', tcolor(ilap,:)); 
        end
    end
end
title('Chk: Segmentation length (squares)')                 
%MaxTimeE = floor(Fs * min(wheel_len) * ones(1, numLaps-1));   
MaxTimeE = floor(Fs * wheel_len);
SpkWheel_DH = get_high(SpkWheel_lap(:,isIntern==0), MaxTimeE,...
                        trial, color, 'wheel',onlyCorrectTrial);

%% GUI based DataHigh
DataHigh(SpkRun_DH, 'DimReduce');

%% =======================================================================%
%======   Command based GPFA based on DataHigh Library   =================%
%=========================================================================%

bin_size        = 0.02;  %20 ms
zDim            = 20;    % Target latent dimensions
results(1).bin  = bin_size;

D               = SpkRun_DH;

firing_thr      = 0.5 ; % Minimum firing rate find which 
                        % neurons should be kept
m               = mean([D.data],2) * Fs;
keep_neurons    = m >= firing_thr;
fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
            sum(keep_neurons),firing_thr)
% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(keep_neurons,:);
end
yDim            = sum(keep_neurons);
useSqrt         = 1; % square root tranform for pre-processing?    
                                   
  
bin_width       = ceil(bin_size * Fs); % bin size (Seconds * Fs) = samples

%Extrat bins for one trial, since all the trials
%are of the same duration
for ilap = 1 : length(D)
    seq         = [];
    T           = floor(size(D(ilap).data, 2) / bin_width);
    seq.y       = nan(yDim, T);
    for t = 1:T
      iStart        = bin_width * (t-1) + 1;
      iEnd          = bin_width * t;
      seq.y(:,t)    = sum(D(ilap).data(:, iStart:iEnd), 2);
    end
    %normalization with square root transform
    if useSqrt
        seq.y       = sqrt(seq.y);
    end
    D(ilap).data    = seq.y;
    D(ilap).y       = seq.y;
    D(ilap).T       = T;    
end

%preallocating variables
folds           = 4;
mask            = false(1,length(D)); % for cross validation
test_mask       = mask;
cv_trials       = randperm(length(D));
fold_indx       = floor(linspace(1,length(D)+1, folds + 1));


for ifold = 1 : folds  % three-fold cross-validation        
    % prepare masks:
    % test_mask isolates a single fold, train_mask takes the rest
    test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;

    train_mask = ~test_mask;

    train_data = D(train_mask);
    test_data  = D(test_mask);
    %training of the GPFA
    [params, gpfa_traj, ll_tr] = gpfa_mod(train_data,zDim,...
                                             'bin_width', bin_width);

    %Posterior of test data given the trained model
    [traj, ll_te] = exactInferenceWithLL(test_data, params,'getLL',1);
    % orthogonalize the trajectories
    [Xorth, Corth] = orthogonalize([traj.xsm], params.C);
    traj = segmentByTrial(traj, Xorth, 'data');

    %Validation with LNO
    cv_gpfa_cell = struct2cell(cosmoother_gpfa_viaOrth_fast...
                              (test_data,params,2:2:zDim));
    
    mse_fold = [];                      
    for dim = length(cv_gpfa_cell):-1:7
        cvdata          = cell2mat(cv_gpfa_cell(dim,:));
        mse_fold(end+1) = sum(sum((cvdata-[test_data.data]).^2));
    end
    

    mse(ifold, : ) = mse_fold;
    paramsGPFA{ifold} = params;
    fprintf('Trained/validated fold %d\n',ifold)
    clear train_data test_data
end        

% best dimension and best across folds
[~, foldmax_te] = max(sum(ll));
[~, imax_te] = max(ll(:,foldmax_te));
[~, foldmax_tr] = max(sum(ll_train));
[~, imax_tr] = max(ll(:,foldmax_tr));

results(ibin).ll_test = ll;
results(ibin).mse = mse;
results(ibin).dim = imax_te;
results(ibin).GPFA = paramsGPFA;
results(ibin).traj = orth_traje_te;

% Setup figure to sumarize results
close all
figure(ibin)
til = sprintf('Run Section With %3.1f ms bins and %2.1fs spike trains', 1000*bin_sizes(ibin), window);
annotation('textbox', [0 0.9 1 0.1], ...
    'String', til, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Fontsize',18)
set(gcf,'color', 'w', 'position', [100 100 1400 700])

% projection over the three best dimensions (SVD)
for itrial = 1:length(train_data)
    p = orth_traje_tr{imax_tr, foldmax_tr}(itrial).data;
    c = orth_traje_tr{imax_tr, foldmax_tr}(itrial).epochColors;
    subplot(2,3,[1 2 4 5]), grid on
    plot(p(1,:), p(2,:), 'Color', c,...
          'linewidth',2); hold on
end
for itrial = 1:length(test_data)
    p = orth_traje_te{imax_te, foldmax_te}(itrial).data;
    c = orth_traje_te{imax_te, foldmax_te}(itrial).epochColors;
    subplot(2,3,[1 2 4 5]), grid on
    plot(p(1,:), p(2,:), 'Color', 0.5*c,...
          'linewidth',2); hold on
end
xlabel('Eigenvector 1')
ylabel('Eigenvector 2')
zlabel('Eigenvector 3')

% MSE
subplot(2,3,3)
mean_mse = mean(mse,2);
plot(mean_mse,'linewidth',2, 'linestyle','-.'), hold on
plot(mse)
plot(imax_te, mean_mse(imax_te),'r*', 'markersize',10)
xlabel('Latent Dimension')
ylabel('Mean Sqaure Error (LNO)')

% LogLike
subplot(2,3,6)
mean_ll_test  = mean(ll,2);
mean_ll_train = mean(ll_train,2);
offset_te = mean(mean_ll_test);
offset_tr = mean(mean_ll_train);
plot(dims,mean_ll_test,'linewidth',2, 'linestyle','-.', 'color','k'), hold on
plot(dims,mean_ll_train,'linewidth',2, 'linestyle','-.', 'color','b')
plot(dims,ll,'k')
plot(dims,ll_train,'b')
plot(imax_te, mean_ll_test(imax_te),'k*', 'markersize',10)
plot(imax_tr, mean_ll_train(imax_tr),'bs', 'markersize',10)
xlabel('Latent Dimension')
ylabel('Log Likelihood')
title('Train','color','b')
box off
namePNG = sprintf('Results/i01_maze06.002_cells_%d', bin_sizes(ibin));
print(gcf,[basepath namePNG],'-dpng')

figure(ibin + 1)
numDim = 7; % 9 latent variables
til = sprintf('9 Latent Variables');
annotation('textbox', [0 0.9 1 0.1], ...
    'String', til, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Fontsize',18)
set(gcf,'color', 'w', 'position', [100 100 1400 700])
time = linspace(0, window, T); % from he extraction program
F    = orth_traje_tr{numDim, foldmax_tr};
for l= 1:length(F)
   for v = 1:dims(numDim)
       subplot(3, 3, v)
       plot(time, F(l).data(v, :),'color',F(l).epochColors), hold on       
   end
end
namePNG = sprintf('Results/i01_maze06.002_LV_%d', 1000*bin_sizes(ibin));
print(gcf,[basepath namePNG],'-dpng')    



nameFile = sprintf('Results/i01_maze06.002_results_%d', 1000*bin_sizes(ibin));
save([basepath 'Results/i01_maze06.002_run_w' num2str(window)...
        's_results.mat'],'results')
