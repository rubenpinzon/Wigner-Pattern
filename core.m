%% Analysis of Buzsaki database i01_maze06_MS.002


clc, close all;

basepath        = '/media/bigdata/i01_maze06.002/';
obj             = load([basepath 'i01_maze06_MS.002_BehavElectrData.mat']);
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
[spk, spk_lap]  = get_spikes(clusters, obj.Spike.res, laps);
N               = size(spk_lap,2);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
trialcolor      = jet(5);
% Extract spks when the mouse if running and in the wheel to calculate
% neural trajectories. 

for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    idx_lap = laps(lap):laps(lap+1);

    wheelNonZero = find(wheelspeed(idx_lap)~=0) +  laps(lap);
    if ~isempty(wheelNonZero)
        length_wheel(lap) = numel(wheelNonZero)/Fs;
        %last lap does not have wheel
        for neu=1:N
            idx = spk_lap{lap,neu}>=wheelNonZero(1)...
                          & spk_lap{lap,neu}<=wheelNonZero(end);
            SpkWheel_lap{lap,neu} = spk_lap{lap,neu}(idx);
        end
    end
    
    %(b) Runing in the maze. Extracted based on the
    %EnterSection time stamp without considering left-right
    idx_run = [events{lap}(2,1), sum(events{lap}(5:6,2))];
    int_at_maze(lap, :) = idx_run;
    length_run(lap) = (idx_run(2)-idx_run(1))/Fs;
    %sect 1:enter, 6:exit
    for neu=1:N
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx);
    end
    %Type of trial
    trial{lap}          = typetrial{obj.Laps.TrialType(laps(lap))};
    color(lap,:)        = trialcolor(obj.Laps.TrialType(laps(lap)),:);
end

% Only pyramidal cells
window   = 3.0;
MaxTimeE = window * Fs; %seconds   
% SpkSeries=get_series(data, MaxTimeE);

%Data processed for datahigh without interneuorns
SpkRun_DH = get_high(SpkRun_lap(:,isIntern==0), MaxTimeE,...
                     trial, color);

MaxTimeE = 14 * Fs; %seconds %minimum time at the wheel 14.7224  
SpkWheel_DH = get_high(SpkWheel_lap(:,isIntern==0), MaxTimeE,...
                        trial, color);

%% GUI based DataHigh
%DataHigh(SpkRun_DH, 'DimReduce');

%% Command based GPFA based on DataHigh Library
% remove one neuron at random until having only 15 
% leave two laps to test anc calculate ll_test = P(Ytest|model)

bin_sizes       = 0.1; 
dims            = 2:15; % Target latent dimensions
results(length(bin_sizes)).bin = bin_sizes(end);
test_trials     = [1, 2; 3, 4; 5, 6]; % one left and one right, 3 folds


D               = SpkRun_DH;
firing_thr      = 0.2 ; % Minimum firing rate find which 
                        % neurons should be kept
m               = mean([D.data],2) * Fs;
keep_neurons    = m >= firing_thr;
fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
            sum(keep_neurons),firing_thr)
% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(keep_neurons,:);
end
cells           = sum(keep_neurons);
removeCells     = randperm(cells);
mask            = false(1,length(D));

for cnt = 1: cells - 15
    
    
    % Remove a random neuron 
    for itrial = 1:length(D)
        D(itrial).data(removeCells(cnt),:) = [];
    end
    yDim        = size(D(1).data, 1);
    fprintf('Neuron %d-th removed, %d remaining\n',removeCells(cnt),...
                                                yDim)
    bin             = bin_sizes;    
    bin_width       = ceil(bin * Fs); % bin size (Seconds * Fs) = samples
    useSqrt         = 1; % square root tranform?    
    seq             = [];
    
    
    %Extrat bins for one trial, since all the trials
    %are of the same duration

    T           = floor(size(D(1).data, 2) / bin_width);

    for n = 1 : length(D)
        seq(n).y   = nan(yDim, T);
        for t = 1:T
          iStart = bin_width * (t-1) + 1;
          iEnd   = bin_width * t;
          seq(n).y(:,t) = sum(D(n).data(:, iStart:iEnd), 2);
        end
        %normalization with square root transform
        if useSqrt
            seq(n).y = sqrt(seq(n).y);
        end
    end
    [D.data]    = deal(seq.y);
    %prellocating variables
    lat         = []; % Struct for the latent variables
    ll_test     = 0;  % these store the likelihood
    ll_train    = 0;  % these store the likelihood
    mse_fold    = 0;  % and mse for one fold
    ll_train    = zeros(length(dims),length(test_trials)); % keeps a running total of likelihood
    ll_test     = zeros(length(dims),length(test_trials)); % keeps a running total of likelihood
    mse         = zeros(length(dims),length(test_trials)); % and mse
    % 
    paramsGPFA = cell(length(dims), length(test_trials));
    orth_traje = cell(length(dims), length(test_trials));
    
    %Prepare data for the datahigh library
    [D.y] = D.data;
    for i=1:length(D);
        D(i).trialId = i;
        D(i).T = size(D(i).data,2);
    end
    
    for idim = 1 : length(dims)
        for ifold = 1 : length(test_trials)  % three-fold cross-validation        
            % prepare masks:
            % test_mask isolates a single fold, train_mask takes the rest
            test_mask = mask;
            test_mask(test_trials(ifold, :)) = true;
            
            train_mask = ~test_mask;

            train_data = D(train_mask);
            test_data = D(test_mask);
            %training of the GPFA
            [params, gpfa_traj, ll_train] = gpfa_mod(train_data,dims(idim),...
                                                     'bin_width', bin_width);

            %Posterior of test data given the trained model
            [a, ll_test] = exactInferenceWithLL(test_data, params,'getLL',1);

            %Validation
            cv_gpfa_cell = struct2cell(cosmoother_gpfa_viaOrth_fast...
                                      (test_data,params,idim));
            cvdata = cell2mat(cv_gpfa_cell(7,:));
            mse_fold = sum(sum((cvdata-[test_data.data]).^2));

            mse(idim, ifold) = mse_fold;
            ll_test(idim, ifold) = ll_test;
            ll_train(idim, ifold) = ll_train;
            paramsGPFA{idim, ifold} = params;
            orth_traje{idim, ifold} = gpfa_traj;
        end        
    end
    % best dimension and best across folds
    [a, foldmax] = max(sum(ll_test));
    [a, imax] = max(like(:,foldmax));
    
    results(cnt).like = like;
    results(cnt).mse = mse;
    results(cnt).dim = imax;
    % Setup figure to sumarize results
    cla
    figure(cnt)
    til = sprintf('Run Section With %3.1f ms bins and %2.1fs spike trains', 1000*bin, window);
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', til, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center',...
        'Fontsize',18)
    set(gcf,'color', 'w', 'position', [100 100 1400 700])

    % projection over the three best dimensions (SVD)
    for itrial = 1:length(test_data)
        p = orth_traje{imax, foldmax}(itrial).data;
        c = orth_traje{imax, foldmax}(itrial).epochColors;
        subplot(2,3,[1 2 4 5]), grid on
        plot3(p(1,:), p(2,:),p(3,:), 'Color', c,...
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
    plot(imax, mean_mse(imax),'r*', 'markersize',10)
    xlabel('Latent Dimension')
    ylabel('Mean Sqaure Error (LNO)')

    % LogLike
    subplot(2,3,6)
    mean_like = mean(like,2);
    plot(mean_like,'linewidth',2, 'linestyle','-.'), hold on
    plot(like)
    plot(imax, mean_like(imax),'r*', 'markersize',10)
    xlabel('Latent Dimension')
    ylabel('Log Likelihood')
    namePNG = sprintf('Results/i01_maze06.002_cells_%d', yDim);
    print(gcf,[basepath namePNG],'-dpng')
    savefig(gcf, [basepath namePNG '.fig'])
    
end

save([basepath 'Results/i01_maze06.002_run_w' num2str(window)...
        's_results.mat'],'results')
