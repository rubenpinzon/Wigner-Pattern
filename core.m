%% Analysis of Buzsaki database i01_maze06_MS.002


clc, close all; clear all;

basepath        = '/media/bigdata/i01_maze06.002/';
animal          = 'i01_maze06_MS.002';
obj             = load([basepath animal '_BehavElectrData.mat']);
clusters        = obj.Spike.totclu;
laps            = obj.Laps.StartLaps(obj.Laps.StartLaps~=0); %@1250 Hz
%Adding the end of the last lap because Laps only contains the start
laps(end+1)     = obj.Par.SyncOff;
mazesect        = obj.Laps.MazeSection;
wheelspeed      = obj.Laps.WhlSpeedCW;
XT              = obj.Track.X;
YT              = obj.Track.Y;
X               = obj.Spike.X;
Y               = obj.Spike.Y;
events          = obj.Par.MazeSectEnterLeft;
Fs              = obj.Par.SamplingFrequency;
eeg             = obj.Track.eeg;
speed           = obj.Track.speed;
isIntern        = obj.Clu.isIntern;
numLaps         = numel(laps)-1;
time_seg        = linspace(0, length(eeg)/Fs, length(eeg));
[spk, spk_lap, X, Y, X_lap, Y_lap]  = get_spikes(clusters, obj.Spike.res, laps, X, Y);

N               = size(spk_lap,2);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
trialcolor      = hot(5);
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
    idx_run = [events{lap}(2,1), sum(events{lap}(7:8,1))];
    int_at_maze(lap, :) = idx_run;    
    length_run(lap) = (idx_run(2)-idx_run(1))/Fs;
    %sect 1:enter, 6:exit
    for neu=1:N
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx) - idx_run(1);
        X_Run_Lap{lap,neu}  = X_lap{lap, neu}(idx);
        Y_Run_Lap{lap,neu}  = Y_lap{lap, neu}(idx);
    end
    %Type of trial
    trial{lap}          = typetrial{obj.Laps.TrialType(laps(lap))};
    color(lap,:)        = trialcolor(obj.Laps.TrialType(laps(lap)),:);
end
%
%Segment base on spatial coordinates rather than time.
%interpolate the position to the longest time
%the two arms of the T are mirrowed. There are in total 25 grids fo 100x100

longest = max(length_run)*Fs;
Starpos = [XT(int_at_maze(:,1)), YT(int_at_maze(:,1))];
Endpos  = [XT(int_at_maze(:,2)), YT(int_at_maze(:,2))];
xl      = zeros(1,longest) ;yl = zeros(1,longest); %vectors to calculate the mean trajectory
xr      = zeros(1,longest) ;yr = zeros(1,longest); 
r       = 0;l = 0;
x       = []; y = [];

for ilap = 1 : numLaps
    time    = linspace(0,length_run(ilap),longest); 
    %real animal position and interporaltion
    x       = XT(int_at_maze(ilap,1):int_at_maze(ilap,2));
    y       = YT(int_at_maze(ilap,1):int_at_maze(ilap,2));   
    xi      = spline(linspace(0,length_run(ilap),length(x)),x,time);
    yi      = spline(linspace(0,length_run(ilap),length(y)),y,time);        
    if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'errorLeft')
        xr = xr + xi; yr = yr + yi; r  = r + 1;
        %count spikes in the grids of the right arm        
    else
        xl  = xl + xi;yl = yl + yi; l  = l + 1;       
    end
    plot(x, y), hold on
    for icell = 1 : N
        plot(X_Run_Lap{ilap,icell}, Y_Run_Lap{ilap,icell},'x')
    end

end
leftT= [xl; yl]'./l;
rightT= [xr; yr]'./r;

plot(leftT(:,1), leftT(:,2),'k','linewidth', 2), hold on
plot(rightT(:,1), rightT(:,2),'k','linewidth', 2), 
title(sprintf('Animal %s',animal))
xlabel('x')
%script to extract the grids
segments     = 100;
roiDims      = [30 100]; %width and length of ROI
connectgrids = 1;
ctrNeuron    = 30;
show = 0;
gridsL = get_grids(leftT, segments, connectgrids, show, roiDims);
gridsR = get_grids(rightT, segments, connectgrids, show, roiDims);

%count spikes inside grids

count = zeros(N, 43, numLaps); %TODO find the 43 analitically
for ilap = 1 : numLaps
   for icell = 1 : N
       show = 0;
       if icell == ctrNeuron; show = 1; figure(ilap), hold on;end
       x = X_Run_Lap{ilap,icell};
       y = Y_Run_Lap{ilap,icell};
       if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'errorLeft')
          count(icell, :, ilap) = countROIspks(x, y, gridsR, show);
       else
          count(icell, :, ilap) = countROIspks(x, y, gridsL, show); 
       end
       fprintf('Counting spikes Lap %d, neuron %d\n',ilap, icell)
   end    
end

X_pyr = count(~isIntern,:,:);
%%
% Command based GPFA based on DataHigh Library
%DataHigh(D,'DimReduce')
bin_width = 30 ; %30mm

for ilap = 1 : numLaps
    D(ilap ).data = X_pyr(:,:,ilap);
    D(ilap ).condition = trial{ilap}; 
    D(ilap ).epochColors = color(ilap, :);
end

dims            = 3:20; % Target latent dimensions
test_trials     = [1:4; 5:8; 12:15]; % one left and one right, 3 folds

firing_thr      = 0.01 ; % Minimum firing rate find which 
                        % neurons should be kept
m               = mean([D.data],2);
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
yDim            = size(D(1).data, 1);
useSqrt         = 1; % square root tranform?    

%prellocating variables
lat         = []; % Struct for the latent variables
ll_te       = 0;  % these store the likelihood
ll_tr       = 0;  % these store the likelihood
mse_fold    = 0;  % and mse for one fold
folds       = size(test_trials,1);
ll_train    = zeros(length(dims),folds); % keeps a running total of likelihood
ll_test     = zeros(length(dims),folds); % keeps a running total of likelihood
mse         = zeros(length(dims),folds); % and mse    
paramsGPFA  = cell(length(dims), folds);
orth_traje_tr  = cell(length(dims), folds); %training orth_traje
orth_traje_te  = cell(length(dims), folds); %training orth_traje

%Prepare data for the datahigh library
[D.y] = D.data;
for i=1:length(D);
    D(i).trialId = i;
    D(i).T = size(D(i).data,2);
end

for idim = 1 : length(dims)
    for ifold = 1 : folds  % three-fold cross-validation        
        % prepare masks:
        % test_mask isolates a single fold, train_mask takes the rest
        test_mask = mask;
        test_mask(test_trials(ifold, :)) = true;

        train_mask = ~test_mask;

        train_data = D(train_mask);
        test_data = D(test_mask);
        %training of the GPFA
        [params, gpfa_traj, ll_tr] = gpfa_mod(train_data,dims(idim),...
                                                 'bin_width', bin_width);

        %Posterior of test data given the trained model
        [traj, ll_te] = exactInferenceWithLL(test_data, params,'getLL',1);
        % orthogonalize the trajectories
        [Xorth, Corth] = orthogonalize([traj.xsm], params.C);
        traj = segmentByTrial(traj, Xorth, 'data');
        traj = rmfield(traj, {'Vsm', 'VsmGP', 'xsm'});

        %Validation with LNO
        cv_gpfa_cell = struct2cell(cosmoother_gpfa_viaOrth_fast...
                                  (test_data,params,idim));
        cvdata = cell2mat(cv_gpfa_cell(7,:));
        mse_fold = sum(sum((cvdata-[test_data.data]).^2));

        mse(idim, ifold) = mse_fold;
        ll_test(idim, ifold) = ll_te;
        ll_train(idim, ifold) = ll_tr;
        paramsGPFA{idim, ifold} = params;
        orth_traje_tr{idim, ifold} = gpfa_traj;
        orth_traje_te{idim, ifold} = traj;
        fprintf('Dimension %d, fold %d', idim, ifold)
    end        
end
% best dimension and best across folds
[a, foldmax_te] = max(sum(ll_test));
[a, imax_te] = max(ll_test(:,foldmax_te));
[a, foldmax_tr] = max(sum(ll_train));
[a, imax_tr] = max(ll_test(:,foldmax_tr));

results.ll_test = ll_test;
results.ll_train = ll_train;
results.mse = mse;
results.dim = imax_te;
results.GPFA = paramsGPFA;
results.traj_tr = orth_traje_tr;
results.traj_te = orth_traje_te;

% Setup figure to sumarize results
close all
til = sprintf('Spatial binning');
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
mean_ll_test  = mean(ll_test,2);
mean_ll_train = mean(ll_train,2);
offset_te = mean(mean_ll_test);
offset_tr = mean(mean_ll_train);
plot(dims,mean_ll_test,'linewidth',2, 'linestyle','-.', 'color','k'), hold on
plot(dims,mean_ll_train,'linewidth',2, 'linestyle','-.', 'color','b')
plot(dims,ll_test,'k')
plot(dims,ll_train,'b')
plot(imax_te, mean_ll_test(imax_te),'k*', 'markersize',10)
plot(imax_tr, mean_ll_train(imax_tr),'bs', 'markersize',10)
xlabel('Latent Dimension')
ylabel('Log Likelihood')
title('Train','color','b')
box off
namePNG = sprintf('Results/%s_cells',animal);
print(gcf,[basepath namePNG],'-dpng')

figure(iwin + 1)
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
namePNG = sprintf('Results/%s_LV_w%d',animal, window);
print(gcf,[basepath namePNG],'-dpng')    


save([basepath 'Results/' animal '_run_results.mat'],'results')
