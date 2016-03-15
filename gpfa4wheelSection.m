%GPFA4WHEELSECTION This script contains a modularized version of the analysis
%        included in the script branch2.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carried out most of the analysis in the files
%        branch2.m using functions. See branch2.m for further details.
%Version 1.0 Ruben Pinzon@2015

clc, close all; clear all;

basepath        = '/home/ruben/Documents/HAS/HC-5/';
[files, animals, roots]= get_matFiles(basepath);

%========================Variables of Interest===========================
animal          = 6;
fprintf('Loading animal %s\n',animals{animal});
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
time            = linspace(0, length(eeg)/1250,length(eeg));
speed           = data.Track.speed;
wh_speed        = data.Laps.WhlSpeedCW;
isIntern        = data.Clu.isIntern;
numLaps         = length(events);
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
n_cells         = size(spk_lap,2);
n_pyrs          = sum(isIntern==0);
TrialType       = data.Laps.TrialType;
Typetrial_tx    = {'left', 'right', 'errorLeft', 'errorRight'};
clear data
%section in the maze to analyze
in              = 'wheel';
out             = 'wheel';
debug           = false;
namevar         = 'wheel';
%segmentation and filtering of silent neurons
bin_size        = 0.04; %ms
min_firing      = 1.0; %minimium firing rate
% GPFA trainign
n_folds         = 3;
zDim            = 10; %latent dimension
showpred        = false; %show predicted firing rate
train_split      = true; %train GPFA on left/right separately?
name_save_file  = '_trainedGPFA_wheel_bins25to100.mat';
test_lap        = 10;
maxTime         = 6; %maximum segmentation time
filterlaps      = false;
cgergo          = load('colors');
colors          = cgergo.cExpon([2 3 1 1], :);
intervalBins    = [25 105];
% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%
%%
D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType,...
                 wh_speed);

%show raster
if debug
    figure(test_lap)
    set(gcf,'position',[100 100 500*1.62 500])
    time = linspace(0,length(D(test_lap).speed)/Fs,length(D(test_lap).speed));
    raster(D(test_lap).spikes), 
    plot(time,D(test_lap).speed./max(D(test_lap).speed),'k'),hold on
    plot(time,D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')
    xlabel('Time (s)')
end

% ========================================================================%
%==============  (2)  Extract Running Sections    ========================%
%=========================================================================%

R = get_section(D(2:end-1), in, out, debug, namevar); %lap#1: sensor errors 

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%
%load run model and keep the same neurons
% run = load([roots{animal} '_branch2_results40ms.mat']);

[W,keep_neurons]    = segment(R, bin_size, Fs, min_firing,...
                              [namevar '_spike_train'], maxTime);
                          
W                   = filterbins(W,intervalBins(1),intervalBins(2));       % get the last part of the wheel section bins (100-225)    

if filterlaps
    W               = filter_laps(W);
end

%%
% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
M                 = trainGPFA(W, zDim, showpred, n_folds);

if train_split
    [W_left, W_right] = split_trails(W);
    M_left            = trainGPFA(W_left, zDim, showpred, n_folds);
    M_right           = trainGPFA(W_right, zDim, showpred, n_folds);
end

%%
% ========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%
load([roots{animal} name_save_file])

labels = [W.type];
x_orth = show_latent({M},W, colors, labels, debug);                       %Lantent with joint model
%%
%======================================================================== %
%============== (6)    Save data                  ========================%
%=========================================================================%
fprintf('Will save at %s\n',[roots{animal} name_save_file])
save([roots{animal} name_save_file],'M','M_left','M_right','W','R', 'keep_neurons')
%%
%=========================================================================%
%=========(7) Compare mean spike counts              =====================%
%=========================================================================%
figure(7)
set(gcf,'position',[100 100 500*1.62 500],'color','w')
plot(mean([W_left.y],2),'r','displayname','wheel after left')
hold on
plot(mean([W_right.y],2),'b','displayname','wheel after right')
ylabel('Average firing rate')
xlabel('Cell No.')
set(gca,'fontsize',14)
%%
%=========================================================================%
%=========(8) Compute loglike P(wheel|model_wheel)   =====================%
%=========================================================================%

%If model was not trained it can be loaded:
load([roots{animal} name_save_file])

%transformation to W testing
%W           = W(randperm(length(W))); %permutation of laps
%W           = shufftime(W); %time shuffling for each lap

errorTrials = find([W.type] > 2);                                          %erroneous trials wheel events
We          = W(errorTrials);                                              %erroneous trials struct                 

%Classification stats of P(proto_event|model) 
models      = {M_right, M_left};                                           %here models have future run label, 
Xtats       = classGPFA(W, models);
cm          = [Xtats.conf_matrix];
fprintf('Max-min Classifier hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'P(wheel_j after error | models W)';
label.modelA = 'Wheel after rigth alt.';
label.modelB = 'Wheel after left alt.';
label.xaxis = 'j';
label.yaxis = 'P(wheel_j|model)';
compareLogLike(W, Xtats, label)                                           %P(error W | models W)

%XY plot
label.title = '';
label.modelA = 'Wheel after left alt.';
label.modelB = 'Wheel after right alt.';
label.xaxis = 'Log P(wheel|Model_{wheel after left run})';
label.yaxis = 'Log P(wheel|Model_{wheel after right run})';
LDAclass(Xtats, label, cgergo.cExpon([2 3], :))


%=========================================================================%
%=========(9) Compute loglike P(run|model_wheel)     =====================%
%=========================================================================%
%#TODO: Separate this part v in a different script

in              = 'preturn';
out             = 'lat_arm';
maxTime         = 0;
allTrials       = true; %use all trials of running to test since they are 
                        %all unseen to the wheel model

S = get_section(D, in, out, debug, namevar); %lap#1: sensor errors 
R = segment(S, bin_size, Fs, keep_neurons,...
                [namevar '_spike_train'], maxTime);
R = filter_laps(R);
%R = R(randperm(length(R))); 
models      = {M_left, M_right};
Xtats       = classGPFA(R, models,[],allTrials);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'P(run_j | wheel model)';
label.modelA = 'Run rigth alt.';
label.modelB = 'Run left alt.';
label.xaxis = 'j';
label.yaxis = 'P(run_j|wheel model)';
compareLogLike(R, Xtats, label)

%XY plot
label.title = 'Class. with Fisher Disc.';
label.xaxis = 'P(run_j|wheel model right)';
label.yaxis = 'P(run_j|wheel model left)';
LDAclass(Xtats, label)

%%
%=========================================================================% Requires loading the model and data struct and extract the latent
%=========(11) Classification of Xorth point by point ====================% variables: steps (8, and 5) in that order
%=========================================================================% Requires also library prtools for the classifier
                                                                          % An issue is the nonuniform lenght of x_orths
%W=filter_laps(W);                                                        % Run interpolacion to make them uniform

load([roots{animal} name_save_file])                                      % If model was not trained it can be loaded:
label           = [W.type]';                                              %
x_orth          = show_latent({M},W, colors, label, debug);               %Lantent with joint model
%% save latents to mat file and classification
label(label==3) = 1;
label(label==4) = 2;
len_x           = min([W.T]);                                             % min len to cut all trajetories to the same length
lap_int         = 1:52;
n_lat           = 1:10;
for t = 1 : len_x 
    x_fea          = zeros(len(lap_int),len(n_lat));
    for k = 1 : length(lap_int)                                            % Each trial    
        kk         = lap_int(k);        
        x_fea(kk,:) = [x_orth{kk}(n_lat,t)];    
    end
    [rate, folds]   = bayes2c(x_fea,label(lap_int),n_folds,[W.trialId]);
    Fold(t).infoClass = folds;
    classrate(t,:)  = rate; 
end
    
%Data flatted for knn
x_knn = zeros(10,len(x_orth{1}),len(x_orth));
for lap = 1 : len(x_orth)
   x_knn(:,:,lap) = x_orth{lap};    
end
save(sprintf('%sx_orth%s.mat',roots{animal},name_save_file),'x_knn','label')
%%
figure(), hold on                                                         % Show accuracy
%subplot(211)
set(gcf,'color','w')
errorbar(classrate(:,1),classrate(:,2))                                   % Total accuracy
set(gca,'fontsize',14,'fontname','georgia')
grid on
xlabel('Bins'), ylabel('Total Accuracy')
xlim([1, 150])
line([1, 150],[mean(classrate(:,1)) mean(classrate(:,1))],'color','k','linewidth',2)
title(sprintf('Number of laps %d: mean=%3.3f',len(W),mean(classrate(:,1))))
%%
subplot(212)
meanWh_speed_left   = zeros(150*bin_size*Fs,sum([R.type]==1));
meanWh_speed_right  = zeros(150*bin_size*Fs,sum([R.type]==2));
meanWh_speed_right  = zeros(150*bin_size*Fs,len(R));

left                = 1;
right               = 1;
for l = 1 : length(R)
    spe     = R(l).wheel_angularVelo(1:bin_size*Fs:150*bin_size*Fs);      % Wheel speed (donwsampled with the bin size)
    plot(spe,'color',cgergo.cPhase(R(l).type,:),'Displayname',...
        sprintf('W lap %d',l),'marker','*'), hold on
    if R(l).type == 1
        meanWh_speed_left(:, left)   = R(l).wheel_angularVelo(1:150*bin_size*Fs);
        left = left + 1; 
    elseif R(l).type == 2
        meanWh_speed_right(:, right) = R(l).wheel_angularVelo(1:150*bin_size*Fs);
        right = right + 1;
    end
    meanWh(:, l) = R(l).wheel_angularVelo(1:150*bin_size*Fs); 
        
end
xlim([1, 150]), grid on
xlabel('Bins'), ylabel('Speed (mm/s)')
%% 
%=========================================================================% Requires loading the model and data struct and extract the latent
%======= (12) Wheel angular velocity profile colored as trial type =======% Steps (11-8-5)
%=========================================================================%

mu_left     = mean(meanWh_speed_left,2);
sd_left     = std(meanWh_speed_left,1,2);
mu_right    = mean(meanWh_speed_right,2);
sd_right    = std(meanWh_speed_right,1,2);
t       = linspace(0,maxTime,len(mu_left));

check_bin = 69;
output    = [];
trials    = [];
for j = 1 : n_folds
    output  = [output, Fold(check_bin).infoClass(j).ypred];
    trials  = [trials, Fold(check_bin).infoClass(j).trialId];
end

left                = 1;
right               = 1;
for l = 1 : len(R)
   if output(trials==l) == 1 
        mu_left_class(:,left) = meanWh(:, l);
        left = left + 1;
   elseif output(trials==l) == 2
        mu_right_class(:,right) = meanWh(:, l);
        right = right + 1;        
   end
    
end

mu_left_c     = mean(mu_left_class,2);
sd_left_c     = std(mu_left_class,1,2);
mu_right_c    = mean(mu_right_class,2);
sd_right_c    = std(mu_right_class,1,2);

figure, hold on
set(gcf,'color','w')
shadedErrorBar(t,mu_left',sd_left','g')
shadedErrorBar(t,mu_right',sd_right','r')
xlabel('Time(s)'), ylabel('Speed (mm/s)')
grid on

figure, hold on
set(gcf,'color','w')
shadedErrorBar(t,mu_left_c',sd_left_c','g')
shadedErrorBar(t,mu_right_c',sd_right_c','r')
xlabel('Time(s)'), ylabel('Speed (mm/s)')
title(['Color labeled according to classifier, bin#' num2str(check_bin)])
grid on

%%
%=========================================================================%
%=========(13) Studying the speed profile in the wheel====================%
%=========================================================================%
figure
t_max  = nan(length(R),3);
set(gcf,'color','w')
for l = 1 : length(R)
   subplot(3,1,R(l).type),hold on, grid on
   plot(R(l).wheel_06_002_angularVelo, 'color', cgergo.cExpon(R(l).type,:))  
   [a, b]  =  max(R(l).wheel_06_002_angularVelo);
   
   t_max(l, R(l).type) = b;
end

figure
subplot(211)
hist(t_max(:,1))
subplot(212)
hist(t_max(:,2))

%%
%=========================================================================%
%=========(14) Ellipses of the trajectories x_orth    ====================%
%=========================================================================%

n_latents    = [1 2];
lat_x_left   = zeros(length(n_latents)*sum([R.type]==1),150);
lat_x_right  = zeros(length(n_latents)*sum([R.type]==2),150);
left         = 1;
right        = 1;

for l = 1 : length(R)
    if R(l).type == 1
        lat_x_left(1 + (left-1)*(numel(n_latents)):left*(numel(n_latents)),:)   = x_orth{l}(n_latents,:);
        left = left + 1; 
    elseif R(l).type == 2
        lat_x_right(1 + (right-1)*(numel(n_latents)):right*(numel(n_latents)),:) = x_orth{l}(n_latents,:);
        right = right + 1;
    end    
end

ellipseLatent(lat_x_left, lat_x_right)

