%GPFA4WHEELSECTION This script contains a modularized version of the analysis
%        included in the script branch2.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carried out most of the analysis in the files
%        branch2.m using functions. See branch2.m for further details.
%Version 1.0 Ruben Pinzon@2015


clc, close all; clear all;

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
debug           = true;
namevar         = 'wheel';
%segmentation and filtering of silent neurons
bin_size        = 0.04; %ms
min_firing      = 1.0; %minimium firing rate
% GPFA trainign
n_folds         = 3;
zDim            = 10; %latent dimension
showpred        = false; %show predicted firing rate
train_split      = true; %train GPFA on left/right separately?
name_save_file  = '_trainedGPFA_wheel.mat';
test_lap        = 10;
maxTime         = 6; %maximum segmentation time
filterlaps      = false;
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

R = get_section(D(2:end), in, out, debug, namevar); %lap#1: sensor errors 

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%
%load run model and keep the same neurons
% run = load([roots{animal} '_branch2_results40ms.mat']);

[W,keep_neurons]    = segment(R, bin_size, Fs, min_firing,...
                              [namevar '_spike_train'], maxTime);
if filterlaps
    W                   = filter_laps(W);
end

%%
% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
M                 = trainGPFA([W R], zDim, showpred, n_folds);

if train_split
    [W_left, W_right] = split_trails(W);
    M_left            = trainGPFA(W_left, zDim, showpred, n_folds);
    M_right           = trainGPFA(W_right, zDim, showpred, n_folds);
end

%%
% ========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%
cgergo = load('colors');

colors = cgergo.cExpon([2 3 1], :);
labels = [W.type];
figure
Xorth = show_latent({M},W, colors, labels);

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
