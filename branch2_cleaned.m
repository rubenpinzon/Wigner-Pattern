%BRANCH2_CLEANED This script contains a modularized version of the analysis
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
maxTime         = 5; %maximum segmentation time
% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType,...
                 wh_speed);

%show raster
figure(test_lap)
raster(D(test_lap).spikes), hold on
plot(90.*D(test_lap).speed./max(D(test_lap).speed),'k')
plot(90.*D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')

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
W                   = filterlaps(W);
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
%#TODO there is a problem with the crossvalidation when one of the folds 
%includes trials with too few spikes.
%%
% ========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%

Xorth = show_latent(M_right, W_right);

%======================================================================== %
%============== (6)    Save data                  ========================%
%=========================================================================%

save([roots{animal} name_save_file],'M','M_left','M_right','W')
%%
%=========================================================================%
%=========(7) Compute loglike P(wheel|model_wheel)   =====================%
%=========================================================================%

%#TODO: during the classification all the P-trails should be evaluated with
%models that did not include the same trial in the training. W and M is not
%a problem because the trials coincide with the cv_trial vector. W_left and
%W_right on the contrary losses the trialId during the training so that 
%cv_trial is not informative.

load([roots{animal} name_save_file])

%Classification stats of P(proto_event|model) 
models      = {M_left, M_right};
Xtats       = classGPFA(W, debug, models);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*mean(cm(1,1)),...
        100*mean(cm(2,2)))


plot(Xtats.likelihood','o')
xlabel('Trials'), xlim([0 length(W)+1]), ylim([-5600 -4000])
ylabel('LogLikelihood')
for t = 1 : length(W)
  typeAssigned = '2';
  if stats(f).likelihood(1,t) > stats(f).likelihood(2,t)
     typeAssigned = '1'; 
  end
  text(t, -5500, typeAssigned)
end

title('P(wheel(T=3s) | wheel models)')

     