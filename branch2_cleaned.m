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
animal          = 4;
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
namevar         = 'wrun';
%segmentation and filtering of silent neurons
bin_size        = 0.04; %ms
min_firing      = 1.0; %minimium firing rate
% GPFA trainign
n_folds         = 2;
zDim            = 10; %latent dimension
showpred        = false; %show predicted firing rate
train_arms      = true; %train GPFA on arms separately
name_save_file  = '_trainedGPFA.mat';
% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType, wh_speed);

% ========================================================================%
%==============  (2)  Extract Running Sections    ========================%
%=========================================================================%

R = get_section(D, in, out, debug, namevar);

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%

[D,keep_neurons]    = segment(R, bin_size, Fs, min_firing, [namevar '_spike_train']);
D                   = filter_laps(D);

%%
% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
M                 = trainGPFA(D, zDim, showpred, n_folds);
model{1}          = M.params{1}; %one of the three folds
data{1}           = D;  

if train_arms
    [D_left, D_right] = split_trails(D);
    M_left            = trainGPFA(D_left, zDim, showpred, n_folds);
    M_right           = trainGPFA(D_right, zDim, showpred, n_folds);
    model{2}          = M_left.params{1};
    model{3}          = M_right.params{1}; 
    data{2}           = D_left;
    data{3}           = D_right; 
end
%#TODO there is a problem with the crossvalidation when one of the folds 
%includes trials with too few spikes.
%%
% ========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%

Xorth = show_latent(model, data);

%======================================================================== %
%============== (6)    Save data                  ========================%
%=========================================================================%

save([roots{animal} name_save_file],'model','data','-v7.3')
