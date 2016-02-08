%GPFA4WHEELWINDOW This script contains the same procedure in
%                 gpfawheelSection.m for the wheel section of the HC-5
%                 database, however, the activity in segmented every 3 sec,
%                 which is of similar scale as the running section. SLiding
%                 windows with different overlapping are considered. With
%                 this data two GPFA estimators conditioned on wheel events
%                 after left and right alternations is trained. Then the
%                 likelihood of each event which each model is compared to
%                 test the prediction power of the estimator.
%
%See also gpfawheelSection.m, crcns.org/data-sets/hc/hc-5/about-hc-5
%Version 1.0 Ruben Pinzon@2015


clc, close all; clear all;

basepath        = '/home/ruben/Documents/HAS/HC-5/';
[files, animals, roots]= get_matFiles(basepath);


%========================Variables of Interest===========================
animal          = 1;
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
filterlaps      = false;
%Segmenting the wheel events
maxTime         = 3; %maximum segmentation time
%colors
cgergo          = load('colors');

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
    raster(D(test_lap).spikes),hold on 
    plot(90*D(test_lap).speed./max(D(test_lap).speed),'k')
    plot(90*D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')
    xlabel('Samples')
end

% ========================================================================%
%==============  (2)  Extract Wheel Sections      ========================%
%======  Here the wheel section has the full lenght ~ 10s ================%
%=========================================================================%

W = get_section(D(1:end-1), in, out, debug, namevar); %lap#1: sensor errors 

% ========================================================================%
%==== (3)  Segment wheel events with different windows        ============%
%=========================================================================%

[Wh,keep_neurons]    = segment(W, bin_size, Fs, min_firing,...
                              [namevar '_spike_train'], maxTime);
%%
% ========================================================================%
%========= (4)    Train GPFA for each window        ======================%
%=========================================================================%
for vars = 1 : length(Wh) 
    fprintp('Data F%d, training joint estimator\n');
    M                 = trainGPFA(Wh{vars}, zDim, showpred, n_folds);

    if train_split
        fprintp('Data F%d, training joint estimator\n');
        [W_left, W_right] = split_trails(Wh{vars});
        M_left            = trainGPFA(W_left, zDim, showpred, n_folds);
        M_right           = trainGPFA(W_right, zDim, showpred, n_folds);
    end

%=========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%

    colors = cgergo.cExpon([2 3 1], :);
    labels = [Wh{vars}.type];
    figure(vars)
    Xorth = show_latent({M},Wh{vars}, colors, labels);                          
    Model{vars} = {M, M_left, M_right}; 
    
end
