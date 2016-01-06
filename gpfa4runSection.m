%BRANCH2_CLEANED This script contains a modularized version of the analysis
%        included in the script branch2d.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carried out most of the analysis in the files
%        branch2.m using functions. See branch2.m for further details.
%Version 1.0 Ruben Pinzon@2015


clc, close all; clear all;

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);


%========================Paramteres and variables==========================
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
in              = 'mid_arm';
out             = 'lat_arm';
debug           = true;
namevar         = 'wheel';
%segmentation and filtering of silent neurons
bin_size        = 0.04; %ms
min_firing      = 0.8; %minimium firing rate
filterTrails    = false; % filter trails with irregular speed/spike count?
% GPFA trainign
n_folds         = 2;
zDim            = 10; %latent dimension
showpred        = false; %show predicted firing rate
train_split      = true; %train GPFA on left/right separately?
name_save_file  = '_trainedGPFA_run.mat';
test_lap        = 10;
maxTime         = 0; %maximum segmentation time 0 if use all

% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType,...
                 wh_speed);

%show one lap for debug purposes 
figure(test_lap)
raster(D(test_lap).spikes), hold on
plot(90.*D(test_lap).speed./max(D(test_lap).speed),'k')
plot(90.*D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')

% ========================================================================%
%==============  (2)  Extract Running Sections    ========================%
%=========================================================================%

S = get_section(D, in, out, debug, namevar); %lap#1: sensor errors 

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%
%load run model and keep the same neurons
% run = load([roots{animal} '_branch2_results40ms.mat']);

[R,keep_neurons]    = segment(S, bin_size, Fs, min_firing,...
                              [namevar '_spike_train'], maxTime);
%%
% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
M                 = trainGPFA(R, zDim, showpred, n_folds);

if train_split
    [R_left, R_right] = split_trails(R);
    if filterTrails
        R_left            = filter_laps(R_left);
        R_right           = filter_laps(R_right,'bins');
    end

    M_left            = trainGPFA(R_left, zDim, showpred, n_folds);
    M_right           = trainGPFA(R_right, zDim, showpred, n_folds);
end

%%
% ========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%

Xorth = show_latent(M_right, R_right);

%======================================================================== %
%============== (6)    Save data                  ========================%
%=========================================================================%

save([roots{animal} name_save_file],'M','M_left','M_right','W')
%%
%=========================================================================%
%=========(7) Compare mean spike counts              =====================%
%=========================================================================%
figure(7)
set(gcf,'position',[100 100 500*1.62 500],'color','w')
plot(mean([R_left.y],2),'r','displayname','wheel after left')
hold on
plot(mean([R_right.y],2),'b','displayname','wheel after right')
ylabel('Average firing rate')
xlabel('Cell No.')
set(gca,'fontsize',14)
savefig()
%=========================================================================%
%=========(8) Compute loglike P(wheel|model_wheel)   =====================%
%=========================================================================%

load([roots{animal} name_save_file])

%Classification stats of P(proto_event|model) 
models      = {M_left, M_right};
Xtats       = classGPFA(W, models);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

%show likelihood given the models
plot(Xtats.likelihood','-o')
xlabel('Trials'), xlim([0 length(W)+1]), 
ylabel('LogLikelihood')
for t = 1 : length(W)
  typeAssigned = '2';
  
  if Xtats.likelihood(1,t) > Xtats.likelihood(2,t)
     typeAssigned = '1'; 
  end
  c = 'r';
  if Xtats.class_output(t) == Xtats.real_label(t)
      c = 'k';
  end
  text(t, -8700, typeAssigned, 'color',c)
  line([t+0.5 t+0.5],ylim,'linestyle','--','color',[0.6 0.6 0.6])
end
set(gca,'xticklabel',[W.trialId])

     