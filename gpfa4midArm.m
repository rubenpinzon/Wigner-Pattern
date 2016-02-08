%BRANCH2_CLEANED This script contains a modularized version of the analysis
%        included in the script branch2d.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carried out most of the analysis in the files
%        branch2.m using functions. See branch2.m for further details.
%Version 1.0 Ruben Pinzon@2015


clc, close all; clear all;

basepath        = '/home/ruben/Documents/HAS/HC-5/';
[files, animals, roots]= get_matFiles(basepath);


%========================Paramteres and variables==========================
animal          = 2;
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
in              = 'mid_arm';
out             = 'mid_arm';
debug           = true;
namevar         = 'run';
%segmentation and filtering of silent neurons
bin_size        = 0.04; %ms
min_firing      = 1.0; %minimium firing rate
filterTrails    = false; % filter trails with irregular speed/spike count?
% GPFA trainign
n_folds         = 3;
zDim            = 5; %latent dimension
showpred        = false; %show predicted firing rate
train_split      = true; %train GPFA on left/right separately?
name_save_file  = '_trainedGPFA_midArm.mat';
test_lap        = 10;
maxTime         = 0; %maximum segmentation time 0 if use all
%%
% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType,...
                 wh_speed);

%show one lap for debug purposes 
if debug
    figure(test_lap)
    raster(D(test_lap).spikes), hold on
    plot(90.*D(test_lap).speed./max(D(test_lap).speed),'k')
    plot(90.*D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')
end

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

cgergo = load('colors');

colors = cgergo.cExpon([2 3 1], :);
labels = [R.type];
x_orth  = show_latent({M},R,colors, labels);

%======================================================================== %
%============== (6)    Save data                  ========================%
%=========================================================================%
fprintf('Will save at %s\n',[roots{animal} name_save_file])
save([roots{animal} name_save_file],'M','M_left','M_right','R', 'keep_neurons')
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

%%
%=========================================================================%
%=========(8) Compute loglike P(run|model_run)       =====================%
%=========================================================================%

load([roots{animal} name_save_file])                                       %Load trained model
R           = shufftime(R);                                                %shiffling time bins to destroy dynamics
R           = filter_laps(R);
%Classification stats of P(run events|model) 
models      = {M_left, M_right};
Xtats       = classGPFA(R, models);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

%show likelihood given the models
% plot show likelihood given the models
label.title = sprintf('P(run|run mid arm), zDim = %d',zDim);
label.modelA = 'Left alt.';
label.modelB = 'Right alt.';
label.xaxis = 'j';
label.yaxis = 'P(run_j| Models_{left run, right run})';
compareLogLike(R, Xtats, label)

%XY plot
cgergo = load('colors');
label.title = '';
label.xaxis = 'Log P(run|Model_{left run})';
label.yaxis = 'Log P(run|Model_{right run})';
LDAclass(Xtats, label, cgergo.cExpon([2 3], :))
%=========================================================================%
%=========(9) Compute loglike P(wheel|run_model)     =====================%
%=========================================================================%
%#TODO: Separate this part v in a different script

in              = 'wheel';
out             = 'wheel';
maxTime         = 6;
allTrials       = true; %use all trials of running to test since they are 
                        %all unseen to the wheel model

S = get_section(D, in, out, debug, namevar); %lap#1: sensor errors 
W = segment(S, bin_size, Fs, keep_neurons,...
                [namevar '_spike_train'], maxTime);
W = filter_laps(W);
W = W(randperm(length(W))); 

models      = {M_left, M_right};
Xtats       = classGPFA(W, models,[],allTrials);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'P(wheel_j | run model)';
label.modelA = 'Run rigth alt.';
label.modelB = 'Run left alt.';
label.xaxis = 'j';
label.yaxis = 'P(wheel_j|run model)';
compareLogLike(R, Xtats, label)

%XY plot
label.title = 'Class. with Fisher Disc.';
label.xaxis = 'P(wheel_j|run right)';
label.yaxis = 'P(wheel_j|run left)';
LDAclass(Xtats, label)     
%%
%=========================================================================% Requires loading the model and data struct and extract the latent
%=========(10) Classification of Xorth(start,mid,end) ====================% variables: steps (8, and 5) in that order
%=========================================================================% Requires also library prtools for the classifier

for k = 1 : length(x_orth)
 T           = size(x_orth{k},2);
 x_fea(k,:)  = [x_orth{k}(:,1)' x_orth{k}(:,ceil(T/2))' x_orth{k}(:,end)'];%Start, middle, and end points of trajectories as features (5 x 3) = 15 features     
end

label           = [R.type]';                                               %labels for the Bayes classifier (1, -1 => type 2)
label(label==3) = 1;


classrate       = bayes2c(x_fea,label,folds);


