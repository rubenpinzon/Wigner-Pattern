%Proto Events This script contains a modularized version of the analysis
%        included in the script branch2d.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carries out P(proto events |run models)
%Version 1.0 Ruben Pinzon@2015

clc, close all; clear all;

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);


%========================Paramteres and variables==========================
animal          = 6;
Typetrial_tx    = {'left', 'right', 'errorLeft', 'errorRight'};
%segmentation and filtering of silent neurons
bin_size        = 0.04; %ms
min_firing      = 1.0; %minimium firing rate
filterTrails    = false; % filter trails with irregular speed/spike count?
% GPFA trainign
n_folds         = 3;
zDim            = 10; %latent dimension
showpred        = false; %show predicted firing rate
train_split      = true; %train GPFA on left/right separately?
test_lap        = 10;
maxTime         = 0; %maximum segmentation time 0 if use all
allTrials       = true; %use all the proto envetns in testing since they 
                        %unseen by the models
%=========================================================================%
%============== (1)  Load Proto Events from file branch2d  ===============%
%=========================================================================%

load([roots{animal} '_proto_events_v2.mat'])

%=========================================================================%
%=========               (2) Load Run models               ===============%
%=========================================================================%

load([roots{animal} '_trainedGPFA_run.mat'])

%=========================================================================%
%=========            (3) Segment Proto Events             ===============%
%=========================================================================%

PR = segment(P, bin_size/10, Fs, keep_neurons,'data', maxTime);
PR = filter_laps(PR);
%=========================================================================%
%=========            (4) P(proto event |run model)        ===============%
%=========================================================================%

models      = {M_left, M_right};
Xtats       = classGPFA(PR, models,[],allTrials);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'P(proto_j | run model)';
label.modelA = 'Run rigth alt.';
label.modelB = 'Run left alt.';
label.title = 'Class. with Fisher Disc.';
label.xaxis = 'P(proto_j|run right)';
label.yaxis = 'P(proto_j|run left)';
LDAclass(Xtats, label)     


%=========================================================================%
%=======(5) time shuffled P(proto event |run model)        ===============%
%=========================================================================%
PS           = shufftime(P);
models      = {M_left, M_right};
Xtats       = classGPFA(PS, models,[],allTrials);
cm          = [Xtats.conf_matrix];
fprintf('Min-max Classifier hit1: %2.2f%%, hit2: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'P(proto_j | run model)';
label.modelA = 'Run rigth alt.';
label.modelB = 'Run left alt.';
label.title = 'Class. with Fisher Disc.';
label.xaxis = 'P(proto_j (shuf. t)|run right)';
label.yaxis = 'P(proto_j (shuf. t)|run left)';
LDAclass(Xtats, label) 

%=========================================================================%
%=======(6) Reverse time P(proto event |run model)         ===============%
%=========================================================================%

Prev         = reversetime(PR);

models      = {M_left, M_right};
Xtats       = classGPFA(Prev, models,[],allTrials);
cm          = [Xtats.conf_matrix];
fprintf('Min-max Classifier hit1: %2.2f%%, hit2: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'reverse bins P(proto_j | run model)';
label.modelA = 'Run rigth alt.';
label.modelB = 'Run left alt.';
label.title = 'Class. with Fisher Disc.';
label.xaxis = 'P(proto_j (rev. t)|run right)';
label.yaxis = 'P(proto_j (rev. t)|run left)';
LDAclass(Xtats, label) 

