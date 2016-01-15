%GPFA_CLEANED This script contains a modularized version of the analysis
%        included in the script branch2.m, that process the HC-5 database.
%
%        DESCRIPTION: load files from HC-5 database, segment in run and 
%        wheel sections, train a GPFA model with all the data using cross
%        validation, and asses the likelihood of the data.
%
%
%Version 1.0 Ruben Pinzon@2015

clc, close all; clear all;

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);

% ========================================================================%
% ========================     Paramters      ============================%
% ========================================================================%

animal          = 6;                                                       %7 animals in HC-5 database
data            = load(files{animal});                                     %read the processed file *BehavElectrData.mat 
clusters        = data.Spike.totclu;                                       %Vector of spike times
laps            = data.Laps.StartLaps(data.Laps.StartLaps~=0);             %Time of the start of each lap
laps(end+1)     = data.Par.SyncOff;                                        %Total duration of the experiment
mazesect        = data.Laps.MazeSection;                                   %Vector indicating the sections in the maze
events          = data.Par.MazeSectEnterLeft;                              %Time of entry/exit of the animal sections of the maze
Fs              = data.Par.SamplingFrequency;                              %@1250 
X               = data.Track.X;                                            %Animal X tracked position 
Y               = data.Track.Y;                                            %Animal Y tracked position
speed           = data.Track.speed;                                        %Animal speed 
wh_speed        = data.Laps.WhlSpeedCW;                                    %Wheel speed 
isIntern        = data.Clu.isIntern;                                       %Boolean indicating inhibitory cells 
numLaps         = length(events);                                          %Number of Laps/trials 
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);               %Function to extract spikes_per_cell and per_lap
n_cells         = size(spk_lap,2);                                         %Total number of cells
n_pyrs          = sum(isIntern==0);                                        %Number of pyramidal cells 
TrialType       = data.Laps.TrialType;                                     %Trial's type 1:left, 2:right, 3:lefterror, 4righterror 
Typetrial_tx    = {'left', 'right', 'errorLeft', 'errorRight'};            %Name tags for the trial's type 
debug           = true;                                                    %Enable verbose output and diagnostic plots
colors          = [1 0 0; 0 0 1; 0.1 0.1 0.1; 0.1 0.1 0.1];                %Colors for the four types of trials

clear data

% ======= Expecify the sections of interest. options:
%======== [mid_arm, preturn, turn, lat_arm, reward, delay, wheel]
in_wh           = 'wheel';                                                 %Extract from wheel
out_wh          = 'wheel';                                                 %To the end of the wheel section 
namevar_wh      = 'wheel';                                                 %name for the variable where extracted data is stored 
in_run          = 'mid_arm';                                               %for the run section start at mid arm
out_run         = 'lat_arm';                                               %until the end of the lateral arms 
namevar_run     = 'run';                                                   %name of the variable "run"
maxTime_wh      = 6;                                                       %maximum segmentation time seconds
maxTime_run     = 0;                                                       % 0: use all the recording

%============segmentation and filtering of silent neurons
bin_size        = 0.04;                                                    %Segmenting window in seconds
min_firing      = 1.0;                                                     %minimium firing rate to filter cells

%=========== GPFA training
n_folds         = 3;                                                       %Folds for crossvalidation
zDim            = 10;                                                      %latent dimension
showpred        = false;                                                   %show predicted firing rate
name_save_file  = '_trainedGPFA_all.mat';                                  %Name of the file where trained model will be stored                                    
test_lap        = 10;                                                      %Diagnostic lap
filterlaps      = false;                                                   %filter laps by average speed?
%%
% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType,...
                 wh_speed);

%show raster
if debug
    figure(test_lap)
    raster(D(test_lap).spikes), hold on
    plot(90.*D(test_lap).speed./max(D(test_lap).speed),'k')
    plot(90.*D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')
end

% ========================================================================%
%==============  (2)  Extract Sections            ========================%
%=========================================================================%

R = get_section(D, in_run, out_run, debug, namevar_run); 
W = get_section(D, in_wh, out_wh, debug, namevar_wh); 

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%

[R,keep_neurons]    = segment(R, bin_size, Fs, min_firing,...
                              [namevar_run '_spike_train'], maxTime_run);  %segment run trial 
W                   = segment(W, bin_size, Fs, keep_neurons,...
                              [namevar_wh '_spike_train'], maxTime_wh);    %segment wheel trial with the same neurons *kept* in "run" 
                          
if filterlaps
    R               = filter_laps(R);
    W               = filter_laps(W);
end

% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
M                 = trainGPFA([W R], zDim, showpred, n_folds);             %Use all the sections [W R]
%#TODO: The trialId gets confused in the joint of [W R], add to M field 
%eventNo. 
%%
% ========================================================================%
%============== (5)    Show Neural Trajectories   ========================%
%=========================================================================%
labels = [ones(1,length(W)) 2*ones(1,length(R))];                          %Here the labels can be used to show trial type, sections
% labels = [W.type];
colors = cgergo.cExpon;
Xorth  = show_latent({M},[R W], colors, labels);

%======================================================================== %
%============== (6)    Save data                  ========================%
%=========================================================================%
fprintf('Will save at %s\n',[roots{animal} name_save_file])
save([roots{animal} name_save_file],'M','W','R', 'keep_neurons')


models      = {M};
Xtats       = classGPFA(W, models,[],'all');

%======================================================================== %
%============== (7)    Load the saved model       ========================%
%======================================================================== %

load([roots{animal} name_save_file])                                       %Models already trained

%======================================================================== %
%======(8)    Classification based on latent variables    ================%
%======================================================================== %

traj = exactInferenceWithLL(W, M.params{1},'getLL',0);                     %Choose models from fold #1
W    = segmentByTrial(W,orthogonalize([traj.xsm], M.params{1}.C),'Xorth'); %ortogonalize and segment by trial latent vars    
P    = extractPoints(W);                                                   %Extract points of interest from the latent variables (stat, 1/3, 1/2, 2/3, end)


