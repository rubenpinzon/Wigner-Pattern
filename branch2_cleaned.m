% 

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

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
isIntern        = data.Clu.isIntern;
numLaps         = length(events);
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
n_cells         = size(spk_lap,2);
n_pyrs          = sum(isIntern==0);
TrialType       = data.Laps.TrialType;
Typetrial_tx    = {'left', 'right', 'errorLeft', 'errorRight'};

% ========================================================================%
%==============    Extract trials                 ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType);

% ========================================================================%
%==============    Extract Running Sections       ========================%
%=========================================================================%
in      = 'turn';
out     = 'lat_arm';
debug   = true;
namevar = 'run';
R       = get_section(D, in, out, debug, namevar);

% ========================================================================%
%==============    Segment the spike vectors      ========================%
%=========================================================================%
bin_size            = 0.04; %ms
min_firing          = 0.5; %minimium firing rate
[D,keep_neurons]    = segment(R, bin_size, Fs, min_firing, [namevar '_spike_train']);

% ========================================================================%
%==============             Train GPFA            ========================%
%=========================================================================%

zDim        = 10; 
showpred    = false;
M           = trainGPFA(D, zDim, showpred);


