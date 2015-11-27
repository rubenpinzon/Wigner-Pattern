%Simplified version of script BRANCH2.m
%
%
%Ruben Pinzon

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

%% Extract the laps

L = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType);

%% extract the section of interest
name        ='arms_only';
in          = 'preturn';
out         = 'lat_arm';
debug       = true;
S           = get_section(L, in, out, debug, name);

%=========================================================================%
%============       Train GPFA on the theta          =====================%
%=========================================================================%

bin_size        = 0.04;  %20 ms
zDim            = 10;    % Target latent dimensions
min_firing      = 1;
name_field      = 'arms_only_spike_train';
showpred        = true;
mod_tags        = {'_left', '_right', '_both'};
[D,keep_cell]   = segment(S, bin_size, Fs, min_firing, name_field);
[D_left, D_right] = split_trails(D);


for m = 1 : length(mod_tags)
    Data = eval(['D' mod_tags{m}]);
    M    =  trainGPFA(Data, zDim, showpred);
    R(m).model_tag = mod_tags{m};
    R(m).model     = M;
    clear M Data
end

%% Save data


