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
%==============             Extract lap           ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType);

% ========================================================================%
%==============       Find stable pcells          ========================%
%=========================================================================%

debug       = true;
in          = 'turn';
out         = 'lat_arm';
min_speed   = 0;
show_firing = false;
[E, S] = get_pfields(D, in,out, min_speed, 'easy', debug, show_firing);
%#TODO: modify the criteria to remove laps because it depends on the section of the maze

%Get consolitated indixes of stable place cells
stable_E = [];
for e = 1 : length(E)
    stable_E = [stable_E E{e}];    
end
stable_pcells = unique(stable_E); %these are the cells with stable pfield
name          = 'stable-pfields.mat';
save([roots{animal} name],'stable_pcells')


