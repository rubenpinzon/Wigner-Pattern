%% Analysis of Buzsaki database i01_maze05_MS.005
% L/R alternation task, 18 min, 15sec wheel runs, several errors


clear all, clc, close all;

basepath        = '/media/bigdata/i01_maze05.005/';
obj             = load([basepath 'i01_maze05_MS.005_BehavElectrData.mat']);
clusters        = obj.Spike.totclu;
laps            = obj.Laps.StartLaps(obj.Laps.StartLaps~=0); %@1250 Hz
%Adding the end of the last lap because Laps only contains the start
laps(end+1)     = obj.Par.SyncOff;
mazesect        = obj.Laps.MazeSection;
wheelspeed      = obj.Laps.WhlSpeedCW;
X               = obj.Track.X;
Y               = obj.Track.Y;
events          = obj.Par.MazeSectEnterLeft;
Fs              = obj.Par.SamplingFrequency;
speed           = obj.Track.speed;
eeg             = obj.Track.eeg;
numLaps         = numel(laps)-1;
time_seg        = linspace(0, length(eeg)/Fs, length(eeg));
[spk, spk_lap]  = get_spikes(clusters, obj.Spike.res, laps);
frate_mean      = 0.1*meanfrate(spk,0.1*Fs,obj.Par.SyncOff); %has to interpolate to keep #samples
save4attila(obj.Laps.TrialType,[basepath 'i01_maze05.005_trialType'],'')
%%
save4attila([obj.Spike.res/Fs*1000 clusters], [basepath 'spikes'], '%5.2f\t%3d\n')
save4attila(eeg,[basepath 'eeg'],'')
save4attila(mazesect,[basepath 'sections'],'')
save4attila(wheelspeed,[basepath 'wheelspeed'],'')
save4attila(speed,[basepath 'mousespeed'],'')
save4attila([obj.Clu.isIntern], [basepath 'isIntern'], '%d\n')