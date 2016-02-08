%Figures for POSTER IBRO 2016
%
%
%Ruben Pinzon @ 2015

clc, close all; clear all;

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);
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


% ========================================================================%
%==============   (1) Show some place fields      ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType,...
                 wh_speed);

plot(D(1).X,D(1).Y,'color',[0.4 0.4 0.4]), hold on
plot(D(2).X,D(2).Y,'color',[0.4 0.4 0.4])

cell1  = D(1).spike_train(9,:)==1;
plot(D(1).X(cell1),D(1).Y(cell1),'r*')
cell2  = D(1).spike_train(10,:)==1;
plot(D(1).X(cell2),D(1).Y(cell2),'b*')
plot(D(1).X(cell1),D(1).Y(cell1),'r*')
cell2_2  = D(2).spike_train(10,:)==1;
plot(D(2).X(cell2_2),D(2).Y(cell2_2),'b*')
cell3_2  = D(2).spike_train(20,:)==1;
plot(D(2).X(cell3_2),D(2).Y(cell3_2),'b*')
cell4  = D(1).spike_train(30,:)==1;
plot(D(1).X(cell4),D(1).Y(cell4),'y*')

cell5  = D(1).spike_train(40,:)==1;
plot(D(1).X(cell5),D(1).Y(cell5),'g*')