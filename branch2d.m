%Branch 2d IDENTIFYING REPLAY EVENTS AFTER ACTIVITY
%
%
%Here the stopping periods after running are analyzed following the
%procedure described in Foster&Wilson 2006, in which A spike train was 
%constituted from all spikes (from all cellsin the probe sequence) that 
%occurred during stopping periods while the animal faced in the direction
%in which it had just run. This spike train was then broken between every
%pair of successive spikes separated by more than 50 ms, to form a large 
%set of proto-events. Those proto-events in which at least one-third of 
%the cells in the probe sequence fired at least one spike were then 
%selected as events. The few events longer than 500 ms in duration were 
%rejected as a potential source of spurious correlations. 

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);


%========================Variables of Interest===========================
animal          = 5;
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
numLaps         = numel(laps)-1;
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
n_cells         = size(spk_lap,2);
color           = jet(55);
%%
% ========================================================================%
%==============   Extract Stopping section after run =====================%
%=========================================================================%

debug           = true; %to show diganostic plots
speed_th        = 200;
%this is to remove/add the section in the middle arm of the maze
sect_in         = [7, 8]; 
sect_out        = [7, 8];
cnt             = 1;
% Extract spks when the mouse is running and in the wheel to calculate
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    
    idx_run                 = [sum(events{lap}(sect_in,1)), sum(events{lap}(sect_out,2))];
    int_at_maze(lap, :)     = idx_run;
    run_len(lap)            = (idx_run(2)-idx_run(1))/Fs;
    X_lap{lap}              = X(idx_run(1):idx_run(2));
    Y_lap{lap}              = Y(idx_run(1):idx_run(2));
    speed_lap               = speed(idx_run(1):idx_run(2));
    
    %speed below threshold
    speed_lap(speed_lap<speed_th) = 1;
    speed_lap(speed_lap>=speed_th) = 0;
    %extract regions in which the animal is still
    dist        = diff(speed_lap);
    moved       = -find(dist==1);
    stoped      = find(dist==-1);
    period      = [stoped  moved(1:length(stoped))];
    %select those stoppig periods larger than 1s
    winners     = find(sum(period,2) > 1.0*Fs);
    
    for w = 1:length(winners)
        idx_stop    = [-period(winners(w),2) period(winners(w),1)]; 
        
        %spikes
        no_spks   = zeros(1,n_cells);
        for neu=1:n_cells
            times = spk_lap{lap,neu}>=idx_stop(1) & spk_lap{lap,neu}<=idx_stop(end);
            %aligned to the start of the section
            no_spks(neu) = sum(times);
            if no_spks(neu) ~= 0
                data{neu} = spk_lap{lap,neu}(times) - idx_run(1) + 1;
            end
        end
        S(cnt).spk_times   = data;
        S(cnt).NoSpikes    = no_spks;
        S(cnt).TrialId     = lap;
        S(cnt).TrialType   = typetrial{data.Laps.TrialType(laps(lap))};
        S(cnt).TrialNo     = data.Laps.TrialType(laps(lap));
        S(cnt).Interval    = idx_stop;
        S(cnt).Duration    = sum(idx_stop);
        cnt                = cnt + 1;
    end
    
end    

