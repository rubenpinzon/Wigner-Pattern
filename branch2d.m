%Branch 2d IDENTIFYING REPLAY EVENTS AFTER ACTIVITY
%
%
%Here the stopping periods after running are analyzed following the
%procedure described in Foster&Wilson 2006 Neuron 36, A spike train was 
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
n_laps         = numel(laps)-1;
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
n_cells         = size(spk_lap,2);
color           = jet(55);
removeInh       = true;
%%
% ========================================================================%
%==============   Extract Stopping section after run =====================%
%=========================================================================%

debug           = false; %to show diganostic plots
speed_th        = 200;
%this is to remove/add the section in the middle arm of the maze
sect            = [3 4]; %without middle arm 
sect_in         = [7, 8]; 
sect_out        = [7, 8];
cnt             = 1;
% Extract spks when the mouse is running and in the wheel to calculate
for lap = 1:n_laps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    idx_run                 = [sum(events{lap}(sect,1)), sum(events{lap}(5:6,2))];
    if idx_run(1)>idx_run(2)
        fprintf('WARNING: Lap %d animal gor crazy, skipping\n',lap)
        continue
    end
    idx_stop                = [sum(events{lap}(sect_in,1)), sum(events{lap}(sect_out,2))];
    X_lap{lap}              = X(idx_stop(1):idx_stop(2));
    Y_lap{lap}              = Y(idx_stop(1):idx_stop(2));
    speed_lap               = speed(idx_stop(1):idx_stop(2));
    
    %speed below threshold
    speed_lap(speed_lap<speed_th) = 1;
    speed_lap(speed_lap>=speed_th) = 0;
    
    if debug
        figure(lap)
        plot(speed_lap), hold on
        plot(speed(idx_stop(1):idx_stop(2))./max(speed(idx_stop(1):idx_stop(2))),'r')
    end
    
    if debug
        figure(100)
        plot(X_lap{lap}, Y_lap{lap}, 'color', color(lap,:),...
            'displayname',sprintf('Lap %d',lap))
        hold on       
    end
    %extract regions in which the animal is still
    dist        = diff(speed_lap);
    moved       = -find(dist==1);
    stoped      = find(dist==-1);
    period      = [stoped  moved(1:length(stoped))];
    %select those stoppig periods larger than 1s
    winners     = find(sum(period,2) > 1.0*Fs);
    
    for w = 1:length(winners)
        idx_stop    = [-period(winners(w),2) period(winners(w),1)] + idx_stop(1); 
        
        %spikes
        n_spks_stop  = zeros(1,n_cells);
        n_spks_run   = zeros(1,n_cells);
        c_neu        = 0;
        for neu=1:n_cells
            if isIntern(neu) == 0
                c_neu  = c_neu + 1;

                t_stop = spk_lap{lap,neu}>=idx_stop(1) & spk_lap{lap,neu}<=idx_stop(2);
                %aligned to the start of the section
                n_spks_stop(neu) = sum(t_stop);
                if n_spks_stop(neu) ~= 0
                    t_spks_stop{c_neu} = spk_lap{lap,neu}(t_stop) - idx_stop(1) + 1;
                end

                t_run = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
                n_spks_run(neu) = sum(t_run);
                if n_spks_stop(neu) ~= 0
                    %aligned to the start of the section
                    t_spks_run{c_neu} = spk_lap{lap,neu}(t_run) - idx_run(1) + 1;
                end
            end
        end
        S(cnt).t_spks_stop = t_spks_stop;
        S(cnt).t_spks_run  = t_spks_run;
        S(cnt).n_spks_stop = n_spks_stop(isIntern==0);
        S(cnt).n_spks_run  = n_spks_run(isIntern==0);
        S(cnt).TrialId     = lap;
        S(cnt).TrialType   = typetrial{data.Laps.TrialType(laps(lap))};
        S(cnt).TrialTypeNo = data.Laps.TrialType(laps(lap));
        S(cnt).Interval    = idx_stop;
        S(cnt).Dur_stop    = sum(period(winners(w),:));
        S(cnt).Dur_run     = idx_run(end) - idx_run(1);
        S(cnt).Delay       = -period(w,2);
        S(cnt).speed_stop  = speed(idx_stop(1):idx_stop(2));
        S(cnt).speed_run   = speed(idx_run(1):idx_run(2));
        cnt                = cnt + 1;
        clear t_* n_sp*
    end    
end    
%%
%shows move and stop sequences for control purposes. Saves png of first lap 
t_ids      = [S.TrialId];
unique_tr    = unique(t_ids);
for seq = 1 : length(unique_tr)
   seq_r    = unique_tr(seq); 
   n_events = sum(t_ids == unique_tr(seq));
   figure(seq)
   set(gcf, 'position', [2000 500 1000 500], 'color', 'w'), hold on
   
   %theta
   subplot(1,2 + n_events,[1 2]), hold on
   ele   = find(t_ids==seq_r,1);
   raster(S(ele).t_spks_run)
   plot(S(ele).speed_run./10,'r')
   xlabel(sprintf('Running period lap %d::%s',seq_r,S(ele).TrialType))
   %events
   for n = 1 : n_events
      subplot(1,2 + n_events,2+n), hold on
      raster(S(ele+n-1).t_spks_stop)
      plot(S(ele+n-1).speed_stop./10,'r')
      xlabel('Stopping period')
   end
   drawnow 
end
print(figure(1),[roots{animal} 'Example_Run_stop_spks_lap.png'],'-dpng')

%%
%=========================================================================%
%==============   Extract super vector of spikes     =====================%
%=========================================================================%





