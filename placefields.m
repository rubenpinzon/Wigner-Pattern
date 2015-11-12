%Identifying the cells with defined place fields Place and episode fields.
%On each trial each neuron's spike train was convolved with a Gaussian (SD 
%100msec, sampling rate 1250 Hz). firing rate in a wheel as “episode fields”.
%Place fields on the maze and episode fields in the wheel were defined by 
%a minimal peak firing rate of 6.0 Hz or 5Hz and peak firing rate being at
%least 4.5-times or 3-times SD above a mean firing rate, respectively.
%Using these different criteria yielded similar results. The width
%From [Pastalkova 2008] Internally Generated Cell Assembly Sequences in the
%Rat Hippocampus

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
numLaps         = numel(laps)-1;
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
n_cells         = size(spk_lap,2);
n_pyrs          = sum(isIntern==0);
color           = jet(n_pyrs);
conditions      = {'_left', '_right', ''};
% =======================================================================%
%==============   Extract Running Sections        ========================%
%=========================================================================%
debug           = true; %to show diganostic plots
%this is to remove/add the section in the middle arm of the maze
sect_in         = 1; 
sect_out        = [5:6];
kernel          = gausswin(0.1*Fs);

% Extract spks when the mouse is running 
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    
    idx_lap      = [sum(events{lap}(sect_in,1)), sum(events{lap}(sect_out,2))];
    int_at_maze  = idx_lap;
    X_lap        = X(idx_lap(1):idx_lap(2));
    Y_lap        = Y(idx_lap(1):idx_lap(2));
    acc_dst      = cumsum(sqrt((X_lap - X_lap(1)).^2 + (Y_lap - Y_lap(1)).^2));
    speed_lap    = speed(idx_lap(1):idx_lap(2));
    
    %sect 1:enter, 6:exit
    t_lap        = idx_lap(2) - idx_lap(1) + 1;
    cnt          = 0;
    firing       = zeros(n_pyrs, t_lap); 
    for neu=1:n_cells
        if ~isIntern(neu)
            tmp          = zeros(1, t_lap); 
            cnt      = cnt + 1;
            
            idx = spk_lap{lap,neu}>=idx_lap(1) & spk_lap{lap,neu}<=idx_lap(end);
            %aligned to the start of the section            
            spikes{cnt} = spk_lap{lap,neu}(idx) - idx_lap(1) + 1;
            tmp(spikes{cnt}) = 1; 
            %convolve the spike trains with a gauss filter 100 ms
            firing(cnt,:) = conv(tmp,kernel, 'same');
        end
    end
    if debug
        figure(2)
        plot(X_lap, Y_lap, 'color', color(lap,:),...
            'displayname',sprintf('Lap %d',lap))
        hold on       
    end
    
    if debug && lap == 1
       figure(10) 
       for n = 1 : n_pyrs
          plot(firing(n,:)+n*ones(1,t_lap),'color',color(n,:)), hold on
           
       end
    end
    %Type of trial
    D(lap).spikes             = spikes;
    D(lap).X                  = X_lap;
    D(lap).Y                  = Y_lap;
    D(lap).speed              = speed_lap;
    D(lap).times              = idx_lap;
    D(lap).type               = data.Laps.TrialType(laps(lap));
    D(lap).color              = color(data.Laps.TrialType(laps(lap)),:);
    D(lap).acc_dist           = acc_dst;
    D(lap).firing_rate        = firing;
    D(lap).duration           = idx_lap(2) - idx_lap(1);
end    

if debug
   figure(2), title('Position of animal per Lap in section Run')
end

%% temporal aligment for the place fields across laps

psth_w    = 0.02 * Fs;
lap_types = [D.type];
for con = 2 : 2%length(typetrial)
   figure(con) 
   lap_same      = find(lap_types == con);
   lap_durations = [D(lap_same).duration];
   mu_duration   = mean(lap_durations);
   sd_duration   = std(lap_durations);
   
   
   %remove laps away from the mean duration
   remove_lap = find(abs(lap_durations - mu_duration)>1.5*sd_duration);
   lap_same(remove_lap) = [];
   %longest lap
   max_lap = max([D(lap_same).duration]);
   for n = 1 : n_pyrs
       subplot(10,9,n)
       tmp = zeros(1,max_lap);
       for j = 1 : numel(lap_same)
           lap          = lap_same(j);
           spike_train  = D(lap).firing_rate(n,:);
           if D(lap).duration ~= max_lap+1
               spike_train = interp1(spike_train,linspace(0,length(spike_train),max_lap));
           end
           tmp = tmp + spike_train;
           plot(spike_train,'color',color(j,:)), hold on       
       end
       plot(tmp./numel(lap_same), 'k', 'linewidth',2)
       xlim([0 max_lap])
   end
   
end
