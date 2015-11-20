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
animal          = 4;
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
            firing(cnt,:) = Fs*conv(tmp,kernel, 'same');
        end
    end
    if debug
        figure(2)
        plot(X_lap, Y_lap, 'color', color(lap,:),...
            'displayname',sprintf('Lap %d',lap))
        hold on       
    end
    
    %cell ordered in time according to the peak in firing rate
    %For each recurring activity pattern, cell activation onset was defined by
    %the maximum of the first derivative of the smoothed trace. Onsets
    % Cells were then ordered in a sequence according to their median onset
    %delay in each run epoch. Villette, V. et al. ernally recurring hippocampal
    %sequences as a population template of spatiotemporal information.
    %Neuron, 88(2), 357-366.

    f               = diff(firing,1,2);
    [peak, t_peak]  = max(f,[],2);
    [t_peak_s,seq]  = sort(t_peak);
    
    if debug && lap < 5
       figure(10+lap) 
       for n = 1 : n_pyrs
          plot(firing(seq(n),:)+n*ones(1,t_lap),'color',color(seq(n),:)), hold on
           
       end
    end  
    
    
    %Type of trial
    D(lap).onset              = f;
    D(lap).f_rate_order       = seq;
    D(lap).t_peak             = t_peak;
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

close all
psth_w              = 0.02 * Fs;
lap_types           = [D.type];
crit_pastalkova1    = 4 ; %from [Pastalkova 2008] Internally Generated Cell Assembly 
crit_pastalkova2    = 3 ; %from [Pastalkova 2008] Internally Generated Cell Assembly 

for con = 1 : 2%length(typetrial)
   figure(con) 
   lap_same      = find(lap_types == con);
   
   if ~isempty(lap_same)
       lap_durations = [D(lap_same).duration];
       mu_duration   = mean(lap_durations);
       sd_duration   = std(lap_durations);

       %remove laps away from the mean duration
       remove_lap = find(abs(lap_durations - mu_duration)>1.5*sd_duration);
       lap_same(remove_lap) = [];
       %longest lap
       max_lap          = max([D(lap_same).duration]);
       pastalkova_1     = false(1,n_pyrs);  
       pastalkova_2     = false(1,n_pyrs);    
       mu_pfield        = zeros(n_pyrs, max_lap);
       
       for n = 1 : n_pyrs
           subplot(10,9,n)
           %figure(n)
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
           %title(sprintf('n%d',n))
           %ylabel('Hz')
           tmp            = tmp./numel(lap_same); 
           mu_pfield(n,:) = tmp;

           tmp_color = 'k';
           %criteria for place fields
           if  max(tmp) >= crit_pastalkova1
              pastalkova_1(n) = true; 
              tmp_color = 'r';
           end   
           if max(tmp) >= crit_pastalkova2 && max(tmp) > 2*std(tmp(~isnan(tmp))) + mean(tmp(~isnan(tmp)))
              pastalkova_2(n) = true;
              tmp_color = 'm'; 
              if strcmp(tmp_color, 'r')
                tmp_color = 'g'; 
              end
           end

           plot(mu_pfield(n,:), tmp_color, 'linewidth',2)
           xlim([0 max_lap])
           ylim([0 10])
           line(xlim, crit_pastalkova1*[1 1],'color',[0.5 0.5 0.5])
           
       end
       drawnow
       C(con).condition = typetrial{con};
       C(con).mu_pfield = mu_pfield;
       C(con).pastalkova1 = pastalkova_1;
       C(con).pastalkova2 = pastalkova_2;
       C(con).n_pfield = sum(pastalkova_1);  
       C(con).mean_dur_lap = mu_duration;
       C(con).sd_dur_lap = sd_duration;
       clear mu* sd* tmp*
   end

end

stable_cells = C(1).pastalkova1 | C(1).pastalkova2...
    | C(2).pastalkova2 | C(2).pastalkova1;
stable_pcells = find(stable_cells==1)';
%save to tx file the cells with stable pfields
fileID = fopen([roots{animal} '_stable_pfields.txt'],'w');
fprintf(fileID,'%d\n',stable_pcells);
fclose(fileID);
%% plot the firing rates in sequence according to the time of their peak


for lap = 1:numLaps
    firing = D(lap).firing_rate;
    seq    = D(lap).f_rate_order; 
    t_peak = D(lap).t_peak;
    
    figure(10+lap) 
    for n = 1 : n_pyrs
       if stable_cells(seq(n))
            plot(firing(seq(n),:)+n,'color',color(25,:)), hold on  
            plot(t_peak(seq(n)),n,'ok','markersize',6,'markerfacecolor','k')
       end           
    end   
end
