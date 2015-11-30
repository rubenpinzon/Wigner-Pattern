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
isIntern        = data.Clu.isIntern;
numLaps         = numel(laps)-1;
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
n_cells         = size(spk_lap,2);
n_pyrs          = sum(isIntern==0);
color           = hsv(n_pyrs);
conditions      = {'_left', '_right', ''};
% =======================================================================%
%==============   Extract Running Sections        ========================%
%=========================================================================%
debug           = true; %to show diganostic plots
%this is to remove/add the section in the middle arm of the maze
sect_in         = 5:6; 
sect_out        = 5:6;
kernel          = gausswin(0.1*Fs);
color_lap       = hsv(numLaps);
t_lap_max       = 2 * Fs;

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
        if t_lap < t_lap_max
            plot(X_lap, Y_lap, 'color', color_lap(lap,:),...
                'displayname',sprintf('Lap %d',lap))
            hold on   
        end

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
    
    if ~debug 
       figure(10+lap) 
       for n = 1 : n_pyrs
          plot(firing(seq(n),:)+n*ones(1,t_lap),'color',color_lap(lap,:)), hold on
           
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
    D(lap).mu_speed           = mean(speed_lap);
    D(lap).max_speed          = max(speed_lap);
    D(lap).min_speed          = min(speed_lap);
    D(lap).mu_acc             = mean(diff(speed_lap));
    D(lap).max_acc            = max(diff(speed_lap));
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

mean_speed = mean([D.mu_speed]); % m/sample
st_speed   = std([D.mu_speed]);

crit1 = [D.mu_speed] < (mean_speed+st_speed);
crit2 = [D.mu_speed] > (mean_speed-st_speed);
crit3 = [D.min_speed] > 0;


keep_laps = find( crit1 & crit2 & crit3 );
fprintf('Cleaning exp by removing %d out of %d trials with speed out of trend\n',length(D)-length(keep_laps),length(D))
D         = D(keep_laps);
numLaps   = length(D);
if debug
    
    for lap = 1:numLaps
        seq    = D(lap).f_rate_order;
        t_lap  = D(lap).duration +1;
        firing = D(lap).firing_rate;
        figure(1)
        plot(D(lap).X, D(lap).Y, 'color', color_lap(lap,:),...
                'displayname',sprintf('Lap %d',lap))
            hold on
        
        figure(10+lap) 
        for n = 1 : n_pyrs
          plot(firing(seq(n),:)+n*ones(1,t_lap),'color',color_lap(lap,:)), hold on           
        end
    end
end

%% temporal aligment for the place fields across laps

close all
psth_w              = 0.02 * Fs;
lap_types           = [D.type];
crit_pastalkova1    = @(x) max(x)>4 ; %from [Pastalkova 2008] Internally Generated Cell Assembly 
crit_pastalkova2    = @(x) max(x)>3 && max(x) > 2*std(x)+mean(x); %from [Pastalkova 2008] Internally Generated Cell Assembly 
overwrite           = false;

for con = 1 : max(lap_types)
   lap_same      = find(lap_types == con);
   
   if length(lap_same) > 1
       figure(con) 

       lap_durations = [D(lap_same).duration];
       mu_duration   = mean(lap_durations);
       sd_duration   = std(lap_durations);

       
       %longest lap
       max_lap          = max([D(lap_same).duration]);
       pastalkova_1     = false(1,n_pyrs);  
       pastalkova_2     = false(1,n_pyrs);    
       mu_pfield        = zeros(n_pyrs, max_lap);
       
       for n = 1 : n_pyrs
           subplot(10,ceil(n_pyrs/10),n)
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
           if  crit_pastalkova1(tmp(~isnan(tmp)))
              pastalkova_1(n) = true; 
              tmp_color = 'b';
           end   
           if crit_pastalkova2(tmp(~isnan(tmp)))
              pastalkova_2(n) = true;
              tmp_color = 'm'; 
              if strcmp(tmp_color, 'r')
                tmp_color = 'g'; 
              end
           end

           plot(mu_pfield(n,:), tmp_color, 'linewidth',2)
           xlim([0 max_lap])
           ylim([0 10])
           
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

stable_cells = C(1).pastalkova2...
    | C(2).pastalkova2 ;
stable_pcells = find(stable_cells==1)';
%%
%Filter out the cells that are not stable from the probe sequence
for lap = 1:numLaps  
    t_peak                  = D(lap).t_peak;
    t_peak(stable_cells==0) = -1;
    [val, seq]              = sort(t_peak); 
    probe_seq               = seq(val>0);
    D(lap).prob_seq         = probe_seq;
    clear *seq 
end
type        = [D.type];
probe_seq   = [D.prob_seq];
probe_seq_l = probe_seq(:,type==1);
probe_seq_r = probe_seq(:,type==2);

figure()
set(gcf,'color','w')
subplot(1,2,1)
errorbar(mean(probe_seq_l,2), std(probe_seq_l,1,2),'linewidth',2), hold on
for lp_l = 1 : size(probe_seq_l,1)
    plot(lp_l,probe_seq_l(lp_l,:),'bx')
end
xlabel('Cell id in the sequence')
ylabel('Cell Id')
xlim([0, sum(stable_cells)])
set(gca,'fontsize',14)
subplot(1,2,2)
errorbar(mean(probe_seq_r,2), std(probe_seq_r,1,2),'color','r','linewidth',2), hold on
for lp_r = 1 : size(probe_seq_r,1)
    plot(lp_r,probe_seq_r(lp_r,:),'rx')
end
xlabel('Cell id in the sequence')
ylabel('Cell Id')
xlim([0, sum(stable_cells)])
set(gca,'fontsize',14)
annotation('textbox',[0 0.9 1 0.1], 'String',animals{animal},...
    'fontsize',14,'edgecolor','none',...
    'HorizontalAlignment', 'center')

%check time variability
T_peak = [D.t_peak];
T_peak = T_peak(stable_cells==1,:);
T_peak_l = T_peak(:,type==1);
T_peak_r = T_peak(:,type==2);

mu_peak   = mean(T_peak_l,2);
sd_peak   = std(T_peak_l,1,2);
[val, idx]= sort(mu_peak);

herrorbar(mu_peak(idx),1:sum(stable_cells),sd_peak(idx))


%save to tx file the cells with stable pfields
if ~exist([roots{animal} '_stable_pfields.mat'], 'file') || overwrite
    save([roots{animal} '_stable_pfields.mat'],'D','C','stable_cells')
end


%% plot the firing rates in sequence according to the time of their peak
% How stable are the orderings across the laps?%
st_pcells = find(stable_cells==1);


for lap = 1:4%numLaps
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


