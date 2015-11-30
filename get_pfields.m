function [stable, S] = get_pfields(D, in, out, min_speed, method, debug, show_firing)
%GET_PFIELDS This function process a data array of spike times to indetify
%the cells with stable place fields. The criteria used to define a stable
%place field follows the definitions at
%
%Pastalkova, E., Itskov, V., Amarasingham, A., & Buzsáki, G. (2008). 
%Internally generated cell assembly sequences in the rat hippocampus.
%Science, 321(5894), 1322-1327.
%
%Villette, V. et al. ernally recurring hippocampal
%sequences as a population template of spatiotemporal information.
%Neuron, 88(2), 357-366.
%
%Criteria: a minimal peak firing rate of 6.0 Hz or 5Hz and peak firing rate
%being at least 4.5-times or 3-times SD above a mean firing rate
%
%Inputs: D          :HC-5 data in a structure splited by laps
%        in         :refers to the name of the enter section in the maze
%        out        :refers to the name of the exit section in the maze
%        min_speed  :minimum speed allowed for the mouse
%        method     :'easy' for selection of cells based on cosistency of
%                     firig onset. 'pasta' to follow pastalkova criteria
%        debug
%
%@Ruben


numLaps  = length(D);
color    = hsv(numLaps);

sufix    = {'in', 'out'};
for s = 1 : length(sufix)
    sect = eval(sufix{s});
    if strcmp(sect,'mid_arm')
        int = 1;
    elseif strcmp(sect,'preturn')
        int = 2;
    elseif strcmp(sect,'turn')
        int = 3:4;
    elseif strcmp(sect,'lat_arm')
        int = 5:6;
    elseif strcmp(sect,'reward')
        int = 7:8;
    elseif strcmp(sect,'delay')
        int = 9:12;
    elseif strcmp(sect,'wheel')
        int = 13;
    end
    eval(['sect_' sufix{s} '= int;'])
end 

% Extract spks when the mouse is running 
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    
    idx_lap      = [sum(D(lap).sections(sect_in,1)), sum(D(lap).sections(sect_out,2))];
    X_lap        = D(lap).X(idx_lap(1):idx_lap(2));
    Y_lap        = D(lap).Y(idx_lap(1):idx_lap(2));
    acc_dst      = cumsum(sqrt((X_lap - X_lap(1)).^2 + (Y_lap - Y_lap(1)).^2));
    speed_lap    = D(lap).speed(idx_lap(1):idx_lap(2));
    
    t_lap        = idx_lap(2) - idx_lap(1) + 1;
    cnt          = 0;
    firing       = D(lap).firing_rate(:,idx_lap(1):idx_lap(2)); 
    spikes       = D(lap).spike_train(:,idx_lap(1):idx_lap(2));
    if debug
        figure(1)
        subplot(121)
        hold on
        plot(X_lap, Y_lap, 'color', color(lap,:),...
            'displayname',sprintf('Lap %d',D(lap).trialId))
        subplot(122)
        plot(speed_lap, 'color', color(lap,:),'displayname',sprintf('Lap %d',D(lap).trialId))
        hold on   

    end
        
    %cell ordered in time according to the onset of firing activity
    %For each recurring activity pattern, cell activation onset was defined by
    %the maximum of the first derivative of the smoothed trace. Onsets
    %Cells were then ordered in a sequence according to their median onset
    %delay in each run epoch. Villette, V. et al. ernally recurring hippocampal
    %sequences as a population template of spatiotemporal information.
    %Neuron, 88(2), 357-366.

    f               = diff(firing,1,2);
    [peak, t_peak]  = max(f,[],2);
    [t_peak_s,seq]  = sort(t_peak);  
    
    if sum(isnan(speed_lap))~= 0
        disp('NaN in speed found')
    end
    
    %Type of trial
    S(lap).onset              = f;
    S(lap).f_rate_order       = seq;
    S(lap).t_peak             = t_peak;
    S(lap).X                  = X_lap;
    S(lap).Y                  = Y_lap;
    S(lap).speed              = speed_lap;
    S(lap).mu_speed           = mean(speed_lap);
    S(lap).max_speed          = max(speed_lap);
    S(lap).min_speed          = min(speed_lap);
    S(lap).mu_acc             = mean(diff(speed_lap));
    S(lap).max_acc            = max(diff(speed_lap));
    S(lap).firing_rate        = firing;
    S(lap).duration           = idx_lap(2) - idx_lap(1);
    S(lap).type               = D(lap).type;
    S(lap).trialId            = D(lap).trialId;
    S(lap).spike_train        = spikes;
    
end    

%#TODO Used the percentile 95 of the lenght and speed to filter laps
mean_speed = mean([S.mu_speed]); % m/sample
st_speed   = std([S.mu_speed]);

%criteria to remove trials that are out of trend
crit1 = [S.mu_speed] < (mean_speed+st_speed);
crit2 = [S.mu_speed] > (mean_speed-st_speed);
crit3 = [S.min_speed] > min_speed;


keep_laps = find( crit1 & crit2 & crit3 );
fprintf('Cleaning exp by removing %d trials out of %d\n',numLaps-length(keep_laps),numLaps)
S         = S(keep_laps);
numLaps   = length(S);
n_pyrs    = size(S(1).firing_rate,1);

%plot ordered firing rates
if show_firing    
    for lap = 1:numLaps
        seq    = S(lap).f_rate_order;
        t_lap  = S(lap).duration +1;
        firing = S(lap).firing_rate;
        figure(2)        
        plot(S(lap).X, S(lap).Y, 'color', color(lap,:),...
                'displayname',sprintf('Lap %d',lap))
            hold on
        
        figure(10+lap) 
        for n = 1 : n_pyrs
          plot(firing(seq(n),:)+n*ones(1,t_lap),'color',color(lap,:)), hold on           
        end
    end
end

% ========================================================================%
%==============             Get stable cells      ========================%
%=========================================================================%

%easy way is to compute the median of the ofset of the activity, and then
%remove the cells that are not concentrated, that is, large s.d. Second
%methods is by following Lee and Wilson: pyramidal cells firing > spike/lap
%(mean), having a single place field(thus unique place field), and
%that continued firing during the RUN session in more that 2/3 of the trials
%(consistency).
type          = [S.type];
max_var_onset = 350;

if strcmp(method,'easy')      

    %check time variability based on Villette, V. et al. Neuron, 88(2) 2015
    probe_seq = [S.t_peak];
    for t = 1 : max(type)  
       
        T_peak    = probe_seq(:,type==t);
        if length(T_peak)>3
            figure(100+t)
            mu_peak   = mean(T_peak,2);
            sd_peak   = std(T_peak,1,2);
            [val, idx]  = sort(mu_peak);
            ori_idx = idx;
            remove_cell = val<10 | sd_peak(idx) > max_var_onset;
            idx(remove_cell) = [];
            herrorbar(mu_peak(idx),1:numel(idx),sd_peak(idx));            
            title(sprintf('Laps type %d (num trials = %d)',t, sum(type==t)))
            xlabel('Onset time')
            ylabel('Cell ordered (filtered)')
            set(gca,'FontSize',14)
            
            figure(110+t)
            herrorbar(mu_peak(ori_idx),1:numel(ori_idx),sd_peak(ori_idx));            
            title(sprintf('Laps type %d (num trials = %d)',t, sum(type==t)))
            xlabel('Onset time')
            ylabel('Cell ordered')
            set(gca,'FontSize',14)
        end
        
        stable{t} = idx';
    end     
elseif strcmp(method,'lee')
    % place fields according to Lee and Wilson, Neuron(36), 2002 
    trials = [S.type];

    for t = 1 : max(type) 
        trials_t = find(trials==t);
        n_trials = numel(trials_t);        
        
        %criterion A: mean firing of 1 spike/lap
        min_fire  = n_trials; % 
        spike_cnt = sum([S.spike_train],2);
        critA     = find(spike_cnt >= min_fire);
        
        %criterion B: single peak
        max_duration = max([S.duration]);

        
        for cm = 1 : numel(critA)
            c = critA(cm);
            tmp = zeros(n_trials,max_duration);
            for k = 1 : n_trials
                lap = trials_t(k);
                firing_c     = S(lap).firing_rate;
                if S(lap).duration ~= max_duration
                    f_cell = interp1(firing_c(c,:),linspace(0,S(lap).duration,max_duration));
                else
                    f_cell = firing_c(c,1:max_duration);
                end
                tmp(k,:) = f_cell;
                %plot(f_cell), hold on

            end
            
                             
            [dip,xl,xu, ifault] = HartigansDipTest(mean(tmp(:,50:end)));
            

            figure(c)
            plot(mean(tmp(:,50:end)),'k','linewidth',2), hold on
            title (sprintf('Dip stats cell %d = %f',c,dip))

            drawnow
        end
        
    end
end


%plot ordered firing rates
if show_firing    
    for lap = 1:numLaps
        seq      = stable{S(lap).type};
        t_lap    = S(lap).duration +1;
        firing   = S(lap).firing_rate; 
        n_stable = length(seq);
        onset    = S(lap).t_peak;
        
        figure(200+lap) 
        for n = 1 : n_stable
          plot(firing(seq(n),:)+n*ones(1,t_lap),'color',color(lap,:)), hold on 
          t  = onset(seq(n));
          line(t*[1, 1],[n  0.8+n],'color','k','linewidth',3)
        end
        title(sprintf('Lap %d',S(lap).trialId))

    end
end
