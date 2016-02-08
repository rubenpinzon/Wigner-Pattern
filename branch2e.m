%BRANCH 2E Script to find the pyramidal cells with stable place fields according
%       to Pastalkova, E., Itskov, V., Amarasingham, A., & Buzsáki, G. (2008).
%       Internally generated cell assembly sequences in the rat hippocampus.
%       Science, 321(5894), 1322-1327.
%
%
%Ruben Pinzon@2015

clc, close all; clear all;



%========================Variables of Interest===========================

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
isIntern        = data.Clu.isIntern;
numLaps         = length(events);
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
n_cells         = size(spk_lap,2);
n_pyrs          = sum(isIntern==0);
TrialType       = data.Laps.TrialType;
Typetrial_tx    = {'left', 'right', 'errorLeft', 'errorRight'};

% ========================================================================%
%==============       (1)   Extract lap           ========================%
%=========================================================================%

D = extract_laps(Fs,spk_lap,speed,X,Y,events,isIntern, laps, TrialType);

% ========================================================================%
%==============   (2) Find stable pcells          ========================%
%=========================================================================%

debug       = true;
in          = 'turn';
out         = 'lat_arm';
min_speed   = 0;
show_firing = false;
[E, Ev]      = get_pfields(D, in,out, min_speed, 'easy', debug, show_firing);

%Get consolitated indixes of stable place cells
stable_E = [];
for e = 1 : length(E)
    stable_E = [stable_E E{e}];    
end
% stable_pcells = unique(stable_E); %these are the cells with stable pfield
stable_cells = E{1};
n_sta        = length(stable_cells);

%% Get stopping periods

in_reward  = 'reward';
out_reward = 'reward';
Stop       = get_section(D, in_reward, out_reward, debug, 'reward');
speed_th   = 150;
color      = jet(n_sta);
cnt        = 1;

for lap = 1:length(Stop) 
    spd_lap   = Stop(lap).reward_speed;
    speed_lap = Stop(lap).reward_speed;
    spd_lap(spd_lap<speed_th) = 1;
    spd_lap(spd_lap>=speed_th) = 0;
    spd_lap(end) = 0;
    norm_speed = max(speed_lap);
    if debug
        figure(lap)
        plot(norm_speed.*spd_lap), hold on
        plot(speed_lap,'r')
        title(sprintf('Periods of inmobility (<%d), lap %d',speed_th, lap))
    end
    
    %extract regions in which the animal is still
    dist        = diff(spd_lap);
    moved       = find(dist==1);
    stoped      = find(dist==-1);
    per_stop    = [moved(1:length(stoped)) stoped];
    %select those stoppig periods larger than 1s
    winners     = find(per_stop(:,2)-per_stop(:,1) > 1.0*Fs);
    gain        = norm_speed/n_sta;
    for w = 1:length(winners)
        
        idx = [per_stop(winners(w),1) per_stop(winners(w),2)];
        if debug
            plot(per_stop(winners(w),1):per_stop(winners(w),2),...
                speed_lap(idx(1):idx(2)),'k','linewidth',2)
        end
                
        %spikes inside the stoping areas
        n_spks_stop  = zeros(1,n_sta);
        n_spks_run   = zeros(1,n_sta);
        n            = 0;
        all_spks     = []; 
        all_spks_id  = []; 
        t_spks_stop  = cell(1,n_sta);
        t_spks_run   = cell(1,n_sta);
        for n = 1 : n_sta
            neu = stable_cells(n);
            
            t_stop = find(Stop(lap).reward_spike_train(neu,idx(1):idx(2))==1);
            %aligned to the start of the section
            n_spks_stop(n) = length(t_stop);
            t_spks_stop{n} = t_stop;
            all_spks(end+1:end + n_spks_stop(n))    = t_spks_stop{n};
            all_spks_id(end+1:end + n_spks_stop(n)) = neu;

            if debug
                for x = 1 : n_spks_stop(n)
                    line(t_spks_stop{n}(x)*[1 1] + per_stop(winners(w),1),gain*[n n+0.8],'color',color(n,:)) 
                end
            end            
        end        
        if ~isempty(t_spks_stop)        
            Ev(cnt).spks_stop   = t_spks_stop;
            Ev(cnt).spks_run    = t_spks_run;
            Ev(cnt).n_spks_stop = n_spks_stop;
            Ev(cnt).n_spks_run  = n_spks_run;
            Ev(cnt).TrialId     = lap;
            Ev(cnt).TrialType   = Stop(lap).type;
            Ev(cnt).all_spks_stop = [all_spks; all_spks_id];
            Ev(cnt).inlap_pos   = [per_stop(winners(w),1) per_stop(winners(w),2)];
            cnt                = cnt + 1;
        end
        clear t_* n_sp* all*        
    end   
      
    drawnow
end
