%Branch 2d SCRIPT FOR IDENTIFYING PROTO EVENTS
%
%   USAGE: (1) set the value of the global variables. (2) Extract the stopping
%   region after an alternation which corresponds to the reward area and store
%   all critical information (e.g., spikes, speed) in the struct S. (3) For control
%   purposes the raster during the running section and speed of the animal is shown
%   along the stopping periods extracted. (4) The super vector of activity during
%   the stopping periods is plotted. The super vector is constructed by joining
%   the spike times of all neurons in a single vector. Then silent periods of at
%   lest 50ms are detected to extract proto events. (5) Shows the proto events
%   extracted (6) Compute P(proto_event|model) with the model trained during the
%   running section, and classify each proto events according to the maximum likelihood
%   (7) Here shuffling controls are performed to test the values of loglikehood found.
%   the shuffling is made in time (bins) and contribution in the GPFA model (Cell ids).
%   (8) Calculates Loglike(P(run|model)) to have a baseline performance (9) Study the
%   run trials wrongly classified (10) Time shuffling the run trials to assess the
%   time structure of this data. (11) show proto events in latent space of running sections
%   (12) Change scale of GPFA heuristically and evaluate the classification scores.
%
%   DESCRIPTION:
%   Here the stopping periods after running are analyzed following the
%   procedure described in Foster&Wilson 2006 Neuron 36, A spike train was
%   constituted from all spikes (from all cells in the probe sequence) that
%   occurred during stopping periods while the animal faced in the direction
%   in which it had just run. This spike train was then broken between every
%   pair of successive spikes separated by more than 50 ms, to form a large
%   set of proto-events. Those proto-events in which at least one-third of
%   the cells in the probe sequence fired at least one spike were then
%   selected as events. The few events longer than 500 ms in duration were
%   rejected as a potential source of spurious correlations.
%
%Version 1.0 Ruben Pinzon@2015

clc, close all; clear all;

% ========================================================================%
%==================== (1) Variables of Interest ===========================
% ========================================================================%

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);

animal          = 6;
disp(['Processing folder: ' animals{animal}]);
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
n_laps          = numel(laps)-1;
[spk, spk_lap]  = get_spikes(clusters, data.Spike.res,laps);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
n_cells         = size(spk_lap,2);
color           = jet(sum(isIntern==0));
removeInh       = true;
show_run_stop   = true; %diagnostic plot to show running and stoping spikes
TrialType       = data.Laps.TrialType;
debug           = true; %to show diganostic plots
speed_th        = 100;
%this is to remove/add the section in the middle arm of the maze
sect            = 1; %without middle arm
sect_in         = [7, 8];
sect_out        = [7, 8];
cnt             = 1;
debug           = true;
t_window        = 0.05 * Fs;
t_max           = 0.50 * Fs;
n_mincell       = 10; %minimm number of cells active to be

%%
% ========================================================================%
%============ (2) Extract Stopping section after run =====================%
%=========================================================================%

% Extract spks when the mouse is running and in the wheel to calculate
%read the file with the stable place cells if exist
try
    data_pfields    = load([roots{animal} '_stable_pfieldsN.mat']);
    stable_pfields  = data_pfields.stable_cells;
    %probe_seq       = [data_pfields.D.f_rate_order];
    n_pyr           = length(stable_pfields);
    color           = jet(n_pyr); 
    fprintf('Found pfields file with %d cells\n',n_pyr)
catch
    disp('no stable pfields file found. Using all the cells except inter');
    stable_pfields = find(isIntern==0);
    n_pyr           = sum(isIntern==0);

end

for lap = 1:n_laps  

    idx_run                 = [sum(events{lap}(3:4,1)), sum(events{lap}(5:6,2))];

    idx_stop                = [sum(events{lap}(sect_in,1)), sum(events{lap}(sect_out,2))];
    X_lap{lap}              = X(idx_stop(1):idx_stop(2));
    Y_lap{lap}              = Y(idx_stop(1):idx_stop(2));
    speed_lap               = speed(idx_stop(1):idx_stop(2));
    
    %speed below threshold
    speed_lap(speed_lap<speed_th) = 1;
    speed_lap(speed_lap>=speed_th) = 0;
    speed_lap(end) = 0;
    norm_speed = max(speed(idx_stop(1):idx_stop(2)));
    if debug
        figure(lap)
        plot(speed_lap), hold on
        plot(speed(idx_stop(1):idx_stop(2))./norm_speed,'r')
        title(sprintf('Periods of inmobility (<%d), lap %d',speed_th, lap))
    end
    
    
    %extract regions in which the animal is still
    dist        = diff(speed_lap);
    moved       = find(dist==1);
    stoped      = find(dist==-1);
    per_stop    = [moved(1:length(stoped)) stoped];
    %select those stoppig periods larger than 1s
    winners     = find(per_stop(:,2)-per_stop(:,1) > 1.0*Fs);
    
    
    
    for w = 1:length(winners)
        idx_proto    = [per_stop(winners(w),1) per_stop(winners(w),2)] + idx_stop(1); 
        
        if debug
            plot(per_stop(winners(w),1):per_stop(winners(w),2),...
            speed(idx_proto(1):idx_proto(2))./norm_speed,'k','linewidth',2)        
        end
        
        
        %spikes
        n_spks_stop  = zeros(1,n_pyr);
        n_spks_run   = zeros(1,n_pyr);
        all_spks     = []; 
        all_spks_id  = []; 
        t_spks_stop  = cell(1,n_pyr);
        t_spks_run   = cell(1,n_pyr);
        for n=1:n_pyr
            neu = stable_pfields(n);

            t_stop = spk_lap{lap,neu}>=idx_proto(1) & spk_lap{lap,neu}<=idx_proto(2);
            %aligned to the start of the section
            n_spks_stop(n) = sum(t_stop);
            t_spks_stop{n} = spk_lap{lap,neu}(t_stop) - idx_proto(1) + 1;
            all_spks(end+1:end + n_spks_stop(n))    = t_spks_stop{n};
            all_spks_id(end+1:end + n_spks_stop(n)) = n;

            if debug
                for x = 1 : n_spks_stop(n)
                    line(t_spks_stop{n}(x)*[1 1] + per_stop(winners(w),1),[n n+0.8]./n_pyr,'color',color(n,:)) 
                end
            end

            t_run = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
            n_spks_run(n) = sum(t_run);
            t_spks_run{n} = spk_lap{lap,neu}(t_run) - idx_run(1) + 1;
            
        end
        
        if ~isempty(t_spks_stop)        
            S(cnt).spks_stop   = t_spks_stop;
            S(cnt).spks_run    = t_spks_run;
            S(cnt).n_spks_stop = n_spks_stop;
            S(cnt).n_spks_run  = n_spks_run;
            S(cnt).TrialId     = lap;
            S(cnt).TrialType   = typetrial{data.Laps.TrialType(laps(lap))};
            S(cnt).TrialTypeNo = data.Laps.TrialType(laps(lap));
            S(cnt).Interval    = idx_proto;
            S(cnt).Dur_stop    = sum(per_stop(winners(w),:));
            S(cnt).Dur_run     = idx_run(end) - idx_run(1);
            S(cnt).Delay       = -per_stop(w,2);
            S(cnt).speed_stop  = speed(idx_proto(1):idx_proto(2));
            S(cnt).speed_run   = speed(idx_run(1):idx_run(2));
            S(cnt).all_spks_stop = [all_spks; all_spks_id];
            S(cnt).inlap_pos   = [per_stop(winners(w),1) per_stop(winners(w),2)];
            cnt                = cnt + 1;
        end
        clear t_r* t_s* n_sp* all*
    end    
    drawnow
end

fprintf('Number of stoping events (speed<%d) = %d\n',speed_th,length(S))
fprintf('Next: searching for protoevents\n')
%% 
%========================================================================%
%========  (3) Shows move and stop sequences for control purposes========%
%========               Saves png of first lap                  =========%
%========================================================================%
if show_run_stop
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
       raster(S(ele).spks_run)
       plot(S(ele).speed_run./50,'r')
       xlabel(sprintf('Running period lap %d::%s',seq_r,S(ele).TrialType))
       ylim([0 40])
       ylabel('Speed/50')
       %events
       for n = 1 : n_events
          subplot(1,2 + n_events,2+n), hold on
          raster(S(ele+n-1).spks_stop)
          plot(S(ele+n-1).speed_stop./10,'r')
          xlabel('Stopping period')
          ylim([0 40])
       end
       

       drawnow 
    end
    print(figure(1),[roots{animal} 'Example_Run_stop_spks_lap.png'],'-dpng')
end
%%
%=========================================================================%
%============== (4)  Extract super vector of spikes     ==================%
%=========================================================================%

%50 ms window to detect disconnected sequences accoding to Foster & Wilson
%60 ms window according to Diba Buszaki 2007
%at least 10 neurons active 

                                   %replay preplay considered a event
cnt      = 1;
for seq = 1 : length(S)
    %super spike
    [s_spk,idx]    = sort(S(seq).all_spks_stop(1,:));
    s_spk(2,:)     = S(seq).all_spks_stop(2,idx);
    
    spk_seq        = S(seq).spks_stop;
    %show super vector  
    if debug
        figure(S(seq).TrialId)
        for s = 1 : length(s_spk)
           line(s_spk(1,s)*[1 1]+S(seq).inlap_pos(1),[0.6 0.8],'color',color(s_spk(2,s),:),'linewidth',2)         
        end
    end
    % calculate distances between spikes in super vector and show those
    %larger than 50 ms;
    dist_s     = diff(s_spk(1,:));
    dist_proto = [1 find(dist_s>t_window) numel(dist_s)];
    
    if ~isempty(dist_proto)
        
        
        for p = 1 : length(dist_proto)
            proto_int(:,p) = [s_spk(1,dist_proto(p)) s_spk(1,dist_proto(p)+1)] + S(seq).inlap_pos(1);
            if debug
                line(proto_int(:,p),[0.7 0.7],'color','k','linewidth',2) 
                line(proto_int(1,p)*[1 1],[0.5 0.9],'color','k','linewidth',2)
                line(proto_int(2,p)*[1 1],[0.5 0.9],'color','k','linewidth',2)
            end
        end
        %add proto events of the head and tail of the stop periods
        
        
        %proto event has to be withing 500 ms. % odd values are disntace
        %within event
        len_proto = diff(proto_int(:)); 
        %interval that fullfils the criteria
        st_pnt_silent  = proto_int(1:2:end);
        en_pnt_silent  = proto_int(2:2:end);
        p_eve          = find(len_proto(2:2:end)<t_max);   

        for p = 1 : length(p_eve)

            proto(p,:)    = [en_pnt_silent(p_eve(p)) st_pnt_silent(p_eve(p)+1)] - S(seq).inlap_pos(1);
            idx_proto     = find(s_spk(1,:)>=proto(p,1) & s_spk(1,:)<=proto(p,2));
            cell_proto{p} = unique(s_spk(2,idx_proto));   

            k      = 'm';
            size_l = [0.52 0.58];
            n_cell_proto = length(cell_proto{p});
            %Here are the proto events that fullfill the criteria
            if  n_cell_proto > n_mincell
                k = 'r';
                size_l   = [0 1];
                P(cnt).duration = proto(p,2)-proto(p,1);
                P(cnt).interval = proto(p,:);
                P(cnt).trialId  = S(seq).TrialId;
                P(cnt).trialType = S(seq).TrialTypeNo;
                
                spks_times  =  s_spk(1,idx_proto) - proto(p,1)+1;
                spk_pro     = zeros(n_pyr,P(cnt).duration + 1); %proto event
                spk_raw     = {};
                for cc = 1 : n_cell_proto
                    c     = cell_proto{p}(cc);
                    s_idx = spks_times(find(s_spk(2,idx_proto)==c));
                    spk_pro(c,s_idx) = 1;
                    spk_raw{c}       = s_idx;
                end
                P(cnt).data     = spk_pro;
                P(cnt).spk_raw  = spk_raw;
                P(cnt).spk_cnt  = sum(P(cnt).data,2);
                cnt             = cnt + 1;
            end
            if debug
                line([en_pnt_silent(p_eve(p)) st_pnt_silent(p_eve(p)+1)],[0.55 0.55],'color',k,'linewidth',1) 
                line(en_pnt_silent(p_eve(p))*[1 1],size_l,'color',k,'linewidth',2)
                line(st_pnt_silent(p_eve(p)+1)*[1 1],size_l,'color',k,'linewidth',2)
                text(en_pnt_silent(p_eve(p)), 0.50,sprintf('%d',n_cell_proto),'fontsize',9)
            end
        end

        %at least 30% of cells active in proto event to be consider event
        clear proto_int len_proto 
        drawnow
    end
end
%save([roots{animal} '_proto_events_v2.mat'],'P','S')

%% 
%========================================================================%
%========= (5) Shows move and proto events for control purposes =========%
%=========                Saves png of first lap              ===========%
%========================================================================%

if show_run_stop
    t_ids      = [P.trialId];
    unique_tr    = unique(t_ids);
    for seq = 1 : length(unique_tr)
       seq_r    = unique_tr(seq); 
       n_events = sum(t_ids == unique_tr(seq));
       figure(seq)
       set(gcf, 'position', [2000 500 1500 500], 'color', 'w'), hold on

       %theta
       subplot(1,2 + n_events,[1 2]), hold on
       ele   = find(t_ids==seq_r,1);
       raster(S(seq).spks_run)
       plot(S(seq).speed_run./10,'r')
       xlabel(sprintf('Running period lap %d::%s',seq_r,S(seq).TrialType))
       %events
       for n = 1 : n_events
          subplot(1,2 + n_events,2+n), hold on
          raster(P(ele+n-1).spk_raw)
          xlabel('stop. per.')
       end
       drawnow 
    end
    print(figure(1),[roots{animal} 'Example_Run_stop_spks_lap.png'],'-dpng')
end

%%
%=========================================================================%
%=========(6) Load GPFA models from the running laps =====================%
%=========    Compute loglike P(proto_event|model)   =====================%
%=========================================================================%

load([roots{animal} '_branch2_results40ms.mat'])
load([roots{animal} '_proto_events_v2.mat'])

models      = {result_D_left result_D_right};
colors      = hsv(2);
name_field  = 'data';
n_folds     = 3;
debug       = false;
P_seg       = segment(P, 0.004, Fs, keep_cell, name_field, 0); 
n_proto     = length(P_seg);
n_mods      = length(models);

%Classification stats of P(proto_event|model) 
stats       = classGPFA(P_seg, models);
cm          = [stats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))
%save([roots{animal} '_confusionM_protoEvents.mat'],'stats')
%%
%=========================================================================%
%=========    (6.5) Loglike changing binning         =====================%
%=========================================================================%
bin_size = 0.004:0.001:0.07;

for ibin = 1:length(bin_size)
   P_seg   = segment(P(1), bin_size(ibin), Fs, keep_cell, name_field); 
   stats   = classGPFA(P_seg, n_folds, debug, models);
   
   datalike{ibin} = reshape([stats.likelihood],n_mods * n_folds, n_proto);   
   
   
   cm          = [stats.conf_matrix];
   fprintf('bin= %1.3f ave.T = %2.2f | hitA: %2.2f%%, hitB: %2.2f%%\n',bin_size(ibin),...
       mean([P_seg.T]), 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))
   confmat{ibin} = cm;
end

figure(65), hold on
for ibin = 1:length(bin_size)   
    plot(bin_size(ibin)*ones(n_folds,n_proto), datalike{ibin}(1:2:end,:),'.r')
    plot(bin_size(ibin)*ones(n_folds,n_proto), datalike{ibin}(2:2:end,:),'ob')
    drawnow
end

%=========================================================================%
%=========    (6.6) Loglike changing tiem scale GP   =====================%
%=========================================================================%
scalingGP   = linspace(0.01, 5, 200);
P_seg       = segment(P, 0.01, Fs, keep_cell, name_field); 
P_seg       = P_seg(6);
n_proto     = length(P_seg);
colors      = jet(length(scalingGP));
for isca = 1:length(scalingGP)
   stats   = classGPFA(P_seg, n_folds, debug, models, scalingGP(isca));   
   datalike{isca} = reshape([stats.likelihood],n_mods * n_folds, n_proto);
   lv = stats(1).traj{1}.xsm(1,:);
   plot3(1:numel(lv),scalingGP(isca)*ones(1,numel(lv)),lv,'color',colors(isca,:)), hold on  
%    
%    cm          = [stats.conf_matrix];
%    fprintf('Scale= %1.3f | hitA: %2.2f%%, hitB: %2.2f%%\n',scalingGP(isca),...
%        100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))
%    confmat{isca} = cm;
    drawnow
end

figure(68), hold on
for isca = 1:length(datalike)  
    plot(scalingGP(isca), ((datalike{isca}(1,1))),'.r')
    plot(scalingGP(isca), ((datalike{isca}(2,1))),'ob')
    drawnow
end

%%
%=========================================================================%
%=========(7) Shuffle controls time and cell IDs     =====================%
%=========================================================================%

% 500 Suffled cells id proto events
n_perms    = 100;
for s = 1 : n_perms
    [P_perm,idx] = shuffcells(P_seg);
    
    stats = classGPFA(P_perm, n_folds, debug, models);
    
    cm    = [stats.conf_matrix];
    CM(s).conf_matrix = cm;
    CM(s).shuffledidx = idx;
    fprintf('Shuffle %d-th, hitA: %2.2f%%, hitB: %2.2f%%\n',s, 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))
end
save([roots{animal} '_confusionM_protoEvents_cellShuffle.mat'],'CM')

% 100 Shuffled times id proto events
for s = 1 : n_perms

    P_perm = shufftime(P_seg);
    
    stats = classGPFA(P_perm, n_folds, debug, models);
    
    cm    = [stats.conf_matrix];
    CM_time(s).conf_matrix = cm;
    CM_time(s).stats = stats;
    fprintf('Shuffle %d-th, hitA: %2.2f%%, hitB: %2.2f%%\n',s, 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))
    clear P_perm
end
save([roots{animal} '_confusionM_protoEvents_timeShuffle.mat'],'CM_time')

%%
%=========================================================================%
%====== (8) Calculate the LogLike(P(run|model)) for test    ==============%
%======      running trails to get the baseline             ==============%
%=========================================================================%


load([roots{animal} '_branch2_results40ms.mat'])

mod_tags    = {'_left', '_right'};
posterior   = [];
n_folds     = 3;
mask        = false(1,length(D)); % for cross validation


for ifold = 1 : n_folds    
    
    loglike_folds   = [];

    for mod = 1 : length(mod_tags)
        stest_mask      = mask;
        model_data      = eval(['result_D' mod_tags{mod}]);
        cv_trials       = model_data.cv_trials;
        fold_indx       = model_data.foldidx;
        test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
        test_data       = D(test_mask);
        type            = ones(1, length(test_data));
    
        for ilap = 1 : length(test_data) %for each trial independently
            [~, loglike] = exactInferenceWithLL(test_data(ilap), model_data.params{ifold},'getLL',1);
            loglike_folds(mod, ilap) = loglike;    

            if strcmp(test_data(ilap).condition, 'right')
                type(ilap)  = 2;
            end
        end      
        
    end
    
    %Classification rates
    [~, best_mod] = max(loglike_folds);
    TP            = sum(best_mod == 1 & type == 1)/(sum(type == 1));
    FN            = sum(best_mod == 2 & type == 2)/(sum(type == 2));
    FP            = sum(best_mod == 1 & type == 2)/(sum(type == 2));
    TN            = sum(best_mod == 2 & type == 1)/(sum(type == 1));
    
    stats(ifold).conf_matrix    = [TP, FP; TN, FN];
    stats(ifold).class_output   = best_mod;
    stats(ifold).real_label     = type;
    stats(ifold).loglike        = loglike_folds;
    stats(ifold).trailsId       = [test_data.trialId];
end
cm    = [stats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))

%========================================================================%
%======== (9) Study the trials wrongly classified  =======================
%========================================================================%

err_cla     = [];
for ifold = 1 : n_folds
    err_cla =  [err_cla stats(ifold).trailsId(stats(ifold).class_output ~= stats(ifold).real_label)];    
end
err_cla = unique(err_cla);

showPos(X,Y, events, err_cla)
%%
%========================================================================%
%======== (10) Time shuffling the running data (D)  ======================
%========================================================================%

n_perms   = 100;
clear time_shuff best_mod loglike_folds
for sf = 1 : n_perms
    fprintf('Permutation %d out of %d...',sf, n_perms)
    
    for ifold = 1 : n_folds  
        for mod = 1 : length(mod_tags)
            test_mask       = mask;
            model_data      = eval(['result_D' mod_tags{mod}]);
            cv_trials       = model_data.cv_trials;
            fold_indx       = model_data.foldidx;
            test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
            test_data       = D(test_mask);
            type            = ones(1, length(test_data));
            test_shuff      = shufftime(test_data);
            for ilap = 1 : length(test_data) %for each trial independently
                [~, loglike] = exactInferenceWithLL(test_shuff(ilap), model_data.params{ifold},'getLL',1);
                loglike_folds(mod, ilap) = loglike;    

                if strcmp(test_shuff(ilap).condition, 'right')
                      type(ilap)  = 2;
                end
            end
        end
        %Classification rates
        [~, best_mod] = max(loglike_folds);
        TP            = sum(best_mod == 1 & type == 1)/(sum(type == 1));
        FN            = sum(best_mod == 2 & type == 2)/(sum(type == 2));
        FP            = sum(best_mod == 1 & type == 2)/(sum(type == 2));
        TN            = sum(best_mod == 2 & type == 1)/(sum(type == 1));

        stats(ifold).conf_matrix    = [TP, FP; TN, FN];
        stats(ifold).class_output   = best_mod;
        stats(ifold).real_label     = type;
        stats(ifold).loglike        = loglike_folds;
        stats(ifold).trailsId       = [test_data.trialId];
        clear loglike_* best_mod*

    end
    cm    = [stats.conf_matrix];
    fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))

    time_shuff(sf).stats     = stats;
    clear stats loglike*
end

%%
%========================================================================%
%======(11) show proto events in latent space of running sections =======%
%========================================================================%

model = result_D.params{1};
for i = 1 : length(P_seg)
    
    [traj,ll]  = exactInferenceWithLL(P_seg(i), model,'getLL',1);
    Xorth = orthogonalize([traj.xsm], model.C);

    if stats(1).class_output(i) == 1
        c = 'cyan';
    else
        c = 'r';
    end

   % plot3(Xorth(1,:),Xorth(2,:),Xorth(3,:),'color',c)
    plot3(Xorth(1,end),Xorth(2,end),Xorth(3,end),'color',c,'marker','s','markerfacecolor',c)
    plot3(Xorth(1,1),Xorth(2,1),Xorth(3,1),'color',c,'marker','o','markerfacecolor',c, '')

    hold on
end
%%
%========================================================================%
%=======   (12) Change scale of GPFA heuristically and evaluate ===========
%=======              classification performance               ===========
%========================================================================%

scale = linspace(0.01, 10, 100);

for sc = 1 : length(scale) 
    
    stats = classGPFA(P_seg, n_folds, debug, models, scale(sc));
    cm    = [stats.conf_matrix];
    CM_scale(sc).conf_matrix = cm;
    CM_scale(sc).scale = scale(sc);
    fprintf('Scale %f hitA: %2.2f%%, hitB: %2.2f%%\n',scale(sc), 100*mean(cm(1,[1,3,5])),100*mean(cm(2,[2,4,6])))
    
end

%% Show likelihood of classified proto events

load([roots{animal} '_confusionM_protoEvents_timeShuffle.mat'])

% stats for fold # 1
for r = 1 : n_perms
    group     = CM_time(r).stats(1).class_output;
    for p = 1 : length(group)
        loglike_CM(r,p)  = CM_time(r).stats(1).posterior(group(p),p);
    end
    groups(r,:) = group;
    clear group==
end
p         = anova1(loglike_CM(:),groups(:)) 

plot(group, CM, 'x')