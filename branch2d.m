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

%%
% ========================================================================%
%==============   Extract Stopping section after run =====================%
%=========================================================================%

debug           = false; %to show diganostic plots
speed_th        = 50;
%this is to remove/add the section in the middle arm of the maze
sect            = 1; %without middle arm 
sect_in         = [7, 8]; 
sect_out        = [7, 8];
cnt             = 1;
% Extract spks when the mouse is running and in the wheel to calculate

%read the file with the stable place cells if exist
try
    data_pfields    = load([roots{animal} '_stable_pfields.mat']);
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
    if debug
        for w = 1:length(winners)
            idx = [per_stop(winners(w),1) per_stop(winners(w),2)] + idx_stop(1);
            plot(per_stop(winners(w),1):per_stop(winners(w),2),...
                speed(idx(1):idx(2))./norm_speed,'k','linewidth',2)
        end
    end
    
    
    for w = 1:length(winners)
        idx_stop    = [per_stop(winners(w),1) per_stop(winners(w),2)] + idx_stop(1); 
        
        %spikes
        n_spks_stop  = zeros(1,n_pyr);
        n_spks_run   = zeros(1,n_pyr);
        n        = 0;
        all_spks     = []; 
        all_spks_id  = []; 
        t_spks_stop  = cell(1,n_pyr);
        t_spks_run   = cell(1,n_pyr);
        for n=1:n_pyr
            neu = stable_pfields(n);

            t_stop = spk_lap{lap,neu}>=idx_stop(1) & spk_lap{lap,neu}<=idx_stop(2);
            %aligned to the start of the section
            n_spks_stop(n) = sum(t_stop);
            t_spks_stop{n} = spk_lap{lap,neu}(t_stop) - idx_stop(1) + 1;
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
            S(cnt).Interval    = idx_stop;
            S(cnt).Dur_stop    = sum(per_stop(winners(w),:));
            S(cnt).Dur_run     = idx_run(end) - idx_run(1);
            S(cnt).Delay       = -per_stop(w,2);
            S(cnt).speed_stop  = speed(idx_stop(1):idx_stop(2));
            S(cnt).speed_run   = speed(idx_run(1):idx_run(2));
            S(cnt).all_spks_stop = [all_spks; all_spks_id];
            S(cnt).inlap_pos   = [per_stop(winners(w),1) per_stop(winners(w),2)];
            cnt                = cnt + 1;
        end
        clear t_* n_sp* all*
    end    
    drawnow
end

fprintf('Number of stoping events (speed<%d) = %d\n',speed_th,length(S))
fprintf('Next: searching for protoevents\n')
%% 
%========================================================================%
% Shows move and stop sequences for control purposes. Saves png of first lap
%=========================================================================%
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
%==============   Extract super vector of spikes     =====================%
%=========================================================================%

%50 ms window to detect disconnected sequences accoding to Foster & Wilson
%60 ms window according to Diba Buszaki 2007
%at least 10 neurons active 

t_window = 0.05 * Fs;
t_max    = 0.5 * Fs;
n_mincell= 10; %minimm number of cells active to be 
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
    dist_proto = find(dist_s>t_window);
    
    if ~isempty(dist_proto)
        
        for p = 1 : length(dist_proto)
            proto_int(:,p) = [s_spk(1,dist_proto(p)) s_spk(1,dist_proto(p)+1)] + S(seq).inlap_pos(1);
            if debug
                line(proto_int(:,p),[0.7 0.7],'color','k','linewidth',2) 
                line(proto_int(1,p)*[1 1],[0.5 0.9],'color','k','linewidth',2)
                line(proto_int(2,p)*[1 1],[0.5 0.9],'color','k','linewidth',2)
            end
        end
        %proto event has to be withing 500 ms. % odd values are disntace
        %within event
        len_proto = diff(proto_int(:)); 
        %interval that fullfils the criteria
        st_pnt_silent  = proto_int(1:2:end);
        en_pnt_silent  = proto_int(2:2:end);
        p_eve          = find(len_proto(2:2:end)<t_max);   

        for p = 1 : length(p_eve)

            proto(p,:) = [en_pnt_silent(p_eve(p)) st_pnt_silent(p_eve(p)+1)];
            idx_proto  = find(s_spk(1,:)>=proto(p,1) & s_spk(1,:)<=proto(p,2));
            cell_proto{p} = unique(s_spk(2,idx_proto));   

            k = 'm';
            size_l = [0.52 0.58];
            %Here are the proto events that fullfill the criteria
            if length(cell_proto{p}) > n_mincell
                k = 'r';
                size_l   = [0 1];
                P(cnt).duration = proto(p,2)-proto(p,1);
                P(cnt).interval = proto(p,:);
                P(cnt).trialId  = S(seq).TrialId;
                P(cnt).trialType = S(seq).TrialTypeNo;
                
                spk_pro = zeros(n_pyr,P(cnt).duration + 1); %proto event
                spk_raw = {};
                for c = 1 : n_pyr
                    s_idx = find(spk_seq{c}>=proto(p,1) & spk_seq{c}<=proto(p,2));
                    spk_pro(c,s_idx) = 1;
                    spk_raw{c} = s_idx;
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
            end
        end

        %at least 30% of cells active in proto event to be consider event
        clear proto_int len_proto 
%         text(sum(xlim)/4,1,sprintf('Silent periods >=%2.1f ms',1000*t_window/Fs))
%         text(sum(xlim)/2,1,'Events (>30% cells active)','color','r')
%         text(3*sum(xlim)/4,1,'Proto Events (<30% cells active)','color','m')
        drawnow
    end
end
save([roots{animal} '_proto_events.mat'],'P','S')

%% 
%========================================================================%
% Shows move and stop sequences for control purposes. Saves png of first lap
%=========================================================================%

if show_run_stop
    t_ids      = [P.trialId];
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
       plot(S(ele).speed_run./10,'r')
       xlabel(sprintf('Running period lap %d::%s',seq_r,S(ele).TrialType))
       %events
       for n = 1 : n_events
          subplot(1,2 + n_events,2+n), hold on
          raster(P(ele+n-1).spk_raw)
          xlabel('Stopping period')
       end
       drawnow 
    end
    print(figure(1),[roots{animal} 'Example_Run_stop_spks_lap.png'],'-dpng')
end

%%
%=========================================================================%
%============ Load GPFA models from the running laps =====================%
%============ Compute loglike P(proto_event|model)   =====================%
%=========================================================================%

load([roots{animal} '_branch2_results40ms.mat'])
load([roots{animal} '_proto_events.mat'])

mod_tags    = {'_left', '_right'};
colors      = hsv(2);
name_field  = 'data';
P_seg       = segment(P, 0.004, Fs, keep_cell, name_field); 
n_folds     = 3;

for ifold = 1 : folds
    posterior   = zeros(length(mod_tags), length(P));
    proto_traj  = cell(length(mod_tags), length(P));


    for p = 1 : length(P) 
        for m = 1 : length(mod_tags)
            %select the model parameters from the fold#1 
            model = eval(['result_D' mod_tags{m} '.params{ifold};']);
            %model.gamma = model.gamma*10;
            [traj, ll] = exactInferenceWithLL(P_seg(p), model,'getLL',1);        
            posterior(m,p) = ll;
            proto_traj{m,p} = traj;      

        end       

        [~, max_mod]  = max(posterior(:,p));
        loglike_p(p)  = posterior(max_mod, p);

        %get joint neural space
        model = result_D.params{1};
        [traj, ll] = exactInferenceWithLL(P_seg(p), model,'getLL',1);

        figure(ifold)
    %     max_mod = P(p).trialType;
        Xorth = orthogonalize([traj.xsm], model.C);
    %     Xorth = normr(Xorth);

        plot3(Xorth(1,:),Xorth(2,:),Xorth(3,:),'color',colors(max_mod,:)), hold on
        plot3(Xorth(1,1),Xorth(2,1),Xorth(3,1),'color',...
            colors(max_mod,:),'markerfacecolor',colors(max_mod,:),'marker','s')
    %     figure(2)
    %     plot_latent(Xorth, colors(max_mod,:),num2str(p))
        drawnow
        fprintf('Proto event %d\n',p);

    end

    title('Classified according to trial ')
    set(gcf, 'color','w')
    grid on
    xlabel('X1'),ylabel('X2'),zlabel('X3')
    set(gca,'fontsize',14)
    %max log like P(event|model)
    [val, best_mod]  = max(posterior);
    type            = [P.trialType]; %{P(proto|model) , realtag}


    TP            = sum(best_mod == 1 & type == 1)/(sum(type == 1));
    FN            = sum(best_mod == 2 & type == 2)/(sum(type == 2));
    FP            = sum(best_mod == 1 & type == 2)/(sum(type == 2));
    TN            = sum(best_mod == 2 & type == 1)/(sum(type == 1));
    
    stats_proto(ifold).conf_matrix    = [TP, FP; TN, FN];
    stats_proto(ifold).class_output   = best_mod;
    stats_proto(ifold).real_label     = type;
    stats_proto(ifold).posterior      = posterior;

end



%% Cealculate the LogLike(P(run|model)) for test trails to get the baseline
%%of the classifier

load([roots{animal} '_branch2_results40ms.mat'])

mod_tags    = {'_left', '_right'};
posterior   = [];
n_folds     = 3;
mask        = false(1,length(D)); % for cross validation


for ifold = 1 : n_folds    
    test_mask       = mask;
    test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
    test_data       = D(test_mask);
    loglike_folds   = zeros(length(mod_tags), length(test_data));
    type            = ones(1, length(test_data));
    
    for mod = 1 : length(mod_tags)
        model_data = eval(['result_D' mod_tags{mod}]);
        cv_trials       = model_data.cv_trials;
        fold_indx       = model_data.foldidx;
    
        for ilap = 1 : length(test_data) %for each trial independently
            [~, loglike] = exactInferenceWithLL(test_data(ilap), model_data.params{ifold},'getLL',1);
            loglike_folds(mod, ilap) = loglike;
             

            if strcmp(test_data(ilap).condition, 'right') && mod == 1
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



