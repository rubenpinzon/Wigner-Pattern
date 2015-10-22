%Using the trained model with left, right both alternations test the SPWs


clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);


%% =======================================================================%
%======== (1) Load the models from Theta oscillations    =================%
%=========================================================================%
animal          = 6;
data            = load(files{animal});
laps            = data.Laps.StartLaps(data.Laps.StartLaps~=0); %@1250 Hz
laps(end+1)     = data.Par.SyncOff;
isIntern        = data.Clu.isIntern;
[spk, spk_lap]  = get_spikes(data.Spike.totclu, data.Spike.res,laps);
conditions      = {'_left', '_right', ''};
Fs              = 1250;
debug           = false;
load([roots{animal} '_branch2_results40ms.mat']);
[D_left, D_right] = split_trails(D);
ifold           = 1;
n_cells         = size(spk_lap,2);
%% =======================================================================%
%======== (2) Get SPWs in the reward or wheel area       =================%
%=========================================================================%

spw_tag = {'after_left', 'after_right', 'after_left_error', 'after_right_error'};
markers = 1.25*load([roots{animal} '_spws.txt']); %from ms to samples
updateGP = true;

if debug
    figure(20)
    set(gcf, 'position', [0 1 1000 300], 'color', 'w')
    plot(eeg./max(eeg),'color',[0.8, 0.8, 0.8]), hold on
    plot(0.08*mazesect-0.5,'b')
    ylim([-0.6 0.6])
    plot(repmat(markers(:,1),1,2),ylim,'r')
    plot(repmat(markers(:,2),1,2),ylim,'c')
    plot(0.1*data.Laps.TrialType,'r')
end

color   = jet(4);
for sp = 1 : length(markers)   
    spw_type{sp}     = spw_tag{data.Laps.TrialType(ceil(markers(sp,1)))+1};
    color_spw(sp,:)  = color(data.Laps.TrialType(ceil(markers(sp,1)))+1,:);
    %spikes during the SPW    
    idx_run                 = ceil(markers(sp,1):markers(sp,2));
    spw_len(sp)            = (idx_run(end)-idx_run(1))/Fs;
    cnt = 1;
    %sect 1:enter, 6:exit
    for neu=1:n_cells
        if isIntern(neu)==0
            idx = spk{neu}>=idx_run(1) & spk{neu}<=idx_run(end);
            %aligned to the start of the section
            SpkSPW{sp,cnt} = spk{neu}(idx) - idx_run(1) + 1;
            cnt = cnt +  1;
            if debug
                plot(repmat(spk{neu}(idx),1,2)', ...
                    repmat(0.01*[neu neu+0.9],...
                    length(spk{neu}(idx)),1)', 'linewidth',2)
                ylim([0 1.2])
            end
        end
    end
end

bin_size        = 0.004; %ten times smaller than Theta
SpkSPW_DH       = get_high(SpkSPW(:,keep_cell==1), ceil(spw_len*Fs),...
                     spw_type, color_spw, 'spw', 0);
F               = segment(SpkSPW_DH, 0.004, Fs, 0.01); %4ms bin size

%% =======================================================================%
%========        (3) LogLike p(spw|model)                =================%
%=========================================================================%
ll_spw_given_model = zeros(length(F), length(conditions));

for model = 1 : length(conditions)
    params = eval(['result_D' conditions{model} '.params{' num2str(ifold) '};']);
    for spw = 1 : length(F)
        %remove all but one SPW
        X = F(spw);
        [traj, ll]   = exactInferenceWithLL(X, params,'getLL',1);
        ll_spw_given_model(spw, model) = ll;
        fprintf('Condition %s, SPW %d/%d\n', conditions{model}, spw, length(F))
    end
end

[a,b] = max(ll_spw_given_model');
for spw = 1 : length(F)
    fprintf('Max log like for spw %d with model %s\n',spw,conditions{b(spw)})
end

%spws segmented accoring to the max LogLike
for model = 1 : length(conditions)
    eval(sprintf('spws%s = F(b(b==%d));',conditions{model}, model))
end

%plots
for model = 1 : length(conditions)
    figure(model)
    params = eval(['result_D' conditions{model} '.params{' num2str(ifold) '};']);
    eval(sprintf('spws%s = F(b(b==%d));',conditions{model}, model))
    [traj, ll_te]   = exactInferenceWithLL(spws{1}, params,'getLL',1);
    [Xorth, Corth] = orthogonalize([traj.xsm], params.C);

end