%% Analysis of Buzsaki database
clc, close all, clear all;

animals         = {'i01_maze06.002', 'i01_maze05.005', 'i01_maze06.005',...
                    'i01_maze08.001' };
basepath        = '/media/bigdata/';
files           = get_matFiles(animals, basepath);

%========================Variables of Interest===========================
obj             = load([basepath files{2}]);
clusters        = obj.Spike.totclu;
laps            = obj.Laps.StartLaps(obj.Laps.StartLaps~=0); %@1250 Hz
%Adding the end of the last lap because Laps only contains the start
laps(end+1)     = obj.Par.SyncOff;
mazesect        = obj.Laps.MazeSection;
wheelspeed      = obj.Laps.WhlSpeedCW;
X               = obj.Track.X;
Y               = obj.Track.Y;
events          = obj.Par.MazeSectEnterLeft;
Fs              = obj.Par.SamplingFrequency;
eeg             = obj.Track.eeg;
speed           = obj.Track.speed;
isIntern        = obj.Clu.isIntern;
numLaps         = numel(laps)-1;
time_seg        = linspace(0, length(eeg)/Fs, length(eeg));
[spk, spk_lap]  = get_spikes(clusters, obj.Spike.res,...
                                                laps);
N               = size(spk_lap,2);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
tcolor          = hsv(numLaps+1);
%
debug           = 1; %to show diganostic plots

% Extract spks when the mouse is running and in the wheel to calculate
for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    idx_lap = laps(lap):laps(lap+1);

    wheelNonZero        = find(wheelspeed(idx_lap)~=0) + laps(lap);
    wheelSpeed_lap{lap} = wheelspeed(wheelNonZero);
    if ~isempty(wheelNonZero)
        wheel_len(lap)   = numel(wheelNonZero)/Fs;

        if debug
            figure(1)
            x   = linspace(0, wheel_len(lap), numel(wheelSpeed_lap{lap})); 
            plot(x, wheelSpeed_lap{lap}, 'color', tcolor(lap,:)), hold on
            text(10, 100 + 20*lap, ['Lap ' num2str(lap)], 'color', tcolor(lap,:));
        end
        
        %last lap does not have wheel
        for neu=1:N
            idx = spk_lap{lap,neu}>=wheelNonZero(1)...
                          & spk_lap{lap,neu}<=wheelNonZero(end);
            %aligned to the start of the section
            SpkWheel_lap{lap,neu} = spk_lap{lap,neu}(idx) - wheelNonZero(1) + 1;
        end
    end
    
    %(b) Runing in the maze. Extracted based on the
    %EnterSection time stamp without considering left-right
    idx_run                 = [events{lap}(2,1), sum(events{lap}(5:6,2))];
    int_at_maze(lap, :)     = idx_run;
    run_len(lap)            = (idx_run(2)-idx_run(1))/Fs;
    X_at_maze               = X(idx_run(1):idx_run(2));
    Y_at_maze               = Y(idx_run(1):idx_run(2));
    speed_lap               = speed(idx_run(1):idx_run(2));

    %sect 1:enter, 6:exit
    for neu=1:N
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        %aligned to the start of the section
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx) - idx_run(1) + 1;
    end
    if debug
        figure(2)
        plot(X_at_maze, Y_at_maze, 'color', tcolor(lap,:)), hold on
        text(700, 400 + 30*lap, ['Lap ' num2str(lap)], 'color', tcolor(lap,:));
        figure(22)
        plot(speed_lap, 'color', tcolor(lap,:)), hold on
        text(700, 2000 + 200*lap, ['Lap ' num2str(lap)], 'color', tcolor(lap,:));
    end
    %Type of trial
    trial{lap}          = typetrial{obj.Laps.TrialType(laps(lap))};
    color(lap,:)        = tcolor(obj.Laps.TrialType(laps(lap)),:);
end
figure(1), title('Wheel speed per Lap')
figure(2), title('Position of animal per Lap in section Run')
figure(22), title('Wheel speed (encoder)')
if debug
   figure(3)
   plot(1:numLaps,run_len,'-s'), hold on
   plot(1:numLaps-1,wheel_len, '-o')
   xlabel('Lap Num.'), ylabel('time (s)')
   title('Duration of sections per lap')
   legend('Run', 'Wheel')
end
%% time-segmenting with the fastest lap
% Convert to DataHigh format without segmenting, that is, using the whole
% time that the animal spent in the runing section. This implies laps with
% different lenghts. The alternative is to segment based on spatial bins.

MaxTimeE_run     = floor(Fs * min(run_len) * ones(1, numLaps));
%MaxTimeE_run = floor(Fs * run_len);
onlyCorrectTrial = true;
%Data processed for datahigh without interneuorns
Run = get_high(SpkRun_lap(:,isIntern==0), MaxTimeE_run,...
                     trial, color, 'run', onlyCorrectTrial);

if debug
    figure(4)
    for ilap = 1: numLaps
        if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'left')
            idx_run                 = [events{ilap}(2,1), events{ilap}(2,1) + MaxTimeE_run(ilap)];
            X_at_maze               = X(idx_run(1):idx_run(2));
            Y_at_maze               = Y(idx_run(1):idx_run(2));
            plot(X_at_maze, Y_at_maze, 'color', tcolor(ilap,:)), hold on
            plot(X_at_maze(end), Y_at_maze(end), 's', 'color',...
                 tcolor(ilap,:), 'MarkerSize', 8, 'markerfacecolor',tcolor(ilap,:))
            text(700, 400 + 30*ilap, ['Lap ' num2str(ilap)], 'color', tcolor(ilap,:)); 
        end
    end
end

%with the wheel section the difference in laps can be seen in the distance
%spanned by the wheel encoder and it could be standarized to be uniform, or 
%just use time but then distance covered is different.

% MaxTimeE = floor(Fs * min(wheel_len) * ones(1, numLaps-1));   
MaxTimeE_wh    = floor(Fs * 10 * ones(1, numLaps-1));
Wheel          = get_high(SpkWheel_lap(:,isIntern==0), MaxTimeE_wh,...
                        trial, color, 'wheel',onlyCorrectTrial);

%% Command based GPFA based Alexander Aecker OOP implementation.
% Check dependency of the method with the bin size


Laps            = length(Run);
bin_size        = 0.1; 
bin_width       = ceil(bin_size * Fs); % bin size (Seconds * Fs) = samples
dims            = 3 : 2 : 30; % Target latent dimensions
results(1).bin  = bin_size;
firing_thr      = 0.2 ; % Minimum firing rate find which 
                        % neurons should be kept
Run             = filterbin(Run, firing_thr, bin_width);

%preallocating variables
%issues: Very low firing rates per lap. To achieve the estimation all the
%laps have to be taken into account, and therefore the temporal lenght has
%to be uniform which is a limiting factor in the running section. Once the
%model is trained, the projection of each lap data into the Latent space
%would produce the neural trajectories.
for ilap = 1 : Laps
    for idim = 1 : length(dims)   
        
            Y     = [Run.bins];
            Y     = reshape(Y, [55 Run(1).T Laps]);
            
            model = GPFA('Tolerance', 1e-6, 'verbose', 1);

            model = model.fit(sqrt(Y), 10 , 'hist');
            [model, Xest] = model.normFactors(Y);

            
        end
    end
end
