%% Analysis of Buzsaki database i01_maze05_MS.005
% L/R alternation task, 18 min, 15sec wheel runs, several errors


clear all; clc, close all;

basepath        = '/media/bigdata/';
animals          ={'i01_maze06.002', 'i01_maze13.003', 'i01_maze08.001',...
                   'i01_maze15.002', 'i01_maze08.004', 'i01_maze06.005',...
                    'i01_maze05.005'};

files = get_matFiles(animals, basepath);

for rat = 1 : length(animals)
    r_obj           = load([basepath  files{rat}]);
    fprintf('Saving files at %s\n', files{rat})
    clusters        = r_obj.Spike.totclu;
    mazesect        = r_obj.Laps.MazeSection;
    wheelspeed      = r_obj.Laps.WhlSpeedCW;
    events          = r_obj.Par.MazeSectEnterLeft;
    speed           = r_obj.Track.speed;
    eeg             = r_obj.Track.eeg;
    isIntern        = r_obj.Clu.isIntern;
    spike_times     = r_obj.Spike.res;
    Fs              = 1250;
    
    %
    cell_type = zeros(length(clusters),1);
    for clus = 1 : length(clusters)
        cell_type(clus) = 1 + (~isIntern(clusters(clus)));    
    end

    areas = [32, 64, 128]; %EC, right and left CA1
    area_type = zeros(length(clusters),1);
    for clus = 1 : length(clusters)
        area = 3;
        if clusters(clus)<32
            area = 1;
        elseif clusters(clus)<64
            area = 2;
        end
        area_type(clus) = area;
    end

    section_spike = zeros(length(clusters),1);
    wheel_speed_spike = zeros(length(clusters),1); 
    for spi = 1 : length(spike_times)
        section_spike(spi) = mazesect(spike_times(spi));
        wheel_speed_spike(spi) = wheelspeed(spike_times(spi));

    end
    %
    header = [basepath files{rat}(1:end-19)];
    save4attila(r_obj.Laps.TrialType,[header 'trialType'],'')
    save4attila([spike_times/Fs*1000 clusters area_type cell_type section_spike wheel_speed_spike], [header 'spikes_aumented'], '%5.2f\t%3d\t%1d\t%1d\t%2d\t%3.3f\n')
    save4attila(eeg,[header 'eeg'],'')
    save4attila(mazesect,[header 'sections'],'')
    save4attila(wheelspeed,[header 'wheelspeed'],'')
    save4attila(speed,[header 'mousespeed'],'')
    save4attila([r_obj.Clu.isIntern], [header 'isIntern'], '%d\n')
    clear r*
end