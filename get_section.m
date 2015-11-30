function D = get_section(D, in, out, debug, name)
%GET_SECTION add a new field to the struct D with the firing rate of the
%cells in a section specifie by in, out. D is generated from extract_laps.m
%
%
%ruben pinzon

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
    else
        disp('not a recognized section. options are [mid_arm, preturn, turn, lat_arm, reward, delay, wheel]')
    end
    eval(['sect_' sufix{s} '= int;'])
end 

% Extract firing rates in the section given
for lap = 1:numLaps  
    
    idx_lap      = [sum(D(lap).sections(sect_in,1)), sum(D(lap).sections(sect_out,2))];
    X_lap        = D(lap).X(idx_lap(1):idx_lap(2));
    Y_lap        = D(lap).Y(idx_lap(1):idx_lap(2));
    acc_dst      = cumsum(sqrt((X_lap - X_lap(1)).^2 + (Y_lap - Y_lap(1)).^2));
    speed_lap    = D(lap).speed(idx_lap(1):idx_lap(2));
    
    t_lap        = idx_lap(2) - idx_lap(1) + 1;
    firing       = D(lap).firing_rate(:,idx_lap(1):idx_lap(2)); 
    spk_train    = D(lap).spike_train(:,idx_lap(1):idx_lap(2)); 

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
    
    %Type of trial
    eval(['D(lap).' name '_firing=firing;'])
    eval(['D(lap).' name '_spike_train=spk_train;'])
    eval(['D(lap).' name '_interval=idx_lap;'])
    eval(['D(lap).' name '_speed=speed_lap;'])

    
end    