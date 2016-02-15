function R = get_section(D, in, out, debug, name)
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
        disp('WARNING: section not a recognized section. options are [mid_arm, preturn, turn, lat_arm, reward, delay, wheel]')
        return;
    end
    eval(['sect_' sufix{s} '= int;'])
end 


% Extract firing rates in the section given
for lap = 1:numLaps  
    
    idx_lap = [sum(D(lap).sections(sect_in,1)), sum(D(lap).sections(sect_out,2))];
    if int == 13
    %for wheel section extract spikes when the wheel is moving    
        wheelNonZero    = find(D(lap).wh_speed~=0);
        if isempty(wheelNonZero)
            fprintf('Skipped lap %d without wheel run\n',lap)
            return
        end
        idx_lap         = [wheelNonZero(1), wheelNonZero(end)];
        eval(['D(lap).' name '_wheelNonZero=wheelNonZero;'])
    end
    
    X_lap        = D(lap).X(idx_lap(1):idx_lap(2));
    Y_lap        = D(lap).Y(idx_lap(1):idx_lap(2));
    acc_dst      = cumsum(sqrt((X_lap - X_lap(1)).^2 + (Y_lap - Y_lap(1)).^2));
    speed_lap    = D(lap).speed(idx_lap(1):idx_lap(2));

    t_lap        = idx_lap(2) - idx_lap(1) + 1;
    %firing       = D(lap).firing_rate(:,idx_lap(1):idx_lap(2)); 
    spk_train    = D(lap).spike_train(:,idx_lap(1):idx_lap(2)); 

    if debug
        figure(1)
        subplot(121)
        hold on
        if int == 13
%             plot(D(lap).wh_speed, 'color', color(lap,:),...
%                 'displayname',sprintf('Lap %d',D(lap).trialId))
            plot(D(lap).wh_speed(wheelNonZero),'color', color(lap,:),...
                'displayname',sprintf('Lap %d',D(lap).trialId))
            xlabel('Samples'), ylabel('Wheel speed')
        else
            plot(X_lap, Y_lap, 'color', color(lap,:),...
                'displayname',sprintf('Lap %d',D(lap).trialId))
            xlabel('X position'), ylabel('Y position')
        end
        subplot(122)
        plot(speed_lap, 'color', color(lap,:),'displayname',sprintf('Lap %d',D(lap).trialId))
        xlabel('Samples'), ylabel('Animal speed')
        hold on   

    end              
    if int == 13
       eval(['R(lap).' name '_angularVelo=D(lap).wh_speed(wheelNonZero);'])
    end
    %Type of trial
    %eval(['R(lap).' name '_firing=firing;'])
    R(lap).trialId = D(lap).trialId;   
    R(lap).type    = D(lap).type; 
    eval(['R(lap).' name '_spike_train=spk_train;'])
    eval(['R(lap).' name '_interval=idx_lap;'])
    eval(['R(lap).' name '_speed=speed_lap;'])
    eval(['R(lap).' name '_position=[X_lap, Y_lap];'])

    
end