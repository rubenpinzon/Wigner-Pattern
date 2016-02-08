function showPos(X,Y, events, filter_laps)
%SHOWPOS this function shows the position of the animal in the running
%        section of the maze. If a specific vector of laps in given in
%        filter only those laps are shown.
%
%        INPUTS:
%        X:         X position of the animal during all the experiment
%        Y:         Y position of the animal during all the experiment
%        events:    Times of starting ending of each lap
%        filter_laps:    Vector of laps to show in [1 - numLaps]
%
%Ruben Pinzon@2015
figure()
set(gcf,'color','w');
numLaps         = length(filter_laps);
color           = hsv(numLaps);
% Extract spks when the mouse is running 
for l = 1: numLaps 
    lap = filter_laps(l);
    
    idx_lap      = [events{lap}(1,1), events{lap}(end,2)];
    X_lap        = X(idx_lap(1):idx_lap(2));
    Y_lap        = Y(idx_lap(1):idx_lap(2));
    
    plot(X_lap, Y_lap, 'color', color(l,:),...
            'displayname',sprintf('Lap %d',lap))
        hold on
    
end
 