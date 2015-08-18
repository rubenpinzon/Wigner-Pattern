
function [leftT, rightT, failed_trial] = protoLap(XT, YT, length_run, trial, X_Run_Lap, Y_Run_Lap, int_at_maze, Fs, animal, isIntern)
% PROTOLAP Calcualtes the typical (mean) trial of the animal give all the
%          trials succeded.


numLaps = length(length_run);
N       = size(X_Run_Lap,2);
longest = max(length_run)*Fs;
xl      = zeros(1,longest) ;yl = zeros(1,longest); %vectors to calculate the mean trajectory
xr      = zeros(1,longest) ;yr = zeros(1,longest); 
r       = 0;l = 0;

failed_trial = [];
for ilap = 1 : numLaps
    duration    = linspace(0,length_run(ilap),longest); 
    %real animal position and interporaltion
    x       = XT(int_at_maze(ilap,1):int_at_maze(ilap,2));
    y       = YT(int_at_maze(ilap,1):int_at_maze(ilap,2));   
    xi      = spline(linspace(0,length_run(ilap),length(x)),x,duration);
    yi      = spline(linspace(0,length_run(ilap),length(y)),y,duration);        
    if strcmp(trial{ilap}, 'right') 
        xr = xr + xi; yr = yr + yi; r  = r + 1;
        %count spikes in the grids of the right arm        
    elseif strcmp(trial{ilap}, 'left') 
        xl  = xl + xi;yl = yl + yi; l  = l + 1;       
    else
        failed_trial =[failed_trial ilap]; %#ok<AGROW>
    end
    plot(x, y), hold on
    for icell = 1 : N
        if ~isIntern(icell)
            plot(X_Run_Lap{ilap,icell}, Y_Run_Lap{ilap,icell},'x')
        end
    end

end
leftT= [xl; yl]'./l;
rightT= [xr; yr]'./r;

plot(leftT(:,1), leftT(:,2),'k','linewidth', 2), hold on
plot(rightT(:,1), rightT(:,2),'k','linewidth', 2), 
title(sprintf('Animal %s',animal))
xlabel('x')