% Script to create figures of Estimated Intrinsic Dimension (EID)


%window 2 secods
% 
bin_sizes = 0.03:0.01:0.1;
%
figure(1)
set(gcf,'color','w', 'position', [3378, 1, 500, 973])

%window 3 secods
load('i01_maze06.002_run_w3s_results.mat')
% results(9)=[];%erase last cell because itnot from w=3s
for i= 1:length(results)
    subplot(length(results),1,i)
    [~, foldmax] = max(sum(results(i).like));
    [~, imax] = max(results(i).like(:,foldmax));
    meanLike = mean(results(i).like,2);
    stdLike = std(results(i).like,0,2);
    plot(meanLike), hold on, box off
    set(gca, 'ycolor', 'w')
    plot(imax, meanLike(imax), 'bo')
    line([imax imax],[(meanLike(imax) - stdLike(imax))...
                    (meanLike(imax) + stdLike(imax))])
    text(30, min(meanLike),['bin = ' num2str(1000*bin_sizes(i)) 'ms']);
    text(1, min(meanLike),['EID = ' num2str(imax)],'color','b');
end
xlabel('Dimension')
subplot(length(results),1,1)
title('Cross-validated (N=3), EID, win=3s, Sect Maze 2to5/6')
%%
load('i01_maze06.002_run_w2s_results.mat')
results(1)=[];%erase bin=20ms because not runned in w=3s
for i= 1:length(results)
    subplot(length(results),1,i)
    [~, foldmax] = max(sum(results(i).like));
    [~, imax] = max(results(i).like(:,foldmax));
    meanLike = mean(results(i).like,2);
    stdLike = std(results(i).like,0,2);
    plot(meanLike), hold on, box off
    set(gca, 'ycolor', 'w')
    plot(imax, meanLike(imax), 'ro')
    line([imax imax],[(meanLike(imax) - stdLike(imax))...
                    (meanLike(imax) + stdLike(imax))])   
    text(20, min(meanLike),['EID = ' num2str(imax)],'color','r');

end