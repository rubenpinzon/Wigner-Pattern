function [D_left, D_right] = split_trails(D)
%SPLIT_TRIALS Splits the struct D containing all the trials into left and
%             righ alternations

D_left = struct('data',[],'condition','','epochColors',[],'trialId',[],'y',[],'T',[]);
D_right = struct('data',[],'condition','','epochColors',[],'trialId',[],'y',[],'T',[]);
for ilap = 1 : length(D)
    if strcmp(D(ilap).condition, 'left')
        D_left(end+1) = D(ilap);
    else
        D_right(end+1) = D(ilap);
    end    
end

D_left(1) = [];
D_right(1) = [];

