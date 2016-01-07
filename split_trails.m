function [D_left, D_right] = split_trails(D)
%SPLIT_TRIALS Splits the struct D containing all the trials into left and
%             righ alternations
%             version 1.0: build two new structs for left and right models
%             version 2.0: split the extracture D with all the fields
%             included.
%
%ruben


typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};

if isfield(D, 'condition')
    %for backwards compatibility
    disp('Condition field found. Splitting with vesion 1.0 of function.')

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
else
    
   types = [D.type];
   D_left  = D(types==1);
   D_right = D(types==2);
   
   %add colors and name of th e trial
   for d = 1 : length(D_left)       
      D_left(d).condition = typetrial{1};
      D_left(d).epochColors = [1 0 0];
   end
   for d = 1 : length(D_right)       
      D_right(d).condition = typetrial{2};
      D_right(d).epochColors = [0 0 1];
   end     
    
end
    
    

