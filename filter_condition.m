function D = filter_condition(D, condition, minT)
%FILTER_CONDITION extracts from the data struct D the trials with the condition
%           passed as argument.
%
%           INPUTS:
%           D: data struct in the datahigh format including a field condition
%           condition: a string 'left', 'right', 'error_right', or 'error_left'
%           minT: a filter to exclude trials with insufficient time points
%
%
%Ruben Pinzon

trial = [];
for i = 1 : length(D)
   if strcmp(D(i).condition, condition) && D(i).T>minT
      trial(end+1) = i;
   end    
end

D = D(trial);