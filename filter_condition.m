function D = filter_condition(D, condition)


trial = [];
for i = 1 : length(D)
   if strcmp(D(i).condition, condition)
      trial(end+1) = i;
   end    
end

D = D(trial);