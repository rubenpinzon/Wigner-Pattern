function D = filter_condition(D, condition, minT)


trial = [];
for i = 1 : length(D)
   if strcmp(D(i).condition, condition) && D(i).T>minT
      trial(end+1) = i;
   end    
end

D = D(trial);