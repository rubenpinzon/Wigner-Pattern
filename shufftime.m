function P = shufftime(P)
%Function SHUFFTIME randomly permutes the bins in the variable y of
%the input stucture.


for lap = 1 : length(P)
    n_bins    = P(lap).T;
    idx       = randperm(n_bins); 
    
    P(lap).y = P(lap).y(:,idx); 
    
    P(lap).shuffled_bins = idx;
end
    
