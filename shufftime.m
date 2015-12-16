function P = shufftime(P)
%SHUFFTIME auxiliary function to randomly permute the time samples/bins in the variable data of
%           the input structure P.
%
%Ruben Pinzon 2015


for lap = 1 : length(P)
    n_bins    = P(lap).T;
    idx       = randperm(n_bins); 
    
    P(lap).y = P(lap).y(:,idx); 
    
    P(lap).shuffled_bins = idx;
end
    
