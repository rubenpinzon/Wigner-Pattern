function P = reversetime(P)
%REVERSETIME auxiliary function to reverse the time samples/bins in the variable data of
%           the input structure P.
%
%Ruben Pinzon 2015


for lap = 1 : length(P)
    n_bins    = P(lap).T;
    idx       = n_bins:-1:1; 
    
    P(lap).y = P(lap).y(:,idx); 
    
    P(lap).reversed_bins = true;
end
    