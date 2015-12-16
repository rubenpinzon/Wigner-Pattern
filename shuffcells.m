function [P,idx] = shuffcells(P)
%SHUFFCELLS auxiliary function to randomly permute the cell ID in the variable data of
%           the input structure P.
%
%Ruben Pinzon 2015

n_cells   = size(P(1).y,1);
idx       = randperm(n_cells); 
for lap = 1 : length(P)
   P(lap).y = P(lap).y(idx,:); 
end
    
