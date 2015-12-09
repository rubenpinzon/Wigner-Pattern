function [P,idx] = shuffcells(P)
%Function SHUFFCELLS randomly permutes the cells in the variable data of
%the input stucture.

n_cells   = size(P(1).y,1);
idx       = randperm(n_cells); 
for lap = 1 : length(P)
   P(lap).y = P(lap).y(idx,:); 
end
    
