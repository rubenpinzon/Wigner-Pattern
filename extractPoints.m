function P = extractPoints(W)
%EXTARCTPOINTS is an auxiliar function to extract the start, 1/3, 1/2, 2/3
%and end point of the latent variables in the data struct W, with name
%Xorth. P is [(5 x n_dims) x n_laps], where n_laps = len(W)
%
%
%Ruben Pinzon

n_laps = length(W);
n_dims = size(W(1).Xorth,1);
P      = zeros(5*n_dims,n_laps); 

for k = 1:n_laps
    t      = ceil((W(k).T-1) * [0 1/3 1/2 2/3 1] + 1);                     %points of interest
    points = W(k).Xorth(:,t)';
    P(:,k) = points(:);                                                    %values stacked as in: [x1(1), x2(1),...,x1(end), x2(end)]                                                    
end