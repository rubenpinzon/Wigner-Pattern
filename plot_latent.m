function plot_latent(Xorth, color, name)
%PLOTLATENT is an auxiliary function to plot the latent variables
%           obtained by the GPFA procedure
%
%
%see also branch2, branch2_cleaned.m
%Ruben Pinzon @ 2015

for v = 1 : 3;%size(Xorth,1)
   subplot(1,3,v)
   plot(Xorth(v,:), 'color', color, 'displayname',name), hold on   
   xlim([0 size(Xorth,2)])
   %ylim([-10 10])
   ylabel(sprintf('X_%d',v))
end


