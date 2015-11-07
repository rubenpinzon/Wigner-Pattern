function plot_latent(Xorth, color, name)


for v = 1 : 3;%size(Xorth,1)
   subplot(1,3,v)
   plot(Xorth(v,:), 'color', color, 'displayname',name), hold on   
   xlim([0 size(Xorth,2)])
   %ylim([-10 10])
   ylabel(sprintf('X_%d',v))
end


