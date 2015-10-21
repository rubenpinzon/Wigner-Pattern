%print the latent variables find via DataHigh

[numDim, bins] = size(D(1).data);


figure()
til = sprintf('Latent Variables');
annotation('textbox', [0 0.9 1 0.1], ...
    'String', til, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Fontsize',18)
set(gcf,'color', 'w', 'position', [100 100 1400 700])

time = linspace(0, 3.0, bins); % from he extraction program

for l= 1:length(D)
   for v = 1:numDim
       subplot(3, 3, v)
       plot(time, D(l).data(v, :),'color',D(l).epochColors), hold on       
   end
end