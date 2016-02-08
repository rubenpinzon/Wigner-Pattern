function raster(X)
%RASTER a simple function to plot a raster plot of spike trains in X.
%
%
%Ruben Pinzon

t_max = 0;
color = jet(length(X));
for t = 1 : length(X)
    x     = X{t};
    n_evt = length(x);
    if n_evt~=0
        t_max = max(t_max,max(x));
    end
    for n = 1:n_evt
        line(x(n)*[1 1],[t t+0.8],'color',color(t,:),'linewidth',2) 
    end
end
ylabel('Cell No.')
xlim([0 t_max])