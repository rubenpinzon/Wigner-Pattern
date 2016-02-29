function shadedSd(X, color)

mu  = mean(X,2);
sd  = std(X,1,2);

h1 = area(mu+sd);
h2 = area(mu-sd);
set(h1,'FaceColor',color,'linestyle','none')
set(h2,'FaceColor','w','linestyle','none')
plot(mu, 'color', 'k')