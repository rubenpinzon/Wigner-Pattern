function h=title(fig, txt)

set(fig,'NextPlot','add');
axes;
h = title(txt);
set(gca,'Visible','off');
set(h,'Visible','on'); 