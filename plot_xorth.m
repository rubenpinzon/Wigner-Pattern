function plot_xorth(x1,x2,x3, fig,label, color, varargin)

name = '';
if length(varargin) ~=0
   name = varargin{1}; 
end

subplot(3,3,fig), hold on, grid on %X-Y, 
if isempty(x3)
    
    plot(x1,x2, 'Color', color)
    plot(x1(1),x2(1),'sk', 'markersize',8, 'MarkerFaceColor', color)
    plot(x1(end),x2(end),'ok', 'markersize',8, 'MarkerFaceColor', color)
    xlabel(label{1}), ylabel(label{2}),

    c_axis = axis;
    lim = [min(x1) max(x1) min(x2) max(x2)];

    axis([min(lim(1),c_axis(1)) max(lim(2),c_axis(2)) min(lim(3),c_axis(3)) max(lim(4),c_axis(4))])
else
   subplot(3,3,fig), hold on, grid on
   set(gca, 'CameraPosition',[-10.1288 -7.2519 2.0493])
   plot3(x1,x2,x3, 'Color', color,'displayname',name)
   plot3(x1(1),x2(1),x3(1),'sk', 'markersize',8, 'MarkerFaceColor', color)
   plot3(x1(end),x2(end),x3(end),'ok', 'markersize',8, 'MarkerFaceColor', color)
   xlabel('X_1'), ylabel('X_2'), zlabel('X_3')
   c_axis = axis;
   lim = [min(x1) max(x1) min(x2) max(x2) min(x3) max(x3)];

    axis([min(lim(1),c_axis(1)) max(lim(2),c_axis(2)) min(lim(3),c_axis(3)) max(lim(4),c_axis(4))...
        min(lim(5),c_axis(5)) max(lim(6),c_axis(6))])
    
end