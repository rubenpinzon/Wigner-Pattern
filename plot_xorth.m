function plot_xorth(x1,x2,x3, fig,label, color, varargin)
%PLOT_XORTH auxiliary function to show the neural space constructed with the three main vectors in the
%           SVD of the mapping matrix C in the GPFA model.
%
%       INPUTS:
%       x1, x2, x3 : three main vectors of SVD(C), where C is the mapping matrix in the GPFA
%       fig : handles of the parent figure
%       label : name of the projection
%       color : linecolor property
%
%see also branch2, branch2_cleaned.m
%Ruben Pinzon @ 2015


name = '';
if length(varargin) ~=0
   name = varargin{1}; 
end

subplot(3,3,fig), hold on, grid on %X-Y, 
if isempty(x3)
    
    %plot(x1,x2, 'Color', lighthen(color,0.5))
    plot(x1(1),x2(1),'pk', 'markersize',14, 'MarkerFaceColor', color,'linewidth',1,'MarkerEdgeColor', 'none')
    plot(x1(end),x2(end),'ok', 'markersize',8, 'MarkerFaceColor', color,'linewidth',1,'MarkerEdgeColor', 'none')
    xlabel(label{1}), ylabel(label{2}),
    set(gca,'fontname','Georgia','fontSize',14,'linewidth',1.5)
    c_axis = axis;
    lim = [min(x1) max(x1) min(x2) max(x2)];

    axis([min(lim(1),c_axis(1)) max(lim(2),c_axis(2)) min(lim(3),c_axis(3)) max(lim(4),c_axis(4))])
else
   subplot(3,3,fig), hold on, grid on
   set(gca, 'CameraPosition',[-10.1288 -7.2519 2.0493],'fontname','Georgia','fontSize',14,'linewidth',1.5)
   plot3(x1,x2,x3, 'Color', lighthen(color,0),'displayname',name,'linewidth',1,'MarkerEdgeColor', 'none')
   plot3(x1(1),x2(1),x3(1),'pk', 'markersize',14, 'MarkerFaceColor', color,'MarkerEdgeColor', 'none')
   plot3(x1(end),x2(end),x3(end),'ok', 'markersize',8, 'MarkerFaceColor', color,'MarkerEdgeColor', 'none')
   xlabel('{\itx}_0'), ylabel('{\itx}_1'), zlabel('{\itx}_2')
   c_axis = axis;
   lim = [min(x1) max(x1) min(x2) max(x2) min(x3) max(x3)];

    axis([min(lim(1),c_axis(1)) max(lim(2),c_axis(2)) min(lim(3),c_axis(3)) max(lim(4),c_axis(4))...
        min(lim(5),c_axis(5)) max(lim(6),c_axis(6))])
    
end

end

function c = lighthen(c, per)
    c = c + max(c)*per*[1 1 1];
    c(c>1) = 1;
end