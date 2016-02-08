function LDAclass(Xtats, label, color)
% Perform LDA classification of input data which has a struct with fields
%           conf_matrix, class_output, real_label, likelihood: [2xN double] 
%           where N is the number of samples.
%
%Ruben Pinzon@2015
twoModels = size(Xtats.likelihood,1)>1;   


mod1 = Xtats.likelihood(1,:);
mod2 = Xtats.likelihood(2,:);

group       = Xtats.real_label';

figure,hold on
set(gcf, 'color','w')
plot(mod1(group==1), mod2(group==1), 'o', 'color', color(1,:),'linewidth',2, 'markersize',9)
plot(mod1(group==2), mod2(group==2), 'v', 'color', color(2,:),'linewidth',2, 'markersize',9)

legend(label.modelA,label.modelB,...
       'Location','NW')
   
line([-15000 15000], [-15000 15000], 'linestyle','--','color',[0.1 0.1 0.1])
grid on



set(gca,'fontsize',14)
xlabel(label.xaxis,'fontname','Georgia')
ylabel(label.yaxis,'fontname','Georgia')
title(sprintf('%s',label.title))
xlim(minmax(mod1))
minmax(mod1)
ylim(minmax(mod2))
minmax(mod2)
set(gca,'dataaspectratio',[1 1 1])
