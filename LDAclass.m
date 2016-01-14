function LDAclass(Xtats, label)
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
h1 = gscatter(mod1,mod2,group,'rb','vo',[],'off');
set(h1,'LineWidth',2)
legend(label.modelA,label.modelB,...
       'Location','NW')
   
line(xlim, ylim, 'linestyle','--','color',[0.1 0.1 0.1])
grid on



set(gca,'fontsize',14)
xlabel(label.xaxis,'fontname','Georgia')
ylabel(label.yaxis,'fontname','Georgia')
title(sprintf('%s',label.title))


