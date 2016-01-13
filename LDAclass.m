function LDAclass(Xtats, label)
% Perform LDA classification of input data which has a struct with fields
%           conf_matrix, class_output, real_label, likelihood: [2xN double] 
%           where N is the number of samples.
%
%Ruben Pinzon@2015
twoModels = size(Xtats.likelihood,1)>1;   


mod1 = Xtats.likelihood(1,:);
mod2 = Xtats.likelihood(2,:);


x_lim = [min(mod1) max(mod1)];
y_lim = [min(mod2) max(mod2)];

%LDA
[gridX,gridY] = meshgrid(linspace(x_lim(1),x_lim(2)),linspace(y_lim(1),y_lim(2)));
gridX = gridX(:);
gridY = gridY(:);


group       = Xtats.real_label';
training    = Xtats.likelihood';
sample      = [gridX gridY];
[C,err,P,logp,coeff] = classify(sample,training,group,'linear');
fprintf('Classification error with LDA =%f\n',err)

figure,hold on
set(gcf, 'color','w')
h1 = gscatter(mod1,mod2,group,'rb','vo',[],'off');
set(h1,'LineWidth',2)
legend(label.modelA,label.modelB,...
       'Location','NW')
   
gscatter(gridX,gridY,C,'rb','.',1,'off');hold on;
K = coeff(1,2).const;
L = coeff(1,2).linear;
% Function to compute K + L*v + 
f = @(x,y) K + [x y]*L;



h2 = ezplot(f,[x_lim y_lim]);

set(h2,'Color','m','LineWidth',2)
set(gca,'fontsize',14)
axis([[x_lim y_lim]])
xlabel(label.xaxis,'fontname','Georgia')
ylabel(label.yaxis,'fontname','Georgia')
title(sprintf('%s, error= %f',label.title,err))


