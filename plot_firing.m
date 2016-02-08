function plot_firing(cvdata, true_data, T, varargin)
%PLOTFIRING is an auxiliary function to plot the predicted firing
%           and the real firing rate of the cells during the GPFA
%           training procedure.
%
%
%see also branch2, branch2_cleaned.m
%Ruben Pinzon @ 2015

%[val, order]  = sort(sum((cvdata-true_data).^2,2));
plot_rates = false;
if length(varargin) ~= 0
    rates = varargin{1};
    plot_rates = true;
end
    
order = 1:size(cvdata,1);
for j = 1 : ceil(numel(order)/4)
   figure(j)
   for i = 1 : 4
      idx = (j-1)*4 + i;
      if idx < numel(order)
          subplot(2,2,i)
          y         = true_data(order(idx),:);
          y_pred    = cvdata(order(idx),:);
          
          plot(y,'b'), hold on  
          plot(y_pred,'r', 'linewidth',2)
          
          if plot_rates
              plot(linspace(0,length(y),length(rates)), rates(j,:),'-.k')
          end
          title(sprintf('Cell %d',order(idx)))
          plot(repmat(T,2,1),ylim,'color',[0.8 0.8 0.8])
          
%           subplot(2,4,2*i)
%           histfit(y)
%           set(gca,'view',[90 -90])
      end
   end
end