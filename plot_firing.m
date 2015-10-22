function plot_firing(cvdata, true_data, T)

%[val, order]  = sort(sum((cvdata-true_data).^2,2));
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
          title(sprintf('Cell %d',order(idx)))
          plot(repmat(T,2,1),ylim,'color',[0.8 0.8 0.8])
          
%           subplot(2,4,2*i)
%           histfit(y)
%           set(gca,'view',[90 -90])
      end
   end
end