function D = filterbin(D, firing_thr, bin_width)
% FILTERBIN filter neurons based on the minimum firing rate and bin data
%
%
%

Laps            = length(D);
m               = mean([D.data],2) * 1250;
keep_neurons    = m >= firing_thr;
fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
            sum(keep_neurons),firing_thr)
% Remove low firing rate neurons
for itrial = 1:Laps
    D(itrial).data = D(itrial).data(keep_neurons,:);
end
yDim            = sum(keep_neurons);
 
                                   
%Extrat bins for each trial one by one
for ilap = 1 : Laps
    seq         = [];
    T           = floor(size(D(ilap).data, 2) / bin_width);
    seq         = nan(yDim, T);
    for t = 1:T
      iStart        = bin_width * (t-1) + 1;
      iEnd          = bin_width * t;
      seq  (:,t)    = sum(D(ilap).data(:, iStart:iEnd), 2);
    end
    
    D(ilap).bins    = seq;
    D(ilap).T       = T;    
end