function D=segment(D, bin_size, Fs, min_firing)
%SEGMENT remove low firing rate neurons and segments in non-overlapping
%        windows
%
%

firing_thr      = min_firing ; % Minimum firing rate find which 
                        % neurons should be kept
m               = mean([D.data],2) * Fs;
keep_neurons    = m >= firing_thr;
fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
            sum(keep_neurons),firing_thr)
% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(keep_neurons,:);
end
yDim            = sum(keep_neurons);
useSqrt         = 1; % square root tranform for pre-processing?    
                                   
  


bin_width       = ceil(bin_size * Fs); % bin size (Seconds * Fs) = samples

%Extrat bins for one trial, since all the trials
%are of the same duration
for ilap = 1 : length(D)
    seq         = [];
    T           = floor(size(D(ilap).data, 2) / bin_width);
    seq.y       = nan(yDim, T);
    for t = 1:T
      iStart        = bin_width * (t-1) + 1;
      iEnd          = bin_width * t;
      seq.y(:,t)    = sum(D(ilap).data(:, iStart:iEnd), 2);
    end
    %normalization with square root transform
    if useSqrt
        seq.y       = sqrt(seq.y);
    end
    D(ilap).y       = seq.y;
    D(ilap).T       = T;    
end
