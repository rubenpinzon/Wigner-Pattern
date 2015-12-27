function [F,keep_neurons]=segment(D, bin_size, Fs, keep_neurons, name_var, maxTime)
%SEGMENT remove low firing rate neurons and segments in non-overlapping
%        windows
%        namevar: is the name of the field in the structure D which is to be
%        segmented.
%        Fs: sampling frequency
%        keep_neurons: a vector indicating neurons to include '1' and exclude '0'
%        bin_size: size of the segmentation bin
%        maxTime: maximum segmenting time s 
%
%Ruben Pinzon@2015

if length(keep_neurons)==1
    min_firing      = keep_neurons;
    firing_thr      = min_firing; % Minimum firing rate find which  neurons should be kept
    m               = mean(eval(['D.' name_var]),2) * Fs;
    keep_neurons    = m >= firing_thr;          
else
    disp('Vector of neurons to remove provided')
    firing_thr      = NaN;
end

fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
                sum(keep_neurons),firing_thr)
          


% Remove low firing rate neurons
for itrial = 1:length(D)
    Temp(itrial).data = eval(['D(itrial).' name_var '(keep_neurons,:);']);
end
yDim            = sum(keep_neurons);
useSqrt         = 1; % square root tranform for pre-processing?    
                                   

bin_width       = ceil(bin_size * Fs); % bin size (Seconds * Fs) = samples

%Extrat bins for one trial, since all the trials
%are of the same duration
for ilap = 1 : length(Temp)
    seq         = [];
    T           = floor(size(Temp(ilap).data, 2) / bin_width);
    if maxTime ~= 0
       T_requested = floor(maxTime * Fs /bin_width); 
       if T_requested > T
           disp('ERROR: Requested time larger than lenght of trial')
           return
       end
       T = T_requested;
    end
    
    seq.y       = nan(yDim, T);
    for t = 1:T
      iStart        = bin_width * (t-1) + 1;
      iEnd          = bin_width * t;
      seq.y(:,t)    = sum(Temp(ilap).data(:, iStart:iEnd), 2)./bin_size;
    end
    %normalization with square root transform
    if useSqrt
        seq.y       = sqrt(seq.y);
    end
    
    F(ilap).trialId = D(ilap).trialId;
    F(ilap).type    = D(ilap).type;
    F(ilap).y       = seq.y;        
    F(ilap).T       = T;    
    
end
