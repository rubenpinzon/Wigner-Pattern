function SpkSeries=get_series(data, MaxTimeE)
%function to create a matrix with dims: Laps x Cells X MaxTimeE
% MaxTimeE is the desired length of the spike series


[laps cells] = size(data);
disp(['Extracting spikes series with window ' num2str(MaxTimeE)])
%(1) find the minimum time stamp to create uniform spike vectors
MinTime     = 10000000*ones(1, laps);
SpkSeries   = zeros(laps, cells, MaxTimeE);
for l=1:laps
   for c=1:cells
       spk =  data{l,c};
       if ~isempty(spk)
            MinTime(l) = min(MinTime(l),min(spk));
       end
   end
   
   for c=1:cells
       spk =  data{l,c}(data{l,c}<MaxTimeE+MinTime(l));
       SpkSeries(l, c, spk-MinTime(l)+1) = 1;       
   end
   
   % Tracking progress
    percent = l/laps*100;
    w       = 50; % Width of progress bar
    perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
    disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1,...
            round(percent*w/100)), '>',...
            repmat(' ', 1, w - round(percent*w/100)), ']']);
end
