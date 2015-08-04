function [D, Dave]=get_high(data, MaxTimeE, trial, color)
%function to create the data struct needed for DataHigh
% data is a cell with dims: Laps x Cells
% MaxTimeE is the desired length of the spike trians


[laps, cells] = size(data);
%(1) find the maximum time stamp to create uniform spike vectors
MaxTime = zeros(1, laps);
MinTime = 10000000*ones(1, laps);
disp('Getting high: converting data to process with DataHigh')
Dave = zeros(cells, MaxTimeE);  
for l=1:laps
   for c=1:cells
       spk =  data{l,c};
       if ~isempty(spk)
            MaxTime(l) = max(MaxTime(l),max(spk));
            MinTime(l) = min(MinTime(l),min(spk));
       end
   end
   D(l).data = zeros(cells, MaxTimeE);
   D(l).condition = trial{l}; 
   D(l).epochColors = color(l, :);
   for c=1:cells
       spk =  data{l,c}(data{l,c}<MaxTimeE+MinTime(l));
       D(l).data(c, spk-MinTime(l)+1) = 1;       
   end
   Dave = Dave + D(l).data;
   % Tracking progress
    percent = l/laps*100;
    w       = 50; % Width of progress bar
    perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
    disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1,...
            round(percent*w/100)), '>',...
            repmat(' ', 1, w - round(percent*w/100)), ']']);
end
Dave = Dave./ cells;


