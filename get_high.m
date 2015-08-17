function D=get_high(data, MaxTimeE, trial, color, label, filterTrial)
%function to create the data struct needed for DataHigh
% data is a cell with dims: Laps x Cells
% MaxTimeE is the desired length of the spike trians

tic
[laps, cells] = size(data);
%(1) find the maximum time stamp to create uniform spike vectors
fprintf('Getting %s high: converting data to process with DataHigh', label)
if ~filterTrial
    for l=1:laps
       D(l).data = zeros(cells, MaxTimeE(l));
       D(l).condition = trial{l}; 
       D(l).epochColors = color(l, :);
       for c=1:cells
           spk =  data{l,c}<MaxTimeE(l);
           D(l).data(c, spk) = 1;
       end
       D(l).trialId    = l;
       % Tracking progress
    end
else
    fprintf('...with removal of wrong alterations')
    cnt = 1;
    for l=1:laps
        if strcmp(trial{l}, 'right') || strcmp(trial{l}, 'left')
           D(cnt).data = zeros(cells, MaxTimeE(l));
           D(cnt).condition = trial{l}; 
           D(cnt).epochColors = color(l, :);
           for c=1:cells
               spk =  data{l,c}<MaxTimeE(l);
               D(cnt).data(c, spk) = 1;
           end
           D(cnt).trialId    = l;

            cnt = cnt + 1;
        end
    end
end


fprintf('Done in %3.3f seconds\n', toc)
