function D=get_high(data, MaxTimeE, trial, color, label, filterTrial)
%function to create the data struct needed for DataHigh
% data is a cell with dims: Laps x Cells
% MaxTimeE is the desired length of the spike trians

tic
[laps, cells] = size(data);
%(1) find the maximum time stamp to create uniform spike vectors
fprintf('Getting %s high: converting data to process with DataHigh', label)
if ~filterTrial
    for lapIdx=1:laps
       D(lapIdx).data = zeros(cells, MaxTimeE(lapIdx));
       D(lapIdx).condition = trial{lapIdx}; 
       D(lapIdx).epochColors = color(lapIdx, :);
       for c=1:cells
           spk =  data{lapIdx,c}<MaxTimeE(lapIdx);
           D(lapIdx).data(c, data{lapIdx,c}(spk)) = 1;
       end
       D(lapIdx).trialId    = lapIdx;
       % Tracking progress
    end
else
    fprintf('...with removal of wrong alterations...')
    cnt = 1;
    for lapIdx=1:laps
        if strcmp(trial{lapIdx}, 'right') || strcmp(trial{lapIdx}, 'left')
           D(cnt).data = zeros(cells, MaxTimeE(lapIdx));
           D(cnt).condition = trial{lapIdx}; 
           D(cnt).epochColors = color(lapIdx, :);
           for c=1:cells
               spk =  data{lapIdx,c}<MaxTimeE(lapIdx);

               D(cnt).data(c, data{lapIdx,c}(spk)) = 1;
               %fprintf('max. time neuron %d, %d, lap %d\n', c, max(data{lapIdx,c}(spk)), lapIdx)

           end
           D(cnt).trialId    = lapIdx;
           cnt = cnt + 1;
        end
    end
end


fprintf('Done in %3.3f seconds\n', toc)
