%% read spike files:

% SPK file
% --------
% Binary file, 16 bits integers.
% 
% Contains the filtered waveforms, for each spike.
% 
% Array format: electrode/sample/spike


fileID = fopen('i01_maze06_MS.002.spk.5','r');
spk = zeros(36,192);
for i=1:36
    spk(i, :) = fread(fileID,192,'int16') ;
end

fclose(fileID)