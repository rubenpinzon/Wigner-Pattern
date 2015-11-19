%script to study the synthetic database created by Balasz
%
%https://bitbucket.org/bbu20/hpc_data
%
%DATA:
% the position of the animal and the spike-times are provided in 25
% simulations
% 	- the length of the data is 10 min with 100 Hz (60 000 data points) 
% 	- the activity of 100 neurons is provided
% 	- only cells with at least 1 place field in the given environment are
%     simulated
% 
% the size of the environment, L={4, 8, 6, 32, 64} and its dimensionality,
% D={1-5} are varied
% 
% pos_Dx_Ly.dat: each row contains the n-dimensional position, time 
% goes from 0.01 s to 600 s
% spikes_Dx_Ly.dat: first column is spike time (s) second is cell-id.
%

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
files           = get_matFiles(basepath,'/spike*','*HPC-Data');
Fs              = 100;
condition       = {'going' 'returning'};
colors          = hsv(3);

for sce = 2 : 2%length(files)
    
   data           = readdat(files{sce});  
   [pos, laps]    = readdat(strrep(files{sce},'spike','pos'));

    
   d_laps         = (laps(:,2)-laps(:,1));
   n_laps         = length(laps);
   [spk, spk_lap] = get_spikes(data(:,2),data(:,1),[laps; 0, 0]./Fs);   
   max_lap        = min(d_laps);
   
   for ilap = 1 : n_laps
%        direction lap
        dir_lap = 2; 
       if pos(laps(ilap,1)) < pos(laps(ilap,2))
          dir_lap = 1; 
       end
       y     = zeros(length(spk), max_lap);
       for cell = 1 : length(spk)
          s          = ceil(100 * spk_lap{ilap,cell}) - laps(ilap,1) + 1; %conver to samples
          s(s>max_lap) = [];
          y(cell, s) = 1;    
       end
       D(ilap).data    = y;
       D(ilap).name    = files{sce};
       D(ilap).spk     = spk_lap;
       D(ilap).pos     = pos;
       D(ilap).trialId = ilap;
       D(ilap).T       = max_lap;
       D(ilap).condition   = condition{dir_lap};
       D(ilap).epochColors = colors(dir_lap,:);
   end
   %DataHigh(D, 'DimReduce');

   
end

