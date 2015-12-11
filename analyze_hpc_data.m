%ANALYZE_HPC_DATA Script to analyze the synthetic database created by Balázs Ujfalussy to study the GPFA model under
%                 controlled conditions including the number of place cells, firing rate, underlying dimension of
%                 dynamics. The outcome (expected) is a comparison of the estimated latent dimensionality of the database
%                 given by validation with the GPFA, and the real dimension in the data.
%                 the source files are hosted at https://bitbucket.org/bbu20/hpc_data
%                 Description of the database:
%                   An R library to generate synthetic hipocampal data with different intrinsic dimensionality
%                   and environment size.
%
%                 USAGE: Set the variable basepath to the location of the folder HPC-Data containing the database
%                 The function get_matFiles will search and load the spike_Dx_Lx.dat files inside the folder.
%
%                 WHAT IT DOES: Reads the spike files in the database and wrap them up in the DataHigh structure to be
%                 processed by the DataHigh GUI, which carries out dimensionality reduction.
%(This script is still under development)
%
%Version 1.0
%Ruben Pinzon Dec@2015

clc, close all; clear all;

basepath        = '/media/bigdata/';
files           = get_matFiles(basepath,'/spike*','*HPC-Data');
Fs              = 100;
condition       = {'going' 'returning'};
colors          = hsv(3);

for sce = 1 : length(files)
    
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
   %open DimReduce GUI to process the wrapped database
   DataHigh(D, 'DimReduce');
end

