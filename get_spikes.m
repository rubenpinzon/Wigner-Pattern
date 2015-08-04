function [spk_per_neuron, spk_per_lap] = get_spikes(clusters, spikes, laps)

for j=1:max(clusters) %for each neuron  
   spk_per_neuron{j}  = spikes(clusters==j);
   
   for k=1:length(laps)-1
        index = spk_per_neuron{j}>=laps(k) & spk_per_neuron{j}<laps(k+1);
        spk_per_lap{k,j} = spk_per_neuron{j}(index);
   end
end