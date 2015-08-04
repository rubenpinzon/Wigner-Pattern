function [mfrate, gfrate] =meanfrate(spk, window, tmax)

bin_min = 0;
bin_max = tmax;
edges = bin_min:window:bin_max;
N=size(spk,2);

freq = zeros(1,length(edges));
for i=1:N
    freq = freq + histc(spk{i},edges)';
end
mfrate = freq./N*window;
