% Script to convert the files in the HC-5/3 database to vectors in the format of Attila's software

% extract spike times and cell numbers in txt file

FileBase = '/media/bigdata/hc-3/ec013.15/ec013.157/ec013.157';

[T,G,Map,Par]=LoadCluRes(FileBase);

save4attila(Map, [FileBase '_mapElec'], '%d\t%d\t%d\n');

%spike times and cell number
save4attila([T/20 G], [FileBase '_spikes'], '%5.2f\t%3d\n');

[data OrigIndex]= LoadBinary([FileBase '.eeg'], 1);
save4attila(data, [FileBase '_eeg'],'%5.2f\n');

whl = load([FileBase '.whl']);
pos = interp1(whl(:,1), linspace(0,length(whl(:,1)), length(data)));
save4attila(pos, [FileBase '_posX'], '%5.2f\n');


%% detect SWRs

plot(data./50), hold on
plot(pos./4 - 40, 'r')
plot(20*imf(7,:)./max(imf(7,:)), 'k')
plot(20*imf(3,:)./max(imf(3,:)),'m')

color = jet(20);
for i = 1:15
    Fs = 1250;
    NFFT = 1024; % Next power of 2 from length of y
    Y = fft(imf(i,:),NFFT);
    f = Fs/2*linspace(0,1,NFFT/2+1);

    % Plot single-sided amplitude spectrum.
    plot(f,2*abs(Y(1:NFFT/2+1)),'displayname',sprintf('imf%d',i), 'color', color(i,:)) , hold on
end
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')