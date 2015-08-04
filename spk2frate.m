function frate = spk2frate(spk, window, Fs, basepath)

sd      = window * Fs; %in samples
gk      = exp(-((-6*sd:6*sd)./sd).^2/2); % Gaussian kernel.

if iscell(spk)
    frate = cell(size(spk));
end

if exist([basepath 'firingrate.mat'], 'file') ~= 2
    for s=1:size(spk,2)
        if ~isempty(spk{s})
            st                        = zeros(1,spk{s}(end)+1); %aligned with time=0
            st(spk{s}+1)              = 1;
            frate{s}                  = conv(st,gk); 
        end

        % Tracking progress
        percent = s/size(spk,2)*100;
        w       = 50; % Width of progress bar
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1,...
                round(percent*w/100)), '>',...
                repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
    save([basepath 'firingrate.mat'],'frate');
else
    disp('Loading existing Firing Rate file')
    frate = load([basepath 'firingrate']);
    frate = frate.frate;
end