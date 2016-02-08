function D = filter_laps(D,varargin)
%Filter the laps in nteh struct D with firing rate count below 1 s.d of
%       the mean of all the laps.
%
%Ruben Pinzon @2015
type = 'spike_cnt';
if nargin > 1
    type = varargin{1};
end

if strcmp(type,'spike_cnt')
    n_trials    = length(D);
    cnt_total   = zeros(1, n_trials);
    for n = 1 : n_trials    
        cnt_total(n) = sum(sum(D(n).y))/D(n).T;        
    end

    %trials with low total firing count
    mu = mean(cnt_total);
    sd = std(cnt_total);

    f_out = cnt_total < (mu-sd); 
    
else 
    T = [D.T];
    f_out = T < (mean(T)-std(T));
end

fprintf('%d trials filter out: [%s]\n', sum(f_out), sprintf('%d ',[D(find(f_out==1)).trialId]))

D = D(~f_out);