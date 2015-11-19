function varargout = readdat(file)

fid    = fopen(file);
header = fgetl(fid);
n_cols = length(strfind(header,' ')) + 1;
fprintf('%d Columns %s\n',n_cols, header)

data   = zeros(1,n_cols);
cnt    = 1;
if ~isempty(str2num(header))
    cnt = 2;
    data(1,:) = str2num(header); 
    disp('Added header to data')
end

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    data(cnt, :) = str2num(tline);
    cnt  = cnt + 1;
end
fclose(fid);
fprintf('file loaded with %d entries\n',length(data))
varargout{1} = data;

if nargout == 2
    varargout{2} =  getlaps(data);
end

end

function laps = getlaps(pos)

%fidn the laps
st_laps        = findpeaks(pos);
en_laps        = findpeaks(-pos);
n_laps         = numel(st_laps)+numel(en_laps)-1;
laps           = zeros(2,n_laps+1);
laps(:,1:2:end)= repmat(st_laps,1,2)';
laps(:,2:2:end)= repmat(en_laps,1,2)';   
laps           = laps(:);
laps([1 end])  = []; 
laps           = reshape(laps,2, n_laps)';

end