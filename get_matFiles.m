function [files, names, roots] = get_matFiles(varargin)


basepath = varargin{1};
if nargin == 1
   pattern = '/*BehavElectrData.mat';
else
   pattern = varargin{2}; 
end
fprintf('Searching files with pattern %s\n', pattern)
animals = dir([basepath 'i01*']);
files = {};
names = {};
roots = {};
for ifolder = 1 : length(animals)
    name = animals(ifolder).name;
    D = dir([basepath name pattern]);
    for ifil = 1 : length(D)
        files{end + 1} = [basepath name '/' D(ifil).name];
        names{end + 1} = name;
        roots{end + 1} = [basepath name '/' name];
    end
    
end 