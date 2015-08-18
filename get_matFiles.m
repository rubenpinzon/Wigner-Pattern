function files = get_matFiles(animals, basepath)

for ifolder = 1 : length(animals)
    D = dir([basepath animals{ifolder} '/*BehavElectrData.mat']);
    files{ifolder} = [animals{ifolder} '/' D.name];
end