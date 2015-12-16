function save4attila(data, name, str_format)
%SAVE4ATTILA is a function to save a matlab variable into a file compatible with Attila's software format
%
%       INPUTS:
%       data: matlab vector, matrix
%       name: name of the file to save including extension.
%       str_format: the format to save the data into the file following fprintf.
%
%see savevectors_hc-3.m for usage examples.
%
%Ruben Pinzon@2015

[f c]=size(data);
if c>f; data=data'; c=f; end;

if ~isempty(str_format) %text file tab separated
     fileID = fopen([name '.txt'],'w');
    if c==1
        fprintf(fileID,str_format,data);
    else
        disp('printing by line')
        for i=1:f
            fprintf(fileID,str_format,data(i,:));
        end        
    end
    fclose(fileID);
    type = 'Text';
else
    fileID = fopen([name '.bin'],'w');
    fwrite(fileID, data,'real*4');
    fclose(fileID);
    type = 'Binary';
end
disp([type ' file saved saved at ' name])