function status = copyMatFile(Filename,SaveFile)

if ~exist('SaveFile','var') || isempty(SaveFile)
    SaveFile = Filename;
end


%% Load in file
load(Filename,'-mat');
Vars = whos(matfile(Filename));
Vars = {Vars(:).name};


%% Save to file
if exist(SaveFile,'file') % save to temporary file first, then rename it and overwrite original
    delimeter = strfind(SaveFile,'.');
    delimeter = delimeter(end);
    tempFile = [SaveFile(1:delimeter-1),'.temp',SaveFile(delimeter:end)];
    save(tempFile, Vars{:}, '-mat', '-v7.3');
    delete(SaveFile);
    status = movefile(tempFile, SaveFile);
else
    save(saveFile, Vars{:}, '-mat', '-v7.3');
end

