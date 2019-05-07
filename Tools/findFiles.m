function OutFiles = findFiles(Files,PTR)

if isempty(Files)
    Files = uigetdir('Choose directory',pwd);
end
if ischar(Files)
    Files = {Files};
end
if ~exist('PTR','var') || isempty(PTR)
    PTR = false;
end

OutFiles = [];
for f = 1:numel(Files)
    if contains(Files{f},'*')
        temp = dir(Files{f});                      % locate files on path
        temp = fullfile(Files{f},{temp(:).name});  % set names to full path
        OutFiles = [OutFiles,temp];                % add files to list
    elseif isfolder(Files{f})
        if PTR
            temp = dir(fullfile(Files{f},PTR));    % locate files on path
        else
            temp = dir(Files{f});                  % locate files in folder
        end
        temp = fullfile(Files{f},{temp(:).name});  % set names to full path
        OutFiles = [OutFiles,temp];                % add files to list
    else
        OutFiles = [OutFiles,Files(f)];            % add file to list
    end
end
OutFiles = OutFiles';

