function saveFile = saveMean(saveFile,Files)

if ~exist('Files','var') || isempty(Files)
    [f,p] = uigetfile({'*.align'},'Choose align files','MultiSelect','on');
    if isnumeric(f)
        saveFile = false;
        return
    end
    Files = fullfile(p,f);
end

if ischar(Files)
    Files = {Files};
end

M = [];
for f = 1:numel(Files)
    load(Files{f},'m','-mat');
    M = cat(3,M,m);
end

save(saveFile,'M');