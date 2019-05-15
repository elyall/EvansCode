function newFiles = fixFilenames(Files,PTR,Ext)

Files = findFiles(Files,PTR); % pull files out from all input directories
numFiles = numel(Files);

if ~exist('Ext','var')
    Ext = [];
end

% Number after '_' (e.g. sdfasd_13.dfad)
first = cellfun(@(x) find(x=='_',1,'last'), Files)';
last = cellfun(@(x) find(x=='.',1,'last'), Files)';
id = arrayfun(@(x,y,z) str2double(Files{x}(y:z)), 1:numFiles, first+1, last-1);
if isempty(Ext) % keep ending
    newFiles = arrayfun(@(x) sprintf('%s%04d%s', Files{x}(1:first(x)), id(x), Files{x}(last(x):end)) , 1:numFiles, 'UniformOutput',false);
else
    newFiles = arrayfun(@(x) sprintf('%s%04d%s', Files{x}(1:first(x)), id(x), Ext) , 1:numFiles, 'UniformOutput',false);
end

% % Numbers between two '_' (e.g. sdfasd_13_sdf.dfad)
% inds = cellfun(@(x) find(x=='_',2,'last'), Files,'UniformOutput',false);
% inds = cat(1,inds{:})';
% last = cellfun(@(x) find(x=='.',1,'last'), Files)';
% id = arrayfun(@(x,y,z) str2double(Files{x}(y:z)), 1:numFiles, inds(1,:)+1, inds(2,:)-1);
% newFiles = arrayfun(@(x) sprintf('%s%04d.%s', Files{x}(1:inds(1,x)), id(x)) , 1:numFiles, 'UniformOutput',false);

for f = 1:numFiles
    movefile(Files{f},newFiles{f});
end
