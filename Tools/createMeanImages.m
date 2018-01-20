function createMeanImages(AlignFile, ZStackFile, N)

if ~isempty(AlignFile)
    load(AlignFile,'m','-mat');
    [p,~,~] = fileparts(AlignFile);
    save2P(fullfile(p,'avg.tif'),m,'CLim',true,'Class','uint16');
end

if exist('ZStackFile','var')
    if ischar(ZStackFile)
        [p,~,~] = fileparts(ZStackFile);
    elseif iscell(ZStackFile)
        [p,~,~] = fileparts(ZStackFile{1});
    end
    avgStack(ZStackFile,N,'save','saveFile',fullfile(p,'zstack.tif'));
end