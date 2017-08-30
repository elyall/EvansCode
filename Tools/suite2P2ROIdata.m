function ROIdata = suite2P2ROIdata(Ffile, varargin)

ImageFile = '';
Depth = [];
ROIindex = [1 inf];
saveOut = false;
saveFile = '';

directory = cd;

%% Process input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ImageFile'
                ImageFile = varargin{index+1};
                index = index + 2;
            case 'Depth'
                Depth = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Ffile','var') || isempty(Ffile)
    [Ffile,p] = uigetfile({'*.mat'},'Choose F file',directory);
    if isnumeric(Ffile)
        return
    end
    Ffile = fullfile(p,Ffile);
end


%% Load in data
load(Ffile,'stat','ops','Fcell','FcellNeu','-mat');

if isempty(ImageFile)
    ImageFile = ops.fsroot{1}.name;
end
if isempty(Depth)
    Depth = str2double(Ffile(end));
end
if ROIindex(end)==inf
    ROIindex = [ROIindex(1:end-1),ROIindex(end-1)+1:numel(stat)];
elseif ischar(ROIindex) || isequal(ROIindex,true) % load in ROI to save from GUI/classifier created '_proc' file
    if isequal(ROIindex,true)
        ROIindex = [Ffile(1:end-4),'_proc.mat'];
    end
    load(ROIindex,'dat','-mat')
    ROIindex = find([dat.stat(:).iscell]);
end

%% Create ROIdata
numROIs = numel(ROIindex);

% Pull out masks
ROIMasks = false(ops.Ly,ops.Lx,numROIs);
for rindex = 1:numROIs
    ROIMasks(stat(ROIindex(rindex)).ypix,stat(ROIindex(rindex)).xpix,rindex) = true;
end

% Save data
ROIdata = createROIdata(ROIMasks, 'ImageFile', ImageFile, 'Depth', Depth, 'data', Fcell{1}(ROIindex,:), 'neuropil', FcellNeu{1}(ROIindex,:));

%% Save output
if saveOut
    if isempty(saveFile)
        saveFile = [Ffile(1:end-3),'rois'];
    end
    if ~exist(saveFile, 'file')
        save(saveFile,'ROIdata','-v7.3');
    else
        save(saveFile,'ROIdata','-append');
    end
end
