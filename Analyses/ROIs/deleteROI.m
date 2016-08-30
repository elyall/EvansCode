function ROIinput = deleteROI(ROIinput, ROIindex, varargin)

saveOut = false;
saveFile = '';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
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

if ~exist('ROIinput', 'var') || isempty(ROIinput)
    [ROIinput,p] = uigetfile({'*.rois;*.segment'}, 'Select ROI file:', directory);
    if isnumeric(ROIinput)
        return
    end
    ROIinput = fullfile(p, ROIinput);
end


%% Load in data
if ischar(ROIinput)
    ROIFile = ROIinput;
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
    [~,~,ext] = fileparts(saveFile);
    switch ext
        case '.segment'
            load(ROIFile,'mask','-mat');
            ROIinput = mask;
        case '.rois'
            load(ROIFile,'ROIdata','-mat');
            ROIinput = ROIdata;
    end
end


%% Remove ROI(s)
if isstruct(ROIinput)
    ROIinput.rois(ROIindex) = [];
elseif isnumeric(ROIinput)
    ROIinput(:,ROIindex) = [];
end


%% Save output to file
if saveOut && ~isempty(saveFile)
    [~,~,ext] = fileparts(saveFile);
    switch ext
        case '.segment'
            if ~isnumeric(ROIinput)
                error('Can''t save to file specified: input needs to be ROImasks');
            end
            mask = ROIinput;
            if ~exist(saveFile, 'file')
                save(saveFile, 'mask', '-mat', '-v7.3');
            else
                save(saveFile, 'mask', '-mat', '-append');
            end
            fprintf('\tROI mask(s) saved to: %s\n',saveFile);
        case '.rois'
            if ~isstruct(ROIinput)
                error('Can''t save to file specified: input needs to be ROIdata');
            end
            ROIdata = ROIinput;     
            if ~exist(saveFile, 'file')
                save(saveFile, 'ROIdata', '-mat', '-v7.3');
            else
                save(saveFile, 'ROIdata', '-mat', '-append');
            end
            fprintf('\tROIdata saved to: %s\n',saveFile);
    end
end

