function [CoM,p,SigChange] = controlChangeWithin(ROIdata)

TrialIndex = [1 inf];   % set running trials
ROIindex = [1 inf];     % set ROIs to compute indices on

stimorder = [];         % set stim order within tuning curves
ControlID = [];         % set which stim ID is a control trial
PWCZ = [2 inf];         % set what curve indices CoM is computed over
DistBtwn = [];          % set distance between stim positions


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'stimorder'
                stimorder = varargin{index+1};
                index = index + 2;
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
            case 'PWCZ'
                PWCZ = varargin{index+1};
                index = index + 2;
            case 'DistBtwn'
                DistBtwn = varargin{index+1};
                index = index + 2;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end


%% Load in data
if ischar(ROIdata)
    load(ROIdata,'ROIdata','-mat');
end

% Determine ROIs to analyze
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois));
end

% Determine PWCZ
if PWCZ(end) == inf
    PWCZ = PWCZ(1):numel(unique(ROIdata.DataInfo.StimIndex));
end


%% Determine trials for each half
if TrialIndex(end) == inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(end-1)+1:max(ROIdata.DataInfo.TrialIndex));
end
TrialIndex = ismember(ROIdata.DataInfo.TrialIndex', TrialIndex);

% split experiment into first vs last half
TrialIndex = find(TrialIndex);
Trials = {TrialIndex(1:ceil(end/2)),TrialIndex(ceil(end/2)+1:end)};


%% compute two tuning curves
Curves = cell(1,2);
ROIs = cell(1,2);
for index = 1:2
    [ROIs{index},Curves{index}] = computeTuningCurve(ROIdata, ROIindex, Trials{index},...
        'ControlID', ControlID,...
        'StimIDs',  stimorder);
end
Curves = cat(3,Curves{1},Curves{2});


%% Calculate metrics

% Preferred Position
CoM = [computeCenterOfMass(Curves(:,:,1),PWCZ,DistBtwn), computeCenterOfMass(Curves(:,:,2),PWCZ,DistBtwn)];
p = ranksum(CoM(:,1),CoM(:,2));

% TMI per position
SigChange = zeros(3,numel(PWCZ));
for rindex = ROIindex
    for pindex = 1:numel(PWCZ)
        pre = ROIs{FileIndex(rindex,1)}.rois(ROIindex(rindex,1)).Raw{PWCZ(pindex)};
        pos = ROIs{FileIndex(rindex,2)}.rois(ROIindex(rindex,2)).Raw{PWCZ(pindex)};
        [~,h] = ranksum(pre, pos);
        if h && mean(pos)-mean(pre)>0
            SigChange(1,pindex) = SigChange(1,pindex)+1;
        elseif h && mean(pos)-mean(pre)<0
            SigChange(2,pindex) = SigChange(2,pindex)+1;
        else
            SigChange(3,pindex) = SigChange(3,pindex)+1;
        end
    end
end
SigChange = SigChange/numel(ROIindex); % turn into a percentage


