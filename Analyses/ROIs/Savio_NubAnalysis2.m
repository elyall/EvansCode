function Savio_NubAnalysis(N)

addpath(genpath('/global/home/users/elyall/Code/Matlab/'))
parpool('local', 20);


% User defined inputs
tc = true;
glm = true;
ld = false;
nul = false;
ho = false;

str = {'_oopsi','_catch','_movprc'};  % save name
secsBefore = 1;
secsAfter = 2;
ControlID = 0;  % for computing significance in tuning curve
TrialStart = 1; % first trial to analyze
SaveDir = '/media/elyall/Data3/dFoFnew/';
SaveDir = '/global/scratch/elyall/dFoFnew/';

% set file to analyze
File = {};
% L2/3
for d = 1:4
    File = [File,{['/global/scratch/elyall/7142/170710/7142_220_002_depth',num2str(d)]}];
    File = [File,{['/global/scratch/elyall/6994/170725/6994_210_000_depth',num2str(d)]}];
    File = [File,{['/global/scratch/elyall/7120/170731/7120_250_003_depth',num2str(d)]}];
    File = [File,{['/global/scratch/elyall/7197/170807/7197_160_001_depth',num2str(d)]}];
end
% L4
File = [File,{'/global/scratch/elyall/7734/180112/7734_338_000'}];
File = [File,{'/global/scratch/elyall/7734/180112/7734_308_001'}];
File = [File,{'/global/scratch/elyall/7736/180117/7736_300_000'}];
File = [File,{'/global/scratch/elyall/7736/180117/7736_265_001'}];
File = [File,{'/global/scratch/elyall/7737/180118/7737_291_000'}];
File = [File,{'/global/scratch/elyall/7737/180118/7737_326_001'}];
% L2/3 anesthetized
File = [File,{'/global/scratch/elyall/9445/181015/9445_180_005'}];
for d = 1:4
    File = [File,{['/global/scratch/elyall/9019/181021/9019_165_000_depth',num2str(d)]}];
    File = [File,{['/global/scratch/elyall/9025/181021/9025_180_002_depth',num2str(d)]}];
end

% for d = 1:4
%     File = [File,{['/media/elyall/Data/7142/170710/7142_220_002_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data/6994/170725/6994_210_000_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data/7120/170731/7120_250_003_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data/7197/170807/7197_160_001_depth',num2str(d)]}];
% end
% File = [File,{'/media/elyall/Data/7734/180112/7734_338_000'}]; % 1st depth
% File = [File,{'/media/elyall/Data/7734/180112/7734_308_001'}]; % 2nd depth
% File = [File,{'/media/elyall/Data/7736/180117/7736_300_000'}]; % 1st depth
% File = [File,{'/media/elyall/Data/7736/180117/7736_265_001'}]; % 2nd depth
% File = [File,{'/media/elyall/Data/7737/180118/7737_291_000'}]; % 1st depth
% File = [File,{'/media/elyall/Data/7737/180118/7737_326_001'}]; % 2nd depth
% 
% File = [File,{'/media/elyall/Data2/9445/181015/9445_180_005'}];
% for d = 1:4
%     File = [File,{['/media/elyall/Data2/9019/181021/9019_165_000_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data2/9025/181021/9025_180_002_depth',num2str(d)]}];
% end


load([File{N},'.rois'], 'ROIdata', 'TrialIndex', '-mat');
ROIs = ROIdata;

for R = 3

% Load ROIdata
ROIdata = ROIs;
fprintf('Analyzing %s\n',[File{N},str{R},'.rois']);

% Load Experiment info
BaseName = closestFile(File{N},'.sbx'); % determine closest file
BaseName = BaseName{1}(1:end-4);
load([BaseName,'.exp'], 'Experiment', 'TrialInfo','AnalysisInfo','frames','-mat');
StimID  = Experiment.StimID;    % IDs of stimuli presented
StimLog = Experiment.stim.stim; % table of pistons presented in each stimulus
StimIndex = TrialInfo.StimID;   % ID of stimulus on each trial

% Compute trial index
if ~exist('TrialIndex','var') || isempty(TrialIndex)
    try
        TrialIndex = determineRunning([BaseName,'.exp'],'TrialIndex',[TrialStart,inf]);
    catch
        TrialIndex = [1,inf];
    end
	tc = true;
	glm = true;
	ld = true;
end

numFramesBefore = ceil(ROIdata.Config.FrameRate/ROIdata.Config.Depth*secsBefore); %trialwise
numFramesAfter = ceil(ROIdata.Config.FrameRate/ROIdata.Config.Depth*secsAfter); % trialwise


% Set save names
[~,fn,~] = fileparts(File{N});
ROIFile  = fullfile(SaveDir,[fn,str{R},'.rois']);
DataFile = fullfile(SaveDir,[fn,str{R},'.data']);
NullFile = fullfile(SaveDir,[fn,str{R},'_null.data']);

if tc
% Compute tuning curves
NeuropilWeight = determineNeuropilWeight(ROIdata); % determine NeuropilWeight
Data = arrayfun(@(x) ROIdata.rois(x).rawdata - NeuropilWeight(x)*ROIdata.rois(x).rawneuropil,1:numel(ROIdata.rois),'UniformOutput',false);
Data = cat(1,Data{:})';
depthID = idDepth(ROIdata.Config,[],'Depth',ROIdata.depth);    % pull out frame indices for depth of current ROIdata struct
if R == 1 % oopsi
    Data = estimateSpikeTiming(Data); % estimate spikes
elseif R == 2 % catch
    numFrames = ceil(ROIdata.Config.FrameRate);
    Baseline = computeBaseline(Data,AnalysisInfo,'type','catch','numFrames',numFrames,'depthID',depthID);
    Data = (Data - Baseline)./Baseline; % dF/F
elseif R == 3 % movprc
    numFrames = round(60*ROIdata.Config.FrameRate/ROIdata.Config.Depth); % compute over 1 min of frames
    numFrames = 2*floor(numFrames/2)+1; % ensure odd
    Baseline = computeBaseline(Data,AnalysisInfo,'type','movprctile','numFrames',numFrames);
    Data = (Data - Baseline)./Baseline; % dF/F
end
[Data,numFramesBefore,numStimFrames] = trialOrganize(Data, AnalysisInfo, depthID,'numFramesBefore',numFramesBefore,'numFramesAfter',numFramesAfter); % organize trials to numTrials by numFrames
ROIdata = distributeROIdata(ROIdata,'dFoF',squeeze(mat2cell(permute(Data,[2,1,3]),size(Data,2),size(Data,1),ones(size(Data,3),1))));
if R == 1
    [Data,Baseline] = computeEvokedSpikes(Data,1:numFramesBefore); % compute evoked response
end
FrameIndex = arrayfun(@(x) numFramesBefore+1:numFramesBefore+numStimFrames(x), 1:numel(numStimFrames), 'UniformOutput', false);
Data = computeTrialMean2(Data,FrameIndex);   % compute stim mean
Data = mat2cell(Data,ones(size(Data,1),1),size(Data,2));
Data = cellfun(@transpose, Data, 'UniformOutput',false);
ROIdata = distributeROIdata(ROIdata,'stimMean',Data);
[ROIdata,~,outliers] = computeTuningCurve(ROIdata, [], TrialIndex,...
        'ControlID', ControlID,...
        'StimIDs',   StimID,...
        'numBoot',   0);
save(ROIFile,'NeuropilWeight','TrialIndex','ROIdata','outliers','Baseline','-v7.3');
fprintf('Completed computeTuningCurve: %s\n', ROIFile);
else
load(ROIFile,'ROIdata','TrialIndex','outliers','-mat');
end


Raw = gatherROIdata(ROIdata,'Raw');
StimMean = gatherROIdata(ROIdata,'stimMean');
numROIs = size(Raw,1);

if glm
% Compute whisker weights
PW_GLM = nan(6,numROIs);
mse = nan(numROIs,1); 
parfor r = 1:numROIs
    StimMat = StimLog(StimIndex+1,:);
    temp = TrialIndex;
    temp(ismember(TrialIndex,find(outliers(:,r)))) = []; % remove outlier trials
    StimMat = StimMat(temp,:);
    Data = StimMean(r,temp)';
    mdl = fitglm(StimMat, Data, 'linear', 'Link', 'identity');
    PW_GLM(:,r) = mdl.Coefficients{:,1};
    mse(r) = mean((mdl.predict-Data).^2); % compute error
end
PW_GLM = PW_GLM';
save(ROIFile,'PW_GLM','mse','-append');
fprintf('Completed GLM weights: %s\n', ROIFile);
end

if ld
% Compute linear difference from sum
RawCorr = Raw;
Baseline = cellfun(@nanmean, Raw(:,1));
for r = 1:size(Raw,1)
	RawCorr(r,:) = cellfun(@(x) x-Baseline(r), Raw(r,:), 'UniformOutput', false);
end
[LD,CI] = LinDiff(RawCorr,StimLog);
save(DataFile,'LD','CI','-v7.3');
fprintf('Completed LinDiff: %s\n', DataFile);
end

if nul
% Compute null tuning and LD
[nullCurves, nullSE, nullP, nullLD, nullCI, order] = nullTuningAndLD(StimMean, StimIndex, StimLog,...
	'TrialIndex', TrialIndex,...
	'Save',false,'saveFile',NullFile);
save(NullFile,'nullCurves','nullSE','nullP','nullLD','nullCI','order','-v7.3'); % forces file overwrite vs within function which appends
fprintf('Completed Null: %s\n', NullFile);
end

if ho
% Compute hold out curves
Raw_train = cell(size(Raw));
Raw_test  = cell(size(Raw));
parfor ind = 1:numel(Raw)
	cv = cvpartition(numel(Raw{ind}),'HoldOut',.5);
	Raw_train{ind} = Raw{ind}(~cv.test);
	Raw_test{ ind} = Raw{ind}( cv.test);
end
Curves_train = cellfun(@nanmean, Raw_train);
Curves_test  = cellfun(@nanmean, Raw_test );
SE_train = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), Raw_train); % compute standard error
SE_test  = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), Raw_test ); % compute standard error
save(NullFile,'Raw_train','Raw_test','Curves_train','Curves_test','SE_train','SE_test','-append');
fprintf('Completed Hold Out: %s\n', NullFile);
end

end
