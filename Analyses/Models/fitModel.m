%% 2481-SBX analysis
clear all

id = 2481;
day = 150721;

Files = {'/home/elyall/Documents/Data/2481/150721/2481_179_000',...
    '/home/elyall/Documents/Data/2481/150721/2481_179_001'};
PWCZ = 4:7;

% mountdir = '/run/user/1000/gvfs';


%% Load data

ROIindex = repmat((1:718)', 1, 2);
[numROIs, numFiles] = size(ROIindex);
FileIndex = repmat(1:numFiles, numROIs, 1);
[~, ~, ~, ~, ~, ROIs] = gatherROIdata(strcat(Files,'.rois'), 'curve', ':', 'none', ROIindex(:), FileIndex(:));

load([Files{1},'.exp'], 'AnalysisInfo', '-mat');


%% Save data

% Data: numSamples x numROIs
Data = reshape([ROIs{1}.rois(:).rawdata], length(ROIs{1}.rois(1).rawdata), numel(ROIs{1}.rois));
Neuropil = reshape([ROIs{1}.rois(:).rawneuropil], length(ROIs{1}.rois(1).rawdata), numel(ROIs{1}.rois));
NeuropilWeight = determineNeuropilWeight(ROIs{1});

% Stim matrix: numStim x numSamples
[StimID, ~, StimIndex] = unique(AnalysisInfo.StimID);
Stim = zeros(numel(StimID), size(Data,1));
for tindex = 1:length(StimIndex)
    Stim(StimIndex(tindex),AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2)) = 1;
end

% Remove nan signal
Data([1,end],:) = [];
Neuropil([1,end],:) = [];
Stim(:,[1,end]) = [];

save('2481data.mat','Data','Stim','Neuropil','NeuropilWeight','-v7.3');


%% Format data (all frames)
ROIid = 220;

% Data vector: 1 x numSamples
numSamples = size(Data,1);

% Subtract off neuropil
Data = Data - NeuropilWeight*ROIs{1}.rois(ROIid).rawneuropil';



% Remove control trials
% if any(StimID==0)
%     Stim(StimID==0,:) = [];
%     numStim = numStim - 1;
% end


%% Format data (pulled from trial parsed data)

% Convert ROIs to matrix
Data = ROIs{1}.rois(220).data; % select ROI
numFramesPerTrial = size(Data, 2);

% Determine stimuli
[StimID, ~, StimIndex] = unique(AnalysisInfo.StimID);
numStim = numel(StimID);
numTrials = length(StimIndex);

% Stim matrix: numStim x numSamples
Stim = zeros(size(Data,1), size(Data,2), numStim);
start = ROIs{1}.DataInfo.numFramesBefore+1;
stop = ROIs{1}.DataInfo.numFramesBefore+ROIs{1}.DataInfo.numStimFrames;
for tindex = 1:numTrials
    Stim(tindex,start:stop(tindex),StimIndex(tindex)) = 1;
end
Stim = permute(Stim, [3, 2, 1]);
Stim = reshape(Stim, numStim, numTrials*numFramesPerTrial);

% Data vector: 1 x numSamples
Data = Data';
Data = Data(:);
numSamples = length(Data);


%% Create lagged stimuli

% Check data vs. stimuli
figure; area(any(Stim,1),'FaceColor',[.9,.9,.9],'EdgeColor',[.9,.9,.9]); hold on; plot(Data./max(Data));

numTimeLags = 100;

LagStim = cat(1,Stim,zeros(numStim*numTimeLags, numSamples));
for index = 1:numTimeLags
    LagStim(index*numStim+1:(index+1)*numStim,:) = cat(2, zeros(numStim, index), Stim(:,1:end-index));
end


%% Compute regression

B=[LagStim',ones(numSamples,1)]\Data; % compute regression with DC coefficient

% m = mean(Data(:));
% B=[LagStim' ones(numSamples,1)]\Data - m; % compute regression with DC coefficient and after zeroing data

% Display coefficients
filter=reshape(B(1:numTimeLags*numStim),numStim,numTimeLags);
figure; imagesc(filter);
xlabel('Time Lag');
ylabel('Stimulus');

% Compute tuning curve
tc = mean(filter,2);
figure; plot(tc);
xlabel('Position');
ylabel('Response (a.u.)');

% Display 1D filters
% figure; plot(filter')

% Compute prediction
pred=[LagStim' ones(numSamples,1)]*B;
% pred=[LagStim' ones(numSamples,1)]*B + m; % compute prediction and add back mean
figure; plot(Data); hold on; plot(pred);
mse = mean((pred-Data).^2);
fprintf('mse = %f\n', mse);
legend('Data','Prediction');


%% Display normal tuning curve
hA = plotTuningCurve(ROIs{1}.rois(220), 1, 'curveColor', 'k', 'Title', '', 'YLim', 'tight');

