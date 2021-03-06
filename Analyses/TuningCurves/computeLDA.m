function [Error,ConfusionMatrices] = computeLDA(ROIdata, varargin)

ROIindex = [1 inf];
StimIndex = [1 inf];
numIter = 10;
% numFolds = 5;

saveOut = false;
saveFile = '';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case {'save','Save'}
                saveOut = true;
                index = index + 1;
            case {'saveFile','SaveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end


%% Load data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = strcat(strtok(ROIFile,'.'),'.lda');
    end
end


%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1),ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIdata.rois);

if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1),StimIndex(end-1)+1:numel(ROIdata.rois(1).Raw)];
end
numStim = numel(StimIndex);


%% Linear-discriminate analysis

% Initialize outputs
Error = nan(numROIs,numIter);
ConfusionMatrices = nan(numROIs,numStim,numStim,numIter);

% Cycle through ROIs computing LDA
fprintf('Computing LDA for %d ROIs (%d iterations)\n',numROIs,numIter);
p = parfor_progress(numROIs);
parfor rindex = 1:numROIs
    
    % Determine minimum number of trials per stimulus
    numTrialsPerStim = cellfun(@numel,ROIdata.rois(ROIindex(rindex)).Raw);
    num = min(numTrialsPerStim(StimIndex));
    
    % Perform multiple iterations, randomly subsampling stimuli with more
    % trials than the minimum
    for iindex = 1:numIter
        
        % Create observation and label vectors
        Data = nan(num*numStim,2);
        for sindex = 1:numStim
            Data((sindex-1)*num+1:sindex*num,1) = randsample(ROIdata.rois(ROIindex(rindex)).Raw{StimIndex(sindex)},num,false);
            Data((sindex-1)*num+1:sindex*num,2) = StimIndex(sindex);
        end
        
        % Perform LDA
        cv = cvpartition(num*numStim,'LeaveOut');                                   % create folds
%         cv = cvpartition(Data(:,2),'KFold',numFolds);                               % create folds
        lda = fitcdiscr(Data(:,1),Data(:,2),'CVPartition',cv);                      % compute LDA
        ldaClass = kfoldPredict(lda); % ldaClass = resubPredict(lda);               % perform prediction
        ConfusionMatrices(rindex,:,:,iindex) = confusionmat(Data(:,2),ldaClass);    % generate confusion matrices
        Error(rindex,iindex) = kfoldLoss(lda); % ldaResubErr = resubLoss(lda);      % calculate error
    end
    
    parfor_progress(p); % update status
end
parfor_progress(p,0);

% Average over iterations
Error = mean(Error,2);
ConfusionMatrices = mean(ConfusionMatrices,4);


%% Save output
if saveOut && ~isempty(saveFile)
    if exist(saveFile,'file')
        save(saveFile,'Error','ConfusionMatrices','ROIindex','-append');
    else
        save(saveFile,'Error','ConfusionMatrices','ROIindex','-mat','-v7.3');
    end
    fprintf('LDA outputs saved to: %s\n',saveFile);
end
    
