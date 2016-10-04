function plotActivityVsRunning(ROIdata, RunData, varargin)

ROIindex = [1 inf];
TrialIndex = [1 inf];
StimID = 4;
hA = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'StimID'
                StimID = varargin{index+1};
                index = index + 2;
            case {'axes','Axes'}
                hA = varargin{index+1};
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


%% Load in Data
if ischar(ROIdata)
    load(ROIdata, 'ROIdata' ,'-mat');
end
if ischar(RunData)
    RunData = gatherRunData(RunData,[],'TrialIndex',TrialIndex);
end


%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1),ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1),TrialIndex(end-1)+1:max(ROIdata.DataInfo.TrialIndex)];
end
if isempty(StimID)
    StimID = unique(ROIdata.DataInfo.StimID);
end


%% Plot correlation
TrialIndex = ismember(ROIdata.DataInfo.TrialIndex,TrialIndex)';
TrialIndex = TrialIndex & ismember(ROIdata.DataInfo.StimID,StimID);
Z = nan(nnz(TrialIndex),numel(ROIindex));
for rindex = 1:numROIs
    Z(:,rindex) = ROIdata.rois(ROIindex(rindex)).stimMean(TrialIndex);
end
X = repmat(nanmean(RunData(TrialIndex,:),2),1,numel(ROIindex));
Y = repmat(nanstd(RunData(TrialIndex,:),[],2),1,numel(ROIindex));

if ~isempty(hA)
    axes(hA);
else
    figure;
end
plot(X(:),Z(:),'.');
% plot3(X(:),Y(:),Z(:),'.');


%% Calculate running mean
x = min(X(:)):1:max(X(:));

% F = fit(X(:),Z(:),'smoothingspline');
% y = feval(F, x);

% [F,ErrorEst] = polyfit(X(:),Z(:),3);
% y = polyval(F,x,ErrorEst);

% F = glmfit(X(:), Z(:), 'normal','link','probit');
% y = glmval(F,x,'probit');

x = X(:,1);
[x,o] = sort(x);
y = median(Z(o,:),2);
% y = smooth(y);
s = std(Z(o,:),[],2)./numel(ROIindex);

hold on; 
plot(x, y);
% plot(x,s);
hold off;




% SEM = std(x)/sqrt(length(x));               % Standard Error
% ts = tinv([0.025  0.975],length(x)-1);      % T-Score
% CI = mean(x) + ts*SEM;                      % Confidence Intervals



