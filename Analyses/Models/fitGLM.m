function [beta,mdl,combinations] = fitGLM(stim,data,varargin)
% data is a vector of length number of trials (T) or a matrix that is T by
% number of time points to be meaned over.
%
% stim is a binary matrix of length number of trials by number of
% conditions, or a dict of length number of stimuli by number of
% conditions, with a TrialIndex passed in to define when each stimulus was
% presented.

InteractionTerms = 1; % order of interactions to include
StimIndex = [];
verbose = false;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case {'InteractionTerms','Terms'}
                InteractionTerms = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if size(data,2)>1
    data = mean(data,2);
end
if size(stim,1)~=size(data,1) % assume stim dict input
    if any(StimIndex==0)
        StimIndex = StimIndex+1;
    end
    stim = stim(StimIndex,:);
end
[numTrials,numConds] = size(stim);


%% Build stim matrix

% Determine combinations
combinations = {};
for index = 1:numel(InteractionTerms)
    currentCombs = combnk(1:numConds,InteractionTerms(index));
    if currentCombs(1)~=1 % combnk sometimes produces the output upside down from what you'd expect
        currentCombs = flipud(currentCombs);
    end
    combinations = [combinations,mat2cell(currentCombs,ones(size(currentCombs,1),1),InteractionTerms(index))'];
end
combinations = [{[]},combinations]; % add zero term (catch trials)
numCombs = numel(combinations);

% Build stimulus matrix
x = false(numTrials,numCombs-1); % doesn't include the zero term
for index = 1:numCombs-1
    x(:,index) = prod(stim(:,combinations{index+1}),2);
end

if verbose
    temp = [{data(~any(x,2))}, cellfun(@(ind) data(ind), mat2cell(x,numTrials,ones(numCombs-1,1)), 'UniformOutput',false)];
    Labels = cellfun(@num2str,combinations,'UniformOutput',false);
    Labels{1} = 'none';
    figure;
    plotSpread(temp,'showMM',1);
    set(gca,'XTick',1:numCombs,'Xticklabel',Labels,'XTickLabelRotation',90);
end

%% Compute GLM

beta = glmfit(x,data,'normal');
mdl = fitglm(x,data);

if verbose
    figure;
    plot(beta);
    set(gca,'XTick',1:numCombs,'Xticklabel',Labels,'XTickLabelRotation',90);
end

%% Compute prediction

pred = [~any(x,2),x]*beta;
mse = mean((pred-data).^2);
fprintf('mse = %f\n', mse);


