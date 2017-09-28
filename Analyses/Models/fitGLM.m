function [beta,mse,mdl,combinations] = fitGLM(stim,data,varargin)
% data is a vector of length number of trials (T) or a matrix that is T by
% number of time points to be meaned over.
%
% stim is a binary matrix of length number of trials by number of stimuli,
% or a dict of length number of conditions by number of stimuli, with a
% TrialIndex passed in to define when each stimulus was presented.

InteractionTerms = 1; % order of interactions to include
StimIndex = [];
verbose = false;
RunSpeed = [];

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
            case 'RunSpeed'
                RunSpeed = varargin{index+1};
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

x = designStimMatrix(stim,InteractionTerms,RunSpeed);


%% Fit GLM

beta = glmfit(x,data,'normal');
mdl = fitglm(x,data);

% Compute prediction
pred = [~any(x,2),x]*beta; 
mse = mean((pred-data).^2);


%% Display output
if verbose
    
    % Display data values per combination
    temp = [{data(~any(x,2))}, cellfun(@(ind) data(ind), mat2cell(x,numTrials,ones(numCombs-1,1)), 'UniformOutput',false)];
    Labels = cellfun(@num2str,combinations,'UniformOutput',false);
    Labels{1} = 'none';
    figure;
    subplot(2,2,1);
    plotSpread(temp,'showMM',1);
    set(gca,'XTick',1:numCombs,'Xticklabel',Labels,'XTickLabelRotation',90);
    
    % Display weights
    subplot(2,2,2);
    plot(beta);
    set(gca,'XTick',1:numCombs,'Xticklabel',Labels,'XTickLabelRotation',90);
    title('Tuning Curve');
    
    % Display prediction vs actual
    subplot(2,2,[3,4]);
    plot(data); hold on; plot(pred);
    legend('Data','Prediction');
    title(sprintf('Prediction (mse = %f)',mse));
    
else
    fprintf('mse = %f\n', mse);
end



