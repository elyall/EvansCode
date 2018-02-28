function [beta,pred,mse,mdl,Filter,combinations,res] = fitGLM(Stim, Data, numLags, varargin)
% Stim is numSamples by numPredictors
% Data is numSamples by N

type = 'fitglm';         % 'glmfit', 'fitglm', or 'fitnlm'
Distribution = 'normal'; % see: glmfit or fitglm
Link = 'identity';       % see: glmfit or fitglm
% Link = @exp_nonlin;      % see: fitnlm
% Link = -1;
beta = [];               % initial coefficients
constant = true;         % booleon specifying whether to include a constant term or not
numNonRegressors = 0;

InteractionTerms = false;% list of order of interactions to include
ExtraVars = [];
verbose = false;
Labels = {};


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'Distribution'
                Distribution = varargin{index+1};
                index = index + 2;
            case 'Link'
                Link = varargin{index+1};
                index = index + 2;
            case 'beta'
                beta = varargin{index+1};
                index = index + 2;
            case 'constant'
                constant = varargin{index+1};
                index = index + 2;
            case 'numNonRegressors'
                numNonRegressors = varargin{index+1};
                index = index + 2;
            case {'InteractionTerms','Terms'}
                InteractionTerms = varargin{index+1};
                index = index + 2;
            case 'ExtraVars'
                ExtraVars = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
            case 'Labels'
                Labels = varargin{index+1};
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

if ~exist('numLags','var') || isempty(numLags)
    numLags = 5;
end


%% Create stim matrix

% Create interaction terms
T = size(Stim,1);
if islogical(InteractionTerms) && isequal(InteractionTerms,true)
    Stim = Stim*Stim';
    combinations = num2cell(1:T);
elseif ~isequal(InteractionTerms,false)
    [Stim,combinations] = designStimMatrix(Stim,InteractionTerms,ExtraVars);
else
    combinations = [];
end
numConds = size(Stim,2) - numNonRegressors;

% Create lagged stimuli
if numLags
    try
        Stim = lagmatrix(Stim,0:numLags); % create lagged stimuli
    catch
        Stim = arrayfun(@(x) cat(1,zeros(x,size(Stim,2)),Stim(1:end-x,:)), 0:numLags, 'UniformOutput', false);
        Stim = cat(2,Stim{:});
    end
end

Stim(isnan(Stim)) = 0;
if constant
    Stim = [ones(T,1),Stim]; % add constant variable (i.e. intercept)
end


%% Initialize outputs
N = size(Data,2);
mdl = cell(1,N);
L = size(Stim,2) - numNonRegressors;
if strcmp(type,'fitnlm')
    L = L+1;
end
if isempty(beta)
    beta = zeros(L,N);
% elseif size(beta,1)~=L
%     error('inital weights must be of size # of predictors (%d currently) by 1 (or by N, %d currently)',L,N);
elseif size(beta,2)~=N
    beta = repmat(beta(:,1),1,N); % use same inital weights for all units
end
pred = nan(T,N);


%% Fit GLM
parfor n = 1:N
    switch type
        
        case 'regress'
            beta(:,n) = regress(Data(:,n),Stim);
            pred(:,n) = Stim*beta(:,n);
            
        case 'glmfit'
            [beta(:,n),Dev,Stats] = glmfit(Stim,Data(:,n),Distribution, 'Link', Link, 'B0', beta(:,n), 'Constant', 'off');
            pred(:,n) = Stim*beta(:,n); % compute prediction
            mdl{n} = {Dev,Stats};
            % pred(:,n) = glmval(beta(:,n),double(Stim),'identity'); % compute prediction (same thing)
            
        case 'fitglm'
            mdl{n} = fitglm(Stim, Data(:,n), 'Link', Link, 'B0', beta(:,n), 'Intercept', false);
            beta(:,n) = mdl{n}.Coefficients{:,1};
            pred(:,n) = mdl{n}.predict;
            
        case 'fitnlm'
            mdl{n} = fitnlm(Stim, Data(:,n)', Link, beta(:,n));
            beta(:,n) = mdl{n}.Coefficients{:,1};
            pred(:,n) = mdl{n}.predict;
            
    end
end
mse = mean((pred-Data).^2); % compute error
res = Data - pred;

switch type
    case 'fitnlm'
        Filter = reshape(beta(constant+1:end-1,:),numConds,numLags+1,N); % reshape weights
    otherwise
        Filter = reshape(beta(constant+1:end,:),numConds,numLags+1,N); % reshape weights
end


%% Display output
if verbose
    for index = 1:N
        if isempty(Labels)
%             Labels = cellfun(@num2str,combinations,'UniformOutput',false);
        end

        figure;
        subplot(2,2,1);
        if numLags == 0 % Display data values per combination
            Stim = logical(Stim);
            temp = [{Data(~any(Stim(:,constant+1:end-numNonRegressors),2),index)}, cellfun(@(ind) Data(ind,index), mat2cell(Stim(:,constant+1:end-numNonRegressors),T,ones(numConds,1)), 'UniformOutput',false)];
            plotSpread(temp,'showMM',1);
            set(gca,'XTick',1:numConds+1,'Xticklabel',[{'none'},Labels],'XTickLabelRotation',90);
            ylabel('Fluorescence');
            title('Data');
        else % Display tuning curve
            plot(mean(Filter(:,:,index),2));
            set(gca,'XTick',1:numConds,'Xticklabel',Labels,'XTickLabelRotation',90);
            xlim([.5,numConds+.5]);
            xlabel('Stimulus');
            ylabel('Weight');
            title('Tuning Curve');
        end
        
        % Display coefficients
        subplot(2,2,2);
        if isvector(Filter(:,:,index))
            plot(Filter(:,:,index));
            set(gca,'XTick',1:numConds,'Xticklabel',Labels,'XTickLabelRotation',90);
            xlim([.5,numConds+.5]);
            xlabel('Stimulus');
            ylabel('Weight');
        else
            imagesc(Filter(:,:,index));
            set(gca,'YTick',1:numConds,'YTickLabel',Labels);
            xlabel('Time Lag (frame)');
            ylabel('Stimulus');
        end
        title('Filter');
        
        % Display prediction vs actual data
        subplot(2,2,[3,4]); plot(Data(:,index)); hold on; plot(pred(:,index));
        legend('Data','Prediction');
        xlabel('Observation');
        ylabel('Signal');
        title(sprintf('Prediction (mse = %f)', mse(index)));
        
    end
% else
%     fprintf('mse = %f\n', mse);
end



