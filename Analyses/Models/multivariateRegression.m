function [beta,pred,mse,Filter,combinations] = multivariateRegression(Stim, Data, numLags, varargin)
% Stim is T by numPredictors
% Data is T by N

InteractionTerms = 1; % order of interactions to include, or true for all possible interactions
RunSpeed = [];
verbose = false;
Labels = {};


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'InteractionTerms','Terms'}
                InteractionTerms = varargin{index+1};
                index = index + 2;
            case 'RunSpeed'
                RunSpeed = varargin{index+1};
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


%% Create interaction terms
T = size(Stim,1);
if islogical(InteractionTerms) && isequal(InteractionTerms,true)
    Stim = Stim*Stim';
    combinations = num2cell(1:T);
else
    [Stim,combinations] = designStimMatrix(Stim,InteractionTerms,RunSpeed);
end
numConds = size(Stim,2);


%% Create lagged stimuli
Stim = lagmatrix(Stim,0:numLags); % create lagged stimuli
Stim(isnan(Stim)) = 0;


%% Compute regression

beta = mvregress([ones(T,1),Stim],Data);


%% Compute outputs

% Compute filter for each unit
N = size(Data,2);
Filter = reshape(beta(2:end,:),numConds,numLags+1,N); % shape filter to be numStim x numLags

% Compute prediction and mean-squared error
pred = [ones(T,1),Stim]*beta; % compute prediction
mse = mean((pred-Data).^2);   % compute mean squared error


%% Display output

if verbose
    for index = 1:N
        if isempty(Labels)
            Labels = cellfun(@num2str,combinations,'UniformOutput',false);
        end

        figure;
        subplot(2,2,1);
        if numLags == 0 % Display data values per combination
            temp = [{Data(~any(Stim,2),index)}, cellfun(@(ind) Data(ind,index), mat2cell(Stim,T,ones(numConds,1)), 'UniformOutput',false)];
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
else
    fprintf('mse = %f\n', mse);
end
