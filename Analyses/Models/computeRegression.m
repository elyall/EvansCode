function [Filter,pred] = computeRegression(Data, Stim, numLags, varargin)
% Data is N x numSamples
% Stim is numStim x numSamples

Labels = {};
interactions = false;
verbose = true;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Labels'
                Labels = varargin{index+1};
                index = index + 2;
            case {'Interactions','interactions'}
                interactions = varargin{index+1};
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

if ~exist('numLags','var') || isempty(numLags)
    numLags = 5;
end


%% Create lagged stimuli

[numCond,numSamples] = size(Stim);
if isempty(Labels)
    Labels = cellstr(num2str((1:numCond)'));
end
LagStim = cat(1, Stim, zeros(numCond*numLags,numSamples));
for index = 1:numLags
    LagStim(index*numCond+1:(index+1)*numCond,:) = cat(2, zeros(numCond,index), Stim(:,1:end-index));
end

if interactions
    LagStim = LagStim'*LagStim;
end

%% Compute regression

B=[LagStim',ones(numSamples,1)]\Data'; % compute regression with DC coefficient
% m = mean(Data(:));                        % compute mean
% B=[LagStim',ones(numSamples,1)]\(Data-m); % compute regression with DC coefficient on zeroed data


%% Compute outputs

% Compute filter for each unit
N = size(Data,1);
Filter=reshape(B(1:numLags*numCond,:),numCond,numLags,N);   % shape filter to be numStim x numLags
tc = squeeze(mean(Filter,2));                               % avg over time to create tuning curves

% Compute prediction and mean-squared error
pred=[LagStim',ones(numSamples,1)]*B;                   % compute prediction
% pred=[LagStim',ones(numSamples,1)]*B + m;             % compute prediction and add back mean
mse = mean((pred-Data').^2);                            % compute mean squared error


%% Plot output

if verbose
    for index = 1:N
        figure;
        
        % Display coefficients
        subplot(2,2,1); imagesc(Filter(:,:,index));
        set(gca,'YTick',1:numCond,'YTickLabel',Labels);
        xlabel('Time Lag (frame)');
        ylabel('Stimulus');
        title('Filter');
        
        % Compute tuning curve
        subplot(2,2,2); plot(tc(:,index));
        set(gca,'XTick',1:numCond,'XTickLabel',Labels);
        xlabel('Stimulus');
        ylabel('Response (a.u.)');
        title('Tuning Curve (avg of filter over time)');
        
        % Display prediction vs actual data
        subplot(2,2,[3,4]); plot(Data(index,:)); hold on; plot(pred(:,index));
        legend('Data','Prediction');
        xlabel('Time (frame)')
        ylabel('Signal')
        title(sprintf('Prediction (mse = %f)', mse(index)));
    end
end
