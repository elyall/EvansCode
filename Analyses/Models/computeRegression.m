function [Filter,pred] = computeRegression(Data, Stim, numLags)
% Data is 1xnumSamples
% Stim is numStimxnumSamples

Labels = {};

%% Create lagged stimuli

[numCond,numSamples] = size(Stim);
if isempty(Labels)
    Labels = cellstr(num2str((1:numCond)'));
end
LagStim = cat(1, Stim, zeros(numCond*numLags,numSamples));
for index = 1:numLags
    LagStim(index*numCond+1:(index+1)*numCond,:) = cat(2, zeros(numCond,index), Stim(:,1:end-index));
end


%% Compute regression

B=[LagStim',ones(numSamples,1)]\Data; % compute regression with DC coefficient
% m = mean(Data(:));                        % compute mean
% B=[LagStim',ones(numSamples,1)]\(Data-m); % compute regression with DC coefficient on zeroed data


%% Compute outputs

Filter=reshape(B(1:numLags*numCond),numCond,numLags);   % shape filter to be numStim x numLags
pred=[LagStim',ones(numSamples,1)]*B;                   % compute prediction
% pred=[LagStim',ones(numSamples,1)]*B + m;               % compute prediction and add back mean


%% Plot output

% Display coefficients
figure; imagesc(Filter);
set(gca,'YTick',1:numCond,'YTickLabel',Labels);
xlabel('Time Lag (frame)');
ylabel('Stimulus');

% Compute tuning curve
tc = mean(Filter,2);
figure; plot(tc);
set(gca,'XTick',1:numCond,'XTickLabel',Labels);
xlabel('Stimulus');
ylabel('Response (a.u.)');

% Display 1D filters
% figure; plot(filter')

% Display prediction vs actual data
figure; plot(Data); hold on; plot(pred);
mse = mean((pred-Data).^2);
fprintf('mse = %f\n', mse);
legend('Data','Prediction');
xlabel('Time (frame)')
ylabel('Fluorescence')

