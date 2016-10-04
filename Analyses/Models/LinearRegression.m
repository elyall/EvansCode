%% Load data
clear all
load('2481data.mat','Data','Stim','Neuropil','NeuropilWeight');
ROIindex = 220;

% Check data vs. stimuli
temp = bsxfun(@minus, Data(:,ROIindex), min(Data(:,ROIindex))); % subtract minimum to zero
temp = bsxfun(@rdivide, temp, max(temp));                       % normalize by maximum
temp = bsxfun(@plus, temp, numel(ROIindex)-1:-1:0);             % offset each trace
figure; area(numel(ROIindex)*any(Stim,1),'FaceColor',[.9,.9,.9],'EdgeColor',[.9,.9,.9]); hold on; plot(temp); axis tight;


%% Create lagged stimuli
numTimeLags = 100;

% Remove control stimuli
% Stim(1,:) = [];

[numStim,numSamples] = size(Stim);
LagStim = cat(1,Stim,zeros(numStim*(numTimeLags-1), numSamples));
for index = 1:numTimeLags-1
    LagStim(index*numStim+1:(index+1)*numStim,:) = cat(2, zeros(numStim, index), Stim(:,1:end-index));
end
% LagStim = LagStim - mean(LagStim(:)); % subtract off mean

%% Adjust data

% Subtract off neuropil
Data = Data - bsxfun(@times, Neuropil, NeuropilWeight');

% Subtact off mean
% Data = Data - mean(Data(:)); % subtract off mean

%% Compute regression
B = [ones(numSamples,1),LagStim']\Data(:,ROIindex); % compute regression with DC coefficient

%% GLM (start)
% B = glmfit(LagStim',Data(:,ROIindex(1))); % compute GLM

%% Figures

% Display coefficients
filter=reshape(B(2:end,:),numStim,numTimeLags,numel(ROIindex));
figure; imagesc(filter(:,:,1));
xlabel('Time Lag');
ylabel('Stimulus');

% Compute tuning curve
tc = squeeze(mean(filter,2));
figure; plot(tc(:,1));
xlabel('Position');
ylabel('Response (a.u.)');

% Display 1D filters
% figure; plot(filter')

% Compute prediction
pred=[ones(numSamples,1),LagStim']*B;

% Plot overlay
temp2 = bsxfun(@minus, pred, min(Data(:,ROIindex)));             % subtract minimum to zero
temp1 = bsxfun(@minus, Data(:,ROIindex), min(Data(:,ROIindex))); % subtract minimum to zero
temp2 = bsxfun(@rdivide, temp2, max(temp1));                      % normalize by maximum
temp1 = bsxfun(@rdivide, temp1, max(temp1));                     % normalize by maximum
temp2 = bsxfun(@plus, temp2, numel(ROIindex)-1:-1:0);             % offset each trace
temp1 = bsxfun(@plus, temp1, numel(ROIindex)-1:-1:0);             % offset each trace
figure; plot(temp1,'r'); hold on; plot(temp2,'b');
mse = mean((pred-Data(:,ROIindex)).^2);
fprintf('mse = %f\n', mse);
% legend('Data','Prediction');


