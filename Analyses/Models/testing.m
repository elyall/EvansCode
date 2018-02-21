% GLM - time

% Gather data
[~,order] = sort(M,'descend');
a = Data(order(1),:);
% a = Data(order(2),:);
a = a - prctile(a,30);

% Create stim matrix
% x = designStimMatrix(StimMat,[1,2]);
x = [zeros(1,5);diff(StimMat,[],1)];
x = x>0;
x = designStimMatrix(x,[1]);
% x = [x,designStimMatrix(StimMat,[1,2])];

for i = 1:size(x,2)
    x(:,i) = conv(x(:,i),filter1,'same');
end

% % Add run speed
% RunSpeed = mean(reshape(frames.RunningSpeed,numDepths(dindex),nFrames),1)';
% RunSpeed(isnan(RunSpeed)) = 0; % set frames acquired pre and post daq running to 0
% x = bsxfun(@times, x, RunSpeed > 100);          % set frames where mouse wasn't running as no stimulus
% x(RunSpeed<=100,:) = nan;
% a(RunSpeed<=100) = nan;
% % x = [x, RunSpeed];                              % add RunSpeed as a regressor
% RunSpeed = RunSpeed > 100;
% % x = [bsxfun(@times, x, ~RunSpeed),bsxfun(@times, x, RunSpeed)]; % separate stim matrix by condition with running and condition without running

% % Add whisker angle
% WhiskerData = squeeze(nanmean(reshape(frames.Whiskers,numDepths(dindex),nFrames,size(frames.Whiskers,2),2),1));
% % WhiskerData = mean(WhiskerData,2); % take mean across all whiskers
% % WhiskerData = reshape(WhiskerData,nFrames,size(WhiskerData,2)*2); % regress on angle and curvature
% WhiskerData = WhiskerData(:,:,1); % regress on angle
% x = [x, WhiskerData];
% a(any(isnan(x),2)) = [];
% x(any(isnan(x),2),:) = [];



% Lag stims
numLags = 100;
n = size(x,2);
try
    x = lagmatrix(x,0:numLags); % create lagged stimuli
catch
    x = arrayfun(@(y) cat(1,zeros(y,size(x,2)),x(1:end-y,:)), 0:numLags, 'UniformOutput', false);
    x = cat(2,x{:});
end
% a(1:10) = [];
% x(1:10,:) = [];
x = [true(size(x,1),1),x]; % add bias


% Regress
% d = a'\[ones(size(x,1),1),x];
% d = a'\x;
d = regress(a',x);

% STRF
% f = reshape(d(2:end),5,11);
f = reshape(d(2:end),n,numLags+1);
figure; imagesc(f);

% Prediction
% p = [ones(size(x,1),1),x]*d;
p = x*d;
% p = x(:,1:2:end)*d(1:2:end);
% p2 = x(:,56:end)*d(56:end);
mse = mean((p-a').^2)
rs = abs(a'-p); % residuals

figure;
% plot(zscore(a));
% hold on;
% plot(zscore(p));
plot(a);
hold on
plot(p);
% plot((frames.Stimulus==1)*30000)
% plot(a,rs,'.'); 

%% 
x_nolag = designStimMatrix(StimMat,[1,2]);
filter1 = f(1,:);
conv(x_nolag(:,1),filter1)

figure; 
hold on
plot(conv(x(:,2),f(1,:)))
plot(conv(x_nolag(:,6),f(6,:)))



%%
beta0 = [zeros(size(d)); 5];
beta0 = [d; 0];


mdl = fitnlm(x, a, @exp_nonlin, beta0);
plot(mdl.predict)

f2 = mdl.Coefficients{:,2};
f2 = reshape(f2(1:end-1),n,numLags+1);
figure; imagesc(f2);

p2 = mdl.predict;
mse2 = sum((p2-a').^2)/var(a)
% y=exp_nonlin(mdl,x);

legend({'actual','regression','nlm'});

%%

Index = Max>1.2 & sigTuned;
InteractionTerms = 1:2;
numLags = 20;

% Gather data
Data = gatherROIdata(ROIs,'rawdata');
Neuropil = gatherROIdata(ROIs,'rawneuropil');
NeuropilWeight = cellfun(@(x) x.DataInfo.NeuropilWeight, ROIs, 'UniformOutput',false);
NeuropilWeight = cat(1,NeuropilWeight{:});
Data = Data - bsxfun(@times, NeuropilWeight, Neuropil);

% Create logical stimulus matrix
load([Files{1},'.exp'],'frames','-mat')
nFrames = numel(frames.Stimulus)/numDepths;
% StimMat = createStimMatrix(frames.Stimulus);                % each condition is a unique stimulus
StimMat = createStimMatrix(frames.Stimulus,'Dict',StimLog); % account for interacting stimuli
StimMat = permute(reshape(StimMat,numDepths,nFrames,size(StimMat,2)),[2,3,1]);
StimMat = any(StimMat,3);

% Determine trials to keep
Trial = reshape(frames.Trial,numDepths,nFrames);
Trial = ismember(Trial,TrialIndex{1});
Trial = any(Trial,1);
% Trial = true(numel(frames.Trial),1);

% % Gather run speed
% RunSpeed = mean(reshape(frames.RunningSpeed,numDepths,nFrames));

% Compute GLM
[beta,~,mse,mdl,Filter,~] = fitGLM(StimMat(Trial,:),Data(Index,Trial)',numLags,...
    'InteractionTerms',InteractionTerms,'verbose'); %,...
%             'RunSpeed',AnalysisInfo.onlineRunSpeed(TrialIndex{1}));

rsq=nan(1,nnz(Index));
for ind = 1:nnz(Index)
    rsq(ind) = mdl{ind}.Rsquared.Adjusted;
end
rsq
mse