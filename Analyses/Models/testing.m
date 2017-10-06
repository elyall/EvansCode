
% Gather data
[~,order] = sort(Max,'descend');
a = Data(order(1),:);
% a = Data(order(2),:);
a = a - prctile(a,30);

% Create stim matrix
x = StimMat;
x = designStimMatrix(StimMat,1:2);

% Add run speed
RunSpeed = mean(reshape(frames.RunningSpeed,numDepths,nFrames))';
x = [x, RunSpeed];
RunSpeed = RunSpeed > 100;
% x = [bsxfun(@times, x, ~RunSpeed),bsxfun(@times, x, RunSpeed)];
x(:,1:end-1) = bsxfun(@times, x(:,1:end-1), RunSpeed);

% Lag stims
numLags = 10;
n = size(x,2);
x = lagmatrix(x,0:numLags);
% a(1:10) = [];
% x(1:10,:) = [];
x(isnan(x)) = 0;

% Regress
% d = a'\[ones(size(x,1),1),x];
% d = a'\x;
d = regress(a',x);

% STRF
% f = reshape(d(2:end),5,11);
f = reshape(d,n,numLags+1);
figure; imagesc(f);

% Prediction
% p = [ones(size(x,1),1),x]*d;
p = x*d;
% p = x(:,1:2:end)*d(1:2:end);
% p2 = x(:,56:end)*d(56:end);
mse = mean((p-a').^2)

figure;
% plot(zscore(a));
% hold on;
% plot(zscore(p));
plot(a);
hold on
plot(p);

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