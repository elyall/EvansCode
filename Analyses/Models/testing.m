
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

