% GLM -trial

% Gather data
[~,order] = sort(M,'descend');
a = Data(order(1),:);
a = a - prctile(a,30);
x = designStimMatrix(StimMat,[1,2]); % Create stim matrix
% good = true(numel(a),1);
x = [true(size(x,1),1),x]; % Add bias

% Remove non-running trials
RunSpeed = AnalysisInfo.onlineRunSpeed;
% good = RunSpeed>100;
% x = x(good,:);
% a = a(good);
% RunSpeed = RunSpeed(good);


%% look at dataz

figure;
plot(RunSpeed,a,'.')

figure;
plotSpread(a,'distributionIdx',AnalysisInfo.StimID(good));


%% Linear fit

d = regress(a',x);
p = x*d;
figure;
plot(d);

d = fitglm(x,a'+2, 'Distribution','gamma','Intercept',false);
p = d.predict - 2;
figure;
plot(d.Coefficients.Estimate);

figure;
plot(a,p,'.');

mse = mean((p-a').^2)


%% Add nonlinearity

beta0 = [zeros(1,size(x,2)), 2];
% beta0 = [d; 0];
mdl = fitnlm(x, a, @div_norm, beta0);

figure
plot(a,mdl.predict,'.')


%% Population
dindex = 4;

Data = cat(2,StimMean{DataIndex==dindex});
Data = Data(TrialIndex{dindex},:);

ID = StimID{dindex}(TrialIndex{dindex});
x = StimLog(ID+1,:);

% b = fitglm(Data,x);

Vec = arrayfun(@(x) mean(Data(ID==x,:)), unique(ID), 'UniformOutput', false);
Vec = cat(1,Vec{:})';

lenC = dot(Vec(:,7),Vec(:,2))/norm(Vec(:,2));  % EDITED, twice
C = lenC*Vec(:,2);

