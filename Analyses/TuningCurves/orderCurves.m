function out = orderCurves(Curves,StimLog,StimMean,StimID,DataIndex,TrialIndex)

numWhiskersPerStim = sum(StimLog,2);
[Max1W,BW] = max(Curves(:,numWhiskersPerStim==1),[],2);

[numROIs,numStim] = size(Curves);
nW = unique(numWhiskersPerStim);
numConditions = numel(nW);

[~,numWhiskers] = size(StimLog);
% Labels = 1:numWhiskers;
% Labels = arrayfun(@(x) Labels(StimLog(x,:)), 1:numStim, 'UniformOutput',false); % gather relevant pieces

out = Curves;
x = [ones(numStim,1),StimLog];
for r = 1:numROIs
    [~,order] = sort(Curves(r,numWhiskersPerStim==1),2,'ascend');
%     weights = StimMean{r}(TrialIndex{DataIndex(r)})\x(StimID{DataIndex(r)}(TrialIndex{DataIndex(r)})+1,:);
%     [~,order] = sort(weights(2:end),'ascend');
    order = order.^3; % ensure each is weighted properly
    for c = nW(2:end-1)'
        [~,curr] = sort(sum(StimLog(numWhiskersPerStim==c,:).*order,2),'descend');
        temp = Curves(r,numWhiskersPerStim==c);
        out(r,numWhiskersPerStim==c) = temp(curr);
    end
end