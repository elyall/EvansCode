Data = zeros(3,3,32);
loc = [2,2;2,1;1,2;3,2;2,3];
Loc = sub2ind([3,3],loc(:,1),loc(:,2));
for s = 1:totalStim
    curr = zeros(3,3);
    curr(Loc(StimLog(s,:))) = 1;
    Data(:,:,s) = curr;
end

Data = repmat(Data,1,1,20); % multiple trials for each stimulus
Data = rand(size(Data)).*Data;

[Phi,ahat] = sparsenet(100,.21,Data,0,0);