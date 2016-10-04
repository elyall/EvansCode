directory = '/home/elyall/Documents/Data/2509/151029/0000/';
[files,p] = uigetfile({'*.sbx'},'Select files',directory,'MultiSelect','on');
files = fullfile(p,files);
numFiles = numel(files);

[imgs, loadObj] = load2P(files);
out = zeros([loadObj.size(1:2),numFiles]);
n = 1;
for findex = 1:numFiles
    c = loadObj.files.Frames;
    out(:,:,findex) = mean(imgs(:,:,1,1,n:n+c-1),5);
    n = n + c - 1;
end