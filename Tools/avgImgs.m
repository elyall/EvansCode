function Images = avgImgs(Images,N,MCdata)

direction = 'descend';
saveOut = true;
saveFile = 'temp.tif';

%% Parse input arguments & load data
if ~exist('Images','var')
    Images = [];
end
if ischar(Images) || iscellstr(Images)
    [Images,loadObj] = load2P(Images,'Frames',[1,inf]);
end

if ~exist('N','var') || isempty(N)
    N = min(1000,size(Images,5));
end

if ~exist('MCdata','var')
    MCdata = [];
end
if ~isempty(MCdata)
    Images = applyMotionCorrection(Images, MCdata, loadObj);
end


%% Compute average
if N~=size(Images,5)
    Images = sort(Images,5,direction);
end
Images = mean(Images(:,:,:,:,1:N),5);


%% Save output
if saveOut
    save2P(saveFile,Images);
end