function [T, Avg, Q] = sbxAlign(Images, varargin)


Frames = [1 inf];
Channel = 1;
Depth = 1;

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Frame','Frames'}
                Frames = varargin{index+1};
                index = index + 2;
            case 'Channel'
                Channel = varargin{index+1};
                index = index + 2;
            case 'Depth'
                Depth = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Images', 'var') || isempty(Images)
    [Images,p] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select images:',directory,'MultiSelect','on');
    if isnumeric(Images)
        return
    elseif iscell(Images)
        for findex = 1:numel(Images)
            Images{findex} = fullfile(p, Images{findex});
        end
    elseif ischar(Images)
        Images = {fullfile(p, Images)};
    end
elseif ischar(Images)
    Images = {Images};
end


%% Load in images
if iscellstr(Images)
    ImageFiles = Images;
    loadImages = true;
else
    loadImages = false;
end


%% Determine frames to analyze
if loadImages
    Config = load2PConfig(ImageFiles);
    totalFrames = sum([Config(:).Frames]);
else
    totalFrames = size(Images, 5);
end

if Frames(end) == inf
    Frames = [Frames(1:end-1), Frames(1:end-1)+1:totalFrames];
end
numFrames = numel(Frames);


%% Get first-order stats
if loadImages
    Images = load2P(ImageFiles, 'Frames', Frames(1));
end
[H, W, D, C, ~] = size(Images);

ms = 0;
vs = 0;
X = [ones(numFrames,1),linspace(-1,1,numFrames)'];
X = bsxfun(@times,X,1./sqrt(sum(X.^2)));
parfor jj = 1:numFrames
    if loadImages
        img = load2P(ImageFiles, 'Frames', Frames(jj), 'Channel', Channel, 'Depth', Depth, 'Double');
    else
        img = Images(:,:,1,1,jj);
    end
    ms = ms + img(:)*X(jj,:);
    vs = vs + img(:).^2;
end
s = reshape(sqrt(1/numFrames*(vs - sum(ms.^2,2))),[H,W]);
thestd = medfilt2(s,[31,31],'symmetric');
gl = X(:,2);
l = reshape(ms(:,2),[H,W]);


%% Alignment first pass
[Avg,~,T] = sbxAlignpar(fname,thestd,gl,l); %Takes about 2.5 minutes

%% Alignment second pass
rgx = (1:size(Avg,2))+45;
rgy = 32 + (1:size(Avg,1));
T0 = T;
for nn = 1:10
    fprintf('Refining alignment... pass %d\n',nn);
    [Avg,~,T] = sbxAligniterative(fname,Avg,rgy,rgx,thestd(rgy,rgx),gl,l);
    dT = sqrt(mean(sum((T0-T).^2,2)));
    T0 = T;
    
    if dT < .25
        break;
    end
    fprintf('delta: %.3f\n',dT);
end

%% Get first-order stats
ms = 0;
vs = 0;
m2 = 0;
v2 = 0;

X = [ones(numFrames,1),linspace(-1,1,numFrames)'];
X = bsxfun(@times,X,1./sqrt(sum(X.^2)));
g = exp(-(-5:5).^2/2/1.6^2);
parfor jj = 1:numFrames
    z = single(sbxreadpacked(fname,jj-1,1));
    z = z./thestd;
    z = circshift(z,T(jj,:));
    z = double(z);

    ms = ms + z(:)*X(jj,:);   
    vs = vs + z(:).^2;
        
    z = conv2(g,g,z,'same');
    m2 = m2 + z(:)*X(jj,:);    
    v2 = v2 + z(:).^2;    
end
ss = sqrt(1/numFrames*(vs - sum(ms.^2,2)));
Avg = reshape(ms(:,1)*X(1,1),size(l));
Var = reshape(ss.^2,size(l));

[Q,~] = qr(ms,0);
[Q2,~] = qr(m2,0);

s2 = reshape(sqrt(1/numFrames*(v2 - sum(m2.^2,2))),size(l));
m2 = reshape(m2(:,1)*X(1,1),size(l));
k = 0;

%% Compute simple stats
parfor jj = 1:numFrames
    z = double(sbxreadpacked(fname,jj-1,1));
    z = z - gl(jj)*l;
    z = circshift(z,T(jj,:));
    z = z./thestd;
        
    z = conv2(g,g,z,'same');    
    z = reshape(z(:) - (Q2*(Q2'*z(:))),size(z));    
    k  =  k + (z./s2).^4;   
end

sm = s2./m2;
k = k/numFrames - 3;

try  
    save([fname '.align'],'Avg','Var','thestd','sm','k','T','Q','-append');
catch
    save([fname '.align'],'Avg','Var','thestd','sm','k','T','Q');
end
