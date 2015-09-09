function [Avg,Var,T] = sbxAlignPar(Images,thestd,gl,l,varargin)


FrameIndex = [1 inf];

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
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

if FrameIndex(end) == inf
    FrameIndex = [FrameIndex(1:end-1), FrameIndex(1:end-1)+1:totalFrames];
end
numFrames = numel(FrameIndex);


p = gcp();
nblocks = 2^floor(log2(p.NumWorkers));              % determine number of CPUs to divide computation over

rg = 0:numFrames-1;
rgs = spliteven(rg,log2(nblocks));

if loadImages
    img = load2P(Images, 'Frames', FrameIndex(1));
else
    img = Images(:,:,1,1,1);
end

rg1 = 33:size(img,1);
rg2 = 46:size(img,2)-45;



thestd = thestd(rg1,rg2);



img = img(rg1,rg2);

l = l(rg1,rg2);



ms = zeros([size(img),nblocks]);

vs = zeros([size(img),nblocks]);

Ts = cell(nblocks,1);



parfor ii = 1:nblocks
    
    subrg = rgs{ii};
    
    [ms(:,:,ii),vs(:,:,ii),Ts{ii}] = sbxalignsub(fname,subrg,rg1,rg2,thestd,gl,l);
    
end



%non-recursive

for nn = 1:log2(nblocks)
    
    nblocksafter = 2^(log2(nblocks)-nn);
    
    msnew = zeros([size(img),nblocksafter]);
    
    vsnew = zeros([size(img),nblocksafter]);
    
    Tsnew = cell(nblocksafter,1);
    
    
    
    for ii = 1:nblocksafter
        
        [u,Var] = fftalign(ms(:,:,ii*2-1)/length(Ts{ii*2-1}), ms(:,:,ii*2  )/length(Ts{ii*2  }));
        
        
        
        Tsnew{ii} = [bsxfun(@plus,Ts{ii*2-1},[u,Var]);
            
        Ts{ii*2}];
    
    Ar = circshift(ms(:,:,ii*2-1),[u, Var]);
    
    msnew(:,:,ii) = (Ar+ms(:,:,ii*2  ));
    
    
    
    Ar = circshift(vs(:,:,ii*2-1),[u, Var]);
    
    vsnew(:,:,ii) = (Ar+vs(:,:,ii*2  ));
    
    
    
    end
    
    
    
    ms = msnew;
    
    vs = vsnew;
    
    Ts = Tsnew;
    
end



Avg = ms/info.max_idx;

Var = vs/info.max_idx;

T = Ts{1};

end



function A = spliteven(idx,ns)

if ns > 0
    idx0 = idx(1:floor(end/2));
    idx1 = idx(floor(end/2)+1 : end);
    idx0s = spliteven(idx0,ns-1);
    idx1s = spliteven(idx1,ns-1);
    A = {idx0s{:},idx1s{:}};
    
else
    A = {idx};
    
end
end %spliteven



function [m,v,T] = sbxalignsub(fname,idx,rg1,rg2,thestd,gl,l)

if(length(idx)==1)
    
    
    
    A = double(sbxreadpacked(fname,idx(1),1));
    
    A = A(rg1,rg2);
    
    A = A - gl(idx+1)*l;
    
    A = A./thestd;
    
    
    
    m = A;
    
    v = A.^2;
    
    T = [0 0];
    
    
    
elseif (length(idx)==2)
    
    
    
    A = double(sbxreadpacked(fname,idx(1),1));
    
    A = A(rg1,rg2);
    
    A = A - gl(idx(1)+1)*l;
    
    A = A./thestd;
    
    
    
    B = double(sbxreadpacked(fname,idx(2),1));
    
    B = B(rg1,rg2);
    
    B = B - gl(idx(2)+1)*l;
    
    B = B./thestd;
    
    
    
    [u,v] = fftalign(A,B);
    
    
    
    Ar = circshift(A,[u,v]);
    
    m = Ar+B;
    
    T = [[u v] ; [0 0]];
    
    
    
    Ar = circshift(A.^2,[u,v]);
    
    v = (Ar+B.^2);
    
    
    
else
    
    
    
    idx0 = idx(1:floor(end/2));
    
    idx1 = idx(floor(end/2)+1 : end);
    
    [A,v1,T0] = sbxalignsub(fname,idx0,rg1,rg2,thestd,gl,l);
    
    [B,v2,T1] = sbxalignsub(fname,idx1,rg1,rg2,thestd,gl,l);
    
    
    
    [u,v] = fftalign(A/length(idx0), B/length(idx1));
    
    
    
    Ar = circshift(A,[u, v]);
    
    m = (Ar+B);
    
    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
    
    
    v1 = circshift(v1,[u, v]);
    
    v = v1+v2;
    
    
    
end

end