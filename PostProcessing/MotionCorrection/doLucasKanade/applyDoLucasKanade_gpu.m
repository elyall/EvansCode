function [out, B, xi, yi] = applyDoLucasKanade_gpu(Images, dpx, dpy, B, xi, yi)


%% Check input arguments
narginchk(3,6);
    
%% Get dimensions of the data
Nbasis = size(dpx, 1)-1;
[Height, Width, Depth, Channels, totalFrames] = size(Images);

%% Create linear b-splines
if ~exist('B', 'var') || isempty(B)
    warning('off','fastBSpline:nomex');
    knots = linspace(1,Height,Nbasis+1);
    knots = [knots(1)-(knots(2)-knots(1)),knots,knots(end)+(knots(end)-knots(end-1))];
    spl = fastBSpline(knots,knots(1:end-2));
    B = spl.getBasis((1:Height)');
    B = gpuArray(full(B));
end

%% Initialize displacement field
if ~exist('xi', 'var') || ~exist('yi', 'var') || isempty(xi) || isempty(yi)
    [xi,yi] = meshgrid(1:Width,1:Height);
end

%% Cycle through each frame and apply transformation
out = gpuArray(zeros(Height, Width, Depth, Channels, totalFrames));
bl = quantile(Images(:),.01);
for f = 1:totalFrames
    Dx = repmat((B*dpx(:,f)),1,Width);
    Dy = repmat((B*dpy(:,f)),1,Width);
    for d = 1:Depth % applies same transformation to all depths
        for c = 1:Channels
            temp = interp2(gpuArray(double(Images(:,:,d,c,f))),xi+Dx,yi+Dy,'linear');
            temp(isnan(temp)) = bl;
            out(:,:,d,c,f) = temp;
        end
    end
end
out = gather(out);
