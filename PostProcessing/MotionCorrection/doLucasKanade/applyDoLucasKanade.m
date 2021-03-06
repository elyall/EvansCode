function [out, B, xi, yi] = applyDoLucasKanade(Images, dpx, dpy, Nbasis, B, xi, yi)


%% Check input arguments
narginchk(3,7);
if ~exist('Nbasis', 'var') || isempty(Nbasis)
    Nbasis = 16;
end

%% Get dimensions of the data
[Height, Width, Depth, Channels, Frames] = size(Images);

%% Create linear b-splines
if ~exist('B', 'var') || isempty(B)
    warning('off','fastBSpline:nomex');
    knots = linspace(1,Height,Nbasis+1);
    knots = [knots(1)-(knots(2)-knots(1)),knots,knots(end)+(knots(end)-knots(end-1))];
    spl = fastBSpline(knots,knots(1:end-2));
    B = spl.getBasis((1:Height)');
    B = full(B);
end

%% Initialize displacement field
if ~exist('xi', 'var') || ~exist('yi', 'var') || isempty(xi) || isempty(yi)
    [xi,yi] = meshgrid(1:Width,1:Height);
end

%% Cycle through each frame and apply transformation
out = zeros(Height, Width, Depth, Channels, Frames);
bl = quantile(Images(:),.01);
for f = 1:Frames
    Dx = repmat((B*dpx(:,f)),1,Width);
    Dy = repmat((B*dpy(:,f)),1,Width);
    for d = 1:Depth % assumes different depths to be captured simultaneously (e.g. Adrian Cheng's technique)
        for c = 1:Channels
            temp = interp2(double(Images(:,:,d,c,f)),xi+Dx,yi+Dy,'linear');
            temp(isnan(temp)) = bl;
            out(:,:,d,c,f) = temp;
        end
    end
end