function [cmap,dim] = colormap2D(Colors,pts,dim)

verbose = false;

if ~exist('Colors','var') || isempty(Colors)
    Colors = [0,0,0;1,0,0;1,1,1;0,0,1]; % black, red, white, blue
end
if ~exist('pts','var') || isempty(pts)
    pts = [0,0;1,0;0,1;1,1]; % top-left, top-right, bottom-left, bottom-right
end
if ~exist('dim','var') || isempty(dim)
    dim = [19,11];
end

%% Create color channels
[X,Y] = meshgrid(0:1/(dim(2)-1):1,0:1/(dim(1)-1):1);
r = scatteredInterpolant(pts(:,1),pts(:,2),Colors(:,1));
R = r(X,Y);
g = scatteredInterpolant(pts(:,1),pts(:,2),Colors(:,2));
G = g(X,Y);
b = scatteredInterpolant(pts(:,1),pts(:,2),Colors(:,3));
B = b(X,Y);

%% Create colormap and shape it to be Nx3
cmap = reshape(cat(3,[R,G,B]),[prod(dim),3]);
cmap(cmap<0) = 0; cmap(cmap>1) = 1; % fix any interpolation errors

%% visualize colormap
if verbose
    img = reshape(1:prod(dim),dim);
    figure; image(img); colormap(cmap);
end