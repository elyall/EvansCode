function cmap = HiLoColormap(cmin, cmax, color_input)
% assumes cmin is less than 0

if ~exist('color_input', 'var')
    top = [1 0 0]; % red
    middle= [1 1 1]; % black
    bottom = [0 0 1]; % blue
    color_input = [bottom;  middle;  top];
end

colornum = 250;

range = cmax - cmin;
mid = round(colornum*(abs(cmin)/range));

cmap = zeros(colornum, 3);

for cindex = 1:3
    cmap(1:mid,cindex) = linspace(color_input(1,cindex), color_input(2,cindex), mid);
%     cmap(mid, cindex) = color_input(2, cindex);
    cmap(mid:colornum,cindex) = linspace(color_input(2,cindex), color_input(3,cindex), colornum-mid+1);
end
    