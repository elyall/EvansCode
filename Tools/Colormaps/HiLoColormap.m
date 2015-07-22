function cmap = HiLoColormap(cmin, cmax, color_input)
% assumes cmin is less than 0

if ~exist('color_input', 'var')
    green_top     = [0 1 0];
    black_middle= [0 0 0];
    red_bottom = [1 0 0];
    color_input = [red_bottom;  black_middle;  green_top];
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
    