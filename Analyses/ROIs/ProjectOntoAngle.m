function [X,Y] = ProjectOntoAngle(Points,theta,center)

theta = 2*pi*theta/360;                         % convert to radians

if ~exist('center','var')
    center = mean(Points);
end
Points = bsxfun(@minus, Points, center);        % center points

X = Points*[cos(theta);sin(theta)];             % project onto axis
Y = Points*[cos(theta+pi/2);sin(theta+pi/2)];   % project onto axis