function [X,Y] = ProjectOntoAngle(Points,theta)

theta = 2*pi*theta/360;                         % convert to radians

Points = bsxfun(@minus, Points, mean(Points));  % center points

X = Points*[cos(theta);sin(theta)];             % project onto axis
Y = Points*[cos(theta+pi/2);sin(theta+pi/2)];   % project onto axis