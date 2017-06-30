function [ProjPoints,orthoProj,newBasis] = ProjectOntoLine(Vector, Points)
% Vector is [x1, y1; x2, y2]
% Points is N x 2

% Normalize vector
v = diff(Vector, [], 1); % determine vector
v = v/norm(v);           % normalize line

% Subtract off origin
centeredPoints = bsxfun(@minus, Points, Vector(1,:));

% Project points onto line
ProjPoints = bsxfun(@times, centeredPoints*v', v);

% Add back origin
ProjPoints = bsxfun(@plus, ProjPoints, Vector(1,:));

% Determine orthogonal projection
orthoProj = Points - ProjPoints;
orthoProj = bsxfun(@plus, orthoProj, Vector(1,:));

newBasis = centeredPoints*null(v(:).');


% v = (V2-V1)/norm(V2-V1); %// normalized vector from V1 to V2
% Q = dot(P-V1,v)*v+V1; %// projection of P onto line from V1 to V2
% dist = norm(P-Q);
% alfa = (Q(1)-V1(1))/(V2(1)-V1(1));