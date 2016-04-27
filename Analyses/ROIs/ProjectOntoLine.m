function ProjPoint = ProjectOntoLine(Vector, Points)

% Normalize vector
Line = diff(Vector, [], 1);     % determine vector
v = Line/norm(Line);            % normalize line

% Subtract off origin
ProjPoint = bsxfun(@minus,Points,Vector(2,:));

% Project points onto line
ProjPoint = newPoints*Line./(Line'*Line)*Line';
% for index = 1:size(ProjPoint,1)
%     ProjPoint(index,:) = dot(ProjPoint(index,:),v)*v;
% end

% Add back origin
ProjPoint = bsxfun(@plus,ProjPoint,Vector(2,:));