function [ProjPoint, ProjPoint_orig] = ProjectOntoLine(Vector, Points)

% subtract off origin
newPoints = bsxfun(@minus, Points, Vector(1,:));
Line = diff(Vector, [], 1)';

% project point onto line
ProjPoint = newPoints*Line./(Line'*Line)*Line';

% add back origin
ProjPoint_orig = bsxfun(@plus, ProjPoint, Vector(1,:));