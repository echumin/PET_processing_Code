function [p, A, x, y, xx, yy] = lga_eng (rint_Cin, Ct, del, istar)

% Version: 0.0.3
% Modified: 20 October 2008
% Author: Marc Normandin (normandin@ieee.org)

% Copyright (C) 2008  Marc D. Normandin
%
% This file is part of yapmat.
%
% yapmat is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% yapmat is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with yapmat.  If not, see <http://www.gnu.org/licenses/>.


% Customize for voxel analysis
% (does this make enough difference in speed to justify loss of modularity?)
y = running_integral (del, Ct) ./ Ct;
x = rint_Cin ./ Ct;
xx = x(istar:end);
yy = y(istar:end);
A = [xx ones(length(xx), 1)];
p = A \ yy;

% y = running_integral (del, Ct) ./ Ct;
% x = running_integral (del, Cr) ./ Ct;
% tm = cumsum (del) - del/2;
% istar = min (find (tm >= tstar));
% if length (istar) == 0
%   fprintf ('\n');
%   error ('Cutoff time t* exceeds the length of the scan!');
% end
% xx = x(istar:end);
% yy = y(istar:end);
% A = [xx ones(length(xx), 1)];
% p = A \ yy;
