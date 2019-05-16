function [p, A, y, x, yy, xx] = lga_voxeng (rint_Cin, Ct, del, istar)
%function [p, istar, A, y, x, yy, xx] = lga_voxeng (rint_Cin, Ct, del, tm, istar)
% change this function name to lga_eng

% tm not used in optimized (rint_Cin and istar passed directly) function
% istar is an input, needn't be returned
% aren't going to care about xx and yy (used for plotting) with voxel analyses

% *** Ci CONFLICTS WITH GSL LIBRARY BINDINGS!!!! *****

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
