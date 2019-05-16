function x = mrtm_eng (Cr, Ct, del, w)

% File: mrtm_eng.m
%
% Description:
%   Computational engine that fits the multilinear reference tissue model to
%   data.  The MRTM is given by 
%                           t_             t_
%                           |              |
%     Ct(t) = R1*Cr(t) + k2*|Cr(u)du - k2a*|Ct(u)du
%                          _|             _|
%                           0              0
%   where Ct is PET in the target region, Cr is PET in the reference region,
%   R1 is the ratio of K1 in the target region to K1 in the reference region,
%   k2 is the target efflux rate in the two-tissue compartment sense (i.e.,
%   equal to R1 times k2 in the reference region) and k2a is the apparent
%   efflux rate of the target region data in the one-tissue model sense (i.e.,
%   k2 = k2a*(1+BP), so BP is given by (k2/k2a)-1).
%
% Usage:
%   X = mrtm_eng (CR, CT, DEL, W)
%
% Inputs:
%   CR   Vector giving the time activity curve from the reference region.
%   CT   Vector giving the time activity curve from the target region.
%   DEL  Vector giving the durations (in minutes) of the PET time frames.
%   W    Vector giving the weights (actually, square root of weights) to be
%        applied to residual in weighted least squares fitting.
%
% Outputs:
%   X    Vector of estimated parameters [R1 k2 k2a]'.  R1 is unitless, k2 and
%        k2a have units of 1/min.  BP can be calculated as (k2/k2a)-1.
%
% Dependencies:
%   running_integral.m
%
% See also:
%   mrtm2_eng.m
%   mrtm.m
%   mrtm_vox.m
%   mrtm2.m
%   mrtm2_vox.m
%
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


rintCt = running_integral (del, Ct);
rintCr = running_integral (del, Cr);
A = [Cr rintCr -rintCt];
W = diag (w);
x = (W * A) \ (W * Ct);
%if (any (isnan (x)) | any (x < 0))
%  options = optimset ('MaxIter', 5, 'Display', 'off');
%  [x,resnorm,residual,exitflag,output] = lsqnonneg (W*A, W*Ct, [], options);
%  fprintf (' %d %d;', output.iterations, 1);
%end
%if (any (isnan (x)) | any (x < 0))
%  x = zeros (size (x));
%end
