function [BP, R1, k2, k2a] = mrtm2_vox (dyn, Cr, ts, del, k2r, lambda, wmeth)

% File: mrtm2_vox.m
%
% Description:
%   Facilitates the application of the multilinear reference tissue model with
%   k2 in the reference region fixed (MRTM2) at each voxel in a dynamic image.
%   For more details on MRTM2, please refer to the documentation for mrtm2.m.
%
% Usage:
%   [BP, R1, K2, K2A] = mrtm2_vox (DYN, CR, TS, DEL, K2R, LAMBDA, WMETH)
%
% Inputs:
%   DYN     4D volume of dynamic image data.  The 4th dimension should
%           represent time, and the first three spatial orientation.
%   CR      Vector giving the time activity curve from the reference region.
%   TS      Vector giving the start times (in minutes) of the PET time frames.
%   DEL     Vector giving the durations (in minutes) of the PET time frames.
%   K2R     The value (1/min) at which to fix k2 in the reference region.
%   LAMBDA  Decay constant (in units of 1/min) for tracer isotope, can be
%           calculated as natural logarithm of 2 divided by the isotope half-
%           life (in Matlab/Octave parlance, 'log(2)/t_half').  HINTS: The
%           half-life of 11C is 20.4 min and of 18F is 109.8 min.
%   WMETH     Weighting method used to weight the residuals during parameter
%             estimation.  Acceptable values are 'unweighted', 'conventional'.
%             Methods to be added are 'reference', 'frame-conv', 'frame-ref',
%             'iter-conv', 'iter-ref'.
%               unweighted     equal weighting of all residuals
%               conventional   weights inversely proportional to variance of
%                              target region data (Poisson noise model)
%               reference      weights inversely proportional to variance of
%                              residual, including input contribution (Poisson)
%               frame-conv     variance determined from total frame counts
%                              (Poisson), weights from target only
%               frame-ref      variance determined from total frame counts
%                              (Poisson), weights account for input noise
%               iter-conv      iteratively reweighted (variance from model fit
%                              at last iteration), target contribution only
%               iter-ref       iteratively reweighted (variance from model fits
%                              at last iteration), account for input noise
%
% Outputs:
%   BP      3D array of binding potential values estimated at each voxel in
%           dynamic image sequence.  Dimensions are determined by the size
%           of the input image volumes (first three dimensions of DYN).
%   R1      3D array of R1 values estimated at each voxel.
%   K2      3D array of k2 values (units of 1/min) estimated at each voxel.
%   K2A     3D array of k2a values (units of 1/min) estimated at each voxel.
%
% Dependencies:
%   compute_weights.m
%   mrtm2_eng.m
%
% See also:
%   mrtm_vox.m
%   mrtm2_eng.m
%   mrtm2.m
%   mrtm_eng.m
%   mrtm.m
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


% allocate memory
dim = size (dyn);
BP = zeros (dim(1), dim(2), dim(3));
R1 = zeros (dim(1), dim(2), dim(3));
k2 = zeros (dim(1), dim(2), dim(3));
k2a = zeros (dim(1), dim(2), dim(3));

% loop over all voxels (can this be vectorized to speed things up?)
for x = 1:dim(1)
  fprintf ('\n  MRTM2 on sagittal slice number %d of %d ...', x, dim(1));
  for y = 1:dim(2)
    for z = 1:dim(3)
      Ct = reshape (dyn(x,y,z,:), dim(4), 1);
% THIS IS A SLOPPY HACK TO IMPLEMENT SWITCH ON FITTING
% NEED TO ADD TO ARGS TO compute_weights:
%   1. wmeth - tells HOW weights are calculated
%   2. xfrm_wts - tells whether or not weights should be transformed to
%                 match transformation (i.e., decay-correction) of PET data
% CHECK IF FRAME-BASED WEIGHTS BEFORE LOOPING
      switch wmeth
        case 'unweighted'
          w = ones (dim(4), 1);
        case 'conventional'
          w = compute_weights (Ct, ts, del, lambda, [0; 1; 0.5]);
          w = w .* exp (-lambda * (ts + del/2));  % transform weights (dk-corr)
        otherwise
         error ('Weighting method not recognized.');
      end
% END HACK
      p = mrtm2_eng (Cr, Ct, del, w, k2r);
      R1(x,y,z) = p(1);
      k2(x,y,z) = p(1) * k2r;
      k2a(x,y,z) = p(2);
      BP(x,y,z) = (p(1) * k2r / p(2)) - 1;
    end
  end
  fprintf (' done.');
end
