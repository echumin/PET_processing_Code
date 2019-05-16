function [BP, int] = lga_vox (fimg1, base, fCin, outdir, tstar)

% File: lga_vox.m
%
% Description:
%   Acts as a wrapper streamlining voxelwise analysis with MRTM2.  First, the 
%   multilinear reference tissue model (MRTM) is applied.  MRTM can be 
%   expressed as
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
%   Next, a global value is determined for k2r (k2 in the reference region, 
%   equal k2 over R1).  Finally, MRTM is re-applied with k2r fixed at this 
%   value (which consitutes MRTM2) on the premise that since there is only one 
%   reference region, there should be only one k2r.  Substituting R1*k2r for 
%   k2, the model now relies on only two floating variables, R1 and k2a:
%                 _           t_      _        t_
%                |            |        |       |
%     Ct(t) = R1*|Cr(t) + k2r*|Cr(u)du | - k2a*|Ct(u)du
%                |_          _|       _|      _|
%                             0                0
%   Reduction in number of parameters improves precision (particularly in low
%   binding regions, where k2 and k2a are highly correlated) at the expense of
%   an increase in bias.  BP is calculated by (R1*k2r/k2a)-1.
%
% Usage:
%   [BP, K2R, R1, K2, K2A] = mrtm2 (IMG1, BASE, CR, OUTDIR, HALFLIFE, WMETH)
%
% Inputs:
%   IMG1      Filename (with path, if necessary) of first frame in the dynamic 
%             image sequence.  Images should be in NIFTI or ANALYZE format, 
%             enumerated by decimal or hexadecimal filename postfixes.
%   BASE      Indicates base of frame enumeration.  Acceptable values are 'hex' 
%             and 'dec'.
%   CR        Filename (with path, if necessary) of reference region data.  The
%             first column should hold the frame start times (in seconds), the
%             second column frame end times (in seconds), and the third column
%             the time activity data for the reference region (same units as the
%             dynamic image 
%   OUTDIR    Directory where all output (images and report files) will be 
%             written.
%   HALFLIFE  The half-life (in minutes) of the tracer isotope.  HINT: 11C has
%             a 20.4 minute half-life and 18F has a 109.8 minute half-life.
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
%   BP        3D array of binding potential values estimated at each voxel in
%             dynamic image sequence.  Dimensions are determined by the size
%             of the input image volumes.
%   K2R       The global value (1/min) used for k2 in the reference region.
%   R1        3D array of R1 values estimated at each voxel.
%   K2 	      3D array of k2 values (units of 1/min) estimated at each voxel.
%   K2A       3D array of k2a values (units of 1/min) estimated at each voxel.
%   
% Dependencies:
%   dlmread  (in Octave, requires Octave-Forge package 'io')
%   voxel_reader.m
%   mrtm_vox.m
%   save_nii.m
%   mrtm2_vox.m
%
% See also:
%   mrtm2_eng.m
%   mrtm2_vox.m
%   mrtm.m
%   mrtm_eng.m
%   mrtm_vox.m
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


warning off; 
tic;
%lambda = log(2) / thalf;

% read and parse reference region, time frame data
fprintf ('\nReading reference region and time frame data ...');
data = dlmread (fCin, '\t');
ts = data(:,1) / 60; 
te = data(:,2) / 60;  
Cin = data(:,3);
del = te - ts;
nframes = length (ts);
fprintf (' done.');

% read image data
fprintf ('\nReading dynamic image data ...');
[dyn, fstump, hdr, nii] = voxel_reader (fimg1, nframes, base);
fprintf (' done.');
dynm = mean (dyn, 4);
dim = size (dyn);

% allocate memory
dim = size (dyn);
BP = zeros (dim(1), dim(2), dim(3));
int = zeros (dim(1), dim(2), dim(3));

% prep some variables outside loop
rint_Cin = running_integral (del, Cin);
tm = cumsum (del) - del/2;
istar = min (find (tm >= tstar));
if length (istar) == 0
  fprintf ('\n');
  error ('Cutoff time t* exceeds the length of the scan!');
end

% do lga at each voxel (can this be vectorized to speed things up?)
fprintf ('\nPerforming Logan graphical analysis at each voxel ...');
for x = 1:dim(1)
  fprintf ('\n  Logan graphical analysis on sagittal slice number %d of %d ...', x, dim(1));
  for y = 1:dim(2)
    for z = 1:dim(3)
%      if dynm (x,y,z) > medval
% fprintf ('\n      Analyzing voxel (%d %d %d).', x, y, z);
        Ct = reshape (dyn(x,y,z,:), dim(4), 1);
        if ~(any (Ct == 0) | any (isnan (Ct)))
%        if x == 2 & y == 70 & z == 58
%          save ('-binary', 'data.obin');
%          keyboard;
%        end
        p = lga_eng (rint_Cin, Ct, del, istar);
      else
% fprintf ('\n      Skipping voxel (%d %d %d).', x, y, z);
%        p = [0; 0];
        p = [NaN; NaN];
      end
      BP(x,y,z) = p(1) - 1;
      int(x,y,z) = p(2);
    end
  end
  fprintf (' done.');
end

% instead of DV, just call it 'm' or 'slope'
%[BP, int] = lga_vox (dyn, Cin, ts, del, tstar);
fprintf (' done.');

% write images to file
fprintf ('\nWriting parametric images to file ...');
nii.img = BP;
save_nii (nii, fullfile (outdir, strcat (fstump, '_LGA_BP.nii')));
nii.img = int;
save_nii (nii, fullfile (outdir, strcat (fstump, '_LGA_int.nii')));
fprintf (' done.');

% write summary report
t = toc / 60;
fid = fopen (fullfile (outdir, strcat (fstump, '_LGA_report.txt')), 'a');
fprintf (fid, '%s\n', date);
fprintf (fid, '%s\n', 'Input');
fprintf (fid, '\t%s\t%s\n', 'Image file:', fimg1);
fprintf (fid, '\t%s\t%s\n', 'Input function file:', fCin);
%fprintf (fid, '\t%s\t%s\n', 'Input type:', ref_or_plas);
fprintf (fid, '%s\n', 'Output');
%fprintf (fid, '\t%s\t%s\n', 'DV image:', ...
%         fullfile (outdir, strcat (fstump, '_LGA_DV.nii')));
fprintf (fid, '\t%s\t%s\n', 'BP image:', ...
         fullfile (outdir, strcat (fstump, '_LGA_BP.nii')));
fprintf (fid, '\t%s\t%s\n', 'int image:', ...
         fullfile (outdir, strcat (fstump, '_LGA_int.nii')));
fprintf (fid, '\t%s%f%s\n', 'Runtime = ', t, ' minutes');
fprintf (fid, '%s\n', 'Environment');
fprintf (fid, '\t%s\t%s\n', 'Architecture:', computer);
fprintf (fid, '\t%s\t%s\n', 'Software version:', version);
fclose (fid);

% wrap things up
fprintf ('\n\nDone! Total run time was %f minutes.\n\n', t);
