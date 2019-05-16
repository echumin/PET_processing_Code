function [BP, k2r, R1, k2, k2a] = mrtm2 (fimg1, base, fCr, outdir, thalf, wmeth, k2ref_mask)

% File: mrtm2.m
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
%   [BP, K2R, R1, K2, K2A] = mrtm2 (IMG1, BASE, CR, OUTDIR, HALFLIFE, WMETH, K2REF_MASK)
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
%   K2REF_MASK  Filename of the binary mask used to determine which voxels to
%               use in the calculation of a global reference region k2. If this
%               input is omitted or the mask file cannot be read, the voxels
%               used in the calculation will be determined automatically by
%               thresholding of the microparameters estimated on the first pass.
%               Should the the threshold method be used, the voxels included in
%               the global reference k2 will be saved out as a binary mask in
%               the image file with "_k2r-calc-mask.nii" postfix.
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
% Version: 0.0.3a2 (rev3)
% Modified: 12 November 2012 by Marc Normandin (normandin@ieee.org)
%           Added subject identifier to saved histogram filename
% Modified: 23 August 2012 by Marc Normandin (normandin@ieee.org)
%           Added option to spatially mask voxels used in global k2ref calcs

% Copyright (c) 2008-2012 Marc D. Normandin
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
lambda = log(2) / thalf;

% read and parse reference region, time frame data
fprintf ('\nReading reference region and time frame data ...');
data = dlmread (fCr, '\t');
ts = data(:,1) / 60; 
te = data(:,2) / 60;  
Cr = data(:,3);
del = te - ts;
nframes = length (ts);
fprintf (' done.');

% read image data
fprintf ('\nReading dynamic image data ...');
[dyn, fstump, hdr, nii] = voxel_reader (fimg1, nframes, base);
fprintf (' done.');
dynm = mean (dyn, 4);
dim = size (dyn);

% do mrtm at each voxel
fprintf ('\nPerforming MRTM at each voxel ...');
[BP, R1, k2, k2a] = mrtm_vox (dyn, Cr, ts, del, lambda, wmeth);
fprintf (' done.');

% write images to file
fprintf ('\nWriting parametric images to file ...');
nii.img = BP;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM_BP.nii')));
nii.img = R1;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM_R1.nii')));
nii.img = k2;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM_k2.nii')));
nii.img = k2a;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM_k2a.nii')));
nii.img = k2 ./ R1;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM_k2r.nii')));
fprintf (' done.');

% calculate a global value for k2r
fprintf ('\nCalculating a global value for k2 in the reference region ...');
nvox = prod (dim(1:3));
BPreshape = reshape (BP, nvox, 1);
R1reshape = reshape (R1, nvox, 1);
k2reshape = reshape (k2, nvox, 1);
k2areshape = reshape (k2a, nvox, 1);
use_mask = 0;
if exist ('k2ref_mask', 'var') 
  try
    nii_mask = load_nii (k2ref_mask);
    mask = nii_mask.img;
    maskreshape = reshape (mask, nvox, 1);
    idx = find (maskreshape == 1);
    use_mask = 1;
    k2ref_method = sprintf ('Mask with file %s', k2ref_mask);
  catch
    fprintf ('\n  *** WARNING *** Couldn''t open file %s as mask for k2ref calculations.', k2ref_mask);
    use_mask = 0;
  end
end
if ~use_mask
  fprintf ('\n  *** WARNING *** Using parameter thresholds to select voxels included in the global calculation.');
  k2ref_method = 'parameter thresholds';
  idx = find (BPreshape>1 & R1reshape>0.5 & k2reshape>0.05 & k2areshape>0.01);
  maskreshape = zeros (size (BPreshape));
  maskreshape(idx) = 1;
  nii.img = reshape (maskreshape, dim(1:3));
  save_nii (nii, fullfile (outdir, sprintf ('%s_MRTM_k2r-calc-mask.nii', fstump)));
end
k2r = median (k2reshape(idx) ./ R1reshape(idx));
fprintf (' done.');

% print histogram of k2r values
figure ('visible', 'off');
hold on;
grid on;
hist (k2reshape(idx) ./ R1reshape(idx), 50); colormap ('white');
h = get (gca);
try
  gtk = graphics_toolkit;
  graphics_toolkit ('gnuplot');
  plot ([k2r, k2r], h.ylim, 'r--', 'LineWidth', 2);
  graphics_toolkit (gtk);
catch
  plot ([k2r, k2r], h.YLim, 'r--', 'LineWidth', 2);
end
plot (k2r, 0, 'ro', 'MarkerEdgeColor', [1,0,0], 'MarkerFaceColor', [1,0,0]);
xlabel ('k_{2}^{ref} (min^{-1})');
ylabel ('Frequency (number of voxels)');
%print ('-dpng', fullfile (outdir, sprintf ('k2r-hist-MRTM-%s.png', wmeth)));
%print ('-depsc', fullfile (outdir, sprintf ('k2r-hist-MRTM-%s.eps', wmeth)));
print ('-dpng', fullfile (outdir, sprintf ('%s_MRTM_k2r-mask-hist.png', fstump)));
print ('-depsc', fullfile (outdir, sprintf ('%s_MRTM_k2r-mask-hist.eps', fstump)));

% do mrtm2 (mrtm with k2r fixed) at each voxel
fprintf ('\nPerforming MRTM (with k2r fixed at %f) at each voxel ...', k2r);
[BP, R1, k2, k2a] = mrtm2_vox (dyn, Cr, ts, del, k2r, lambda, wmeth);
fprintf (' done.');

% write images to file
fprintf ('\nWriting parametric images to file ...');
nii.img = BP;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM2_BP.nii')));
nii.img = R1;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM2_R1.nii')));
nii.img = k2;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM2_k2.nii')));
nii.img = k2a;
save_nii (nii, fullfile (outdir, strcat (fstump, '_MRTM2_k2a.nii')));
fprintf (' done.');

% write summary report
t = toc / 60;
fid = fopen (fullfile (outdir, strcat (fstump, '_MRTM2_report.txt')), 'a');
fprintf (fid, '%s\n', date);
fprintf (fid, '%s\n', 'Input');
fprintf (fid, '\t%s\t%s\n', 'Image file:', fimg1);
fprintf (fid, '\t%s\t%s\n', 'Reference region file:', fCr);
fprintf (fid, '\t%s\t%d\n', 'Number of frames:', nframes);
fprintf (fid, '\t%s\t%s\n', 'Calculation of k2ref:', k2ref_method);
fprintf (fid, '\t%s\t%f\n', 'Tracer half-life (min):', thalf);
fprintf (fid, '\t%s\t%s\n', 'Weighting method:', wmeth);
fprintf (fid, '%s\n', 'Output');
fprintf (fid, '\t%s\t%s\n', 'BP image (MRTM):', ...
         fullfile (outdir, strcat (fstump, '_MRTM_BP.nii')));
fprintf (fid, '\t%s\t%s\n', 'R1 image (MRTM):', ...
         fullfile (outdir, strcat (fstump, '_MRTM_R1.nii')));
fprintf (fid, '\t%s\t%s\n', 'k2 image (MRTM):', ...
         fullfile (outdir, strcat (fstump, '_MRTM_k2.nii')));
fprintf (fid, '\t%s\t%s\n', 'k2a image (MRTM):', ...
         fullfile (outdir, strcat (fstump, '_MRTM_k2a.nii')));
fprintf (fid, '\t%s\t%s\n', 'k2r image (MRTM):', ...
         fullfile (outdir, strcat (fstump, '_MRTM_k2r.nii')));
fprintf (fid, '\t%s%f%s\n', 'Global k2r = ', k2r, '(1/min)');
fprintf (fid, '\t%s\t%s\n', 'BP image (MRTM2):', ...
         fullfile (outdir, strcat (fstump, '_MRTM2_BP.nii')));
fprintf (fid, '\t%s\t%s\n', 'R1 image (MRTM2):', ...
         fullfile (outdir, strcat (fstump, '_MRTM2_R1.nii')));
fprintf (fid, '\t%s\t%s\n', 'k2 image (MRTM2):', ...
         fullfile (outdir, strcat (fstump, '_MRTM2_k2.nii')));
fprintf (fid, '\t%s\t%s\n', 'k2a image (MRTM2):', ...
         fullfile (outdir, strcat (fstump, '_MRTM2_k2a.nii')));
fprintf (fid, '\t%s%f%s\n', 'Runtime = ', t, ' minutes');
fprintf (fid, '%s\n', 'Environment');
fprintf (fid, '\t%s\t%s\n', 'Architecture:', computer);
fprintf (fid, '\t%s\t%s\n', 'Software version:', version);
fclose (fid);

% wrap things up
fprintf ('\n\nDone! Total run time was %f minutes.\n\n', t);
