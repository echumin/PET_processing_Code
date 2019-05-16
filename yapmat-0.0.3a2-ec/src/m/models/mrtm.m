function [BP, R1, k2, k2a, k2r] = mrtm (fCt, fCr, outdir, thalf, wmeth)

% File: mrtm.m
%
% Description:
%   Apply multilinear reference tissue model (MRTM) to region of interest data.
%   The MRTM can be expressed as
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
%   [BP, R1, K2, K2A, K2R] = mrtm (CT, CR, OUTDIR, HALFLIFE, WMETH)
%
% Inputs:      
%   CT        Filename (with path, if necessary) of target region data.  The
%             first column should hold the frame start times (in seconds), the
%             second column frame end times (in seconds), and the third column
%             the time activity data for the target region.
%   CR        Filename of reference region data, formatted in same fashion as
%             target region file.
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
%   BP        Estimated binding potential (unitless).
%   R1        Estimated R1 (unitless).
%   K2        Estimated k2 (1/min) in the two-tissue model sense.
%   K2A       Estimated k2a (1/min), the apparent k2 in the one-tissue model 
%             sense.
%   K2R       Estimated k2 in the reference region (1/min).
%
% Dependencies:
%   dlmread  (in Octave, requires Octave-Forge package 'io')
%   compute_weights.m
%   running_integral.m
%   mrtm_eng.m
%
% See also:
%   mrtm_eng.m
%   mrtm_vox.m
%   mrtm2.m
%   mrtm2_eng.m
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


tic;
lambda = log (2) / thalf;

% read and parse target region, time frame data
disp('Reading target region and time frame data ...');
[~, fname, ~] = fileparts (fCt);
data = dlmread (fCt, '\t');
ts = data(:,1) / 60; 
te = data(:,2) / 60;  
Ct = data(:,3);
del = te - ts;
tm = ts + del/2;
nframes = length (ts);
disp(' done.');

% read and parse reference region, time frame data
disp('Reading reference region and time frame data ...');
data = dlmread (fCr, '\t');
if ~isequal (ts,data(:,1) / 60) | ~isequal (te, data(:,2) / 60)
  error ('Target and reference region data must have same frame times.');
end
Cr = data(:,3);
disp(' done.');

% compute weights
disp('Computing weights ...');
switch wmeth
  case 'unweighted'
    w = ones (nframes, 1);
  case 'conventional'
    w = compute_weights (Ct, ts, del, lambda, [0; 1; 0.5]);
    w = w .* exp (-lambda * (ts + del/2));  % transform weights (dk-corr)
  otherwise
    error ('Weighting method not recognized.');
end
disp(' done.');

% do mrtm on curves
disp('Performing MRTM ...');
p = mrtm_eng (Cr, Ct, del, w);
Ctfit = [Cr running_integral(del, Cr) -running_integral(del, Ct)] * p;
resid = Ct - Ctfit;
wresid = w .* resid;
R1 = p(1);
k2 = p(2);
k2a = p(3);
BP = k2 / k2a - 1;
k2r = k2 / R1;
disp(' done.');

% plot data and model fit
disp('Plotting data and model fit ...');
h1=figure;
plot (tm, Cr, 'ro', tm, Ct, 'bo', tm, Ctfit, 'b');
xlabel ('Time (min)');
ylabel ('Tracer Concentration');
%legend ('C_{R}(t)', 'C_{T}(t)', 'Model fit');
legend ('Cr(t)', 'Ct(t)', 'Model fit');
modFig=fullfile (outdir, strcat (fname, '_MRTM_modelfit.fig'));
savefig(h1,modFig)
%print (modFig, '-dpng');
pause(1)
close all
%fprintf (' done.');

% plot weighted residuals
disp('Plotting weighted residuals ...');
h2=figure;
plot (tm, wresid, 'bo', tm, zeros(nframes, 1), 'k');
xlabel ('Time (min)');
ylabel ('Weighted residuals');
resFig=fullfile (outdir, strcat (fname, '_MRTM_wresid.fig'));
savefig(h2,resFig)
%print (resFig, '-dpng');
pause(1)
close all
%fprintf (' done.');

% write summary report
t = toc;
disp('Writing summary report ...');
fid = fopen (fullfile (outdir, strcat (fname, '_MRTM_report.txt')), 'a');
fprintf (fid, '%s\n', date);
fprintf (fid, '%s\n', 'Input');
fprintf (fid, '\t%s\t%s\n', 'Target region file:', fCt);
fprintf (fid, '\t%s\t%s\n', 'Reference region file:', fCr);
fprintf (fid, '\t%s\t%f\n', 'Tracer half-life (min):', thalf);
fprintf (fid, '\t%s\t%s\n', 'Weighting method:', wmeth);
fprintf (fid, '%s\n', 'Output');
fprintf (fid, '\t%s\t%f\n', 'BP (unitless):', BP);
fprintf (fid, '\t%s\t%f\n', 'R1 (unitless):', R1);
fprintf (fid, '\t%s\t%f\n', 'k2 (1/min):', k2);
fprintf (fid, '\t%s\t%f\n', 'k2a (1/min):', k2a);
fprintf (fid, '\t%s\t%f\n', 'k2r (1/min):', k2r);
fprintf (fid, '\t%s\t%f\n', 'Runtime (sec):', t);
fprintf (fid, '%s\n', 'Environment');
fprintf (fid, '\t%s\t%s\n', 'Architecture:', computer);
fprintf (fid, '\t%s\t%s\n', 'Software version:', version);
fclose (fid);
disp(' done.');

% wrap things up
% fprintf ('\nDone! Total run time was %f seconds.\n', t);
% fprintf ('\nBP = %f', BP);
% fprintf ('\nR1 = %f', R1);
% fprintf ('\nk2 = %f (1/min)', k2);
% fprintf ('\nk2a = %f (1/min)', k2a);
% fprintf ('\nk2r = %f (1/min)\n', k2r);
