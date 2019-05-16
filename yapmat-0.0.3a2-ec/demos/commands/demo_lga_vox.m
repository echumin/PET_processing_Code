fimg1 = '../data/dyndata_01010001.nii';   % first frame of dynamic sequence
base = 'hex';                             % base of frame index ('hex' or 'dec')
fCr = '../data/refdata.txt';              % reference region data file
outdir = '../output/';                    % output directory
tstar = 25;                               % cutoff time (in minutes)
[BP, int] = lga_vox (fimg1, base, fCr, outdir, tstar);  % function call
