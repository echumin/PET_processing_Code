% You can either enter these commands manually, or just run this script at the 
% Matlab/Octave prompt.  We're assuming that you've placed the test data in
% the yapmat-0.0.1/demos/data directory and that your current working directory
% is yapmat-0.0.1/demos/commands

img1 = '../data/wrStrong_1337_5c21_de3_01010001.nii';
Cr = '../data/test_mdata_5443-cereb-1-L-box_w-15_2_15--18_-88_-34_roi.txt';
outdir = '../output';
thalf = 20.4;
[BP, k2r, R1, k2, k2a] = mrtm2 (img1, 'hex', Cr, outdir, thalf, 'unweighted');

figure;
view_coronal_slice (BP.*(BP>0 & BP<6), 90, 'jet');
title ('BP, coronal slice #90');

figure;
view_axial_slice (R1, 60);
title ('R1, axial slice #60');

figure;
view_sagittal_slice (k2a, 70, 'hot');
title ('k2a, sagittal slice #70');
