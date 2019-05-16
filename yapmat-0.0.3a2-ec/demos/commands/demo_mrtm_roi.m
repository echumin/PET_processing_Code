% You can either enter these commands manually, or just run this script at the 
% Matlab/Octave prompt.  We're assuming that you've placed the test data in
% the yapmat-0.0.2/demos/data directory and that your current working directory
% is yapmat-0.0.2/demos/commands

Ct = '../data/tar.txt';
Cr = '../data/ref.txt';
outdir = '../output';
thalf = 20.4;
[BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
fprintf ('\nA summary report and figures showing model fit to data and');
fprintf ('\nweighted residuals were saved in ''yapmat-0.0.2/demos/output''.\n');
