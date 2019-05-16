fCt = '../data/tardata.txt';           % target region data file 
fCr = '../data/refdata.txt';           % reference region data file
outdir = '../output/';                 % output directory
tstar = 25;                            % cutoff time (in minutes)
BP = lga (fCt, fCin, outdir, tstar);   % function call
