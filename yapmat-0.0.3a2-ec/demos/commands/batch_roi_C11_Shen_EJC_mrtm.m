%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/Raclopride_2008/PPE/1Aim_mCT/NTS_mCTAim1/PPE076_MR/Subcort_structural_ROI/ROIs/Meredith_ROI/';
outdir = '/portal01/kkyoder/Raclopride_2008/PPE/1Aim_mCT/NTS_mCTAim1/PPE076_MR/Subcort_structural_ROI/ROIs/Meredith_ROI/';
pid = 'PPE076-';
pid2 = 'PPE076';
thalf = 20.4;
%%%%%%%%%%%%%%%%%%%%%%%%%

%modifying for batch Logan kky
%20090107 modifying for batch MRTM and roi names
%20160212 modifying for Shen functional regions; MEH cerebellum

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T TOUCH STUFF BELOW % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [] = batch_roi (pid, datadir, outdir, thalf);

         
excel_file = fullfile (outdir, strcat (pid, 'batch-spreadsheet_mrtm.xls'));
fid = fopen (excel_file, 'w');
fprintf (fid, '%s\t%s\t%s\t%s\t%s\t%s', ...
         'region', 'BP', 'R1', 'k2 (1/min)', 'k2a (1/min)', 'k2r (1/min)');
              
% fclose (fid);

scan = {'1_'; '2_'};
LR = {'L_'; 'R_'};
region_names = { 'post_DCA'; 'post_DPU'; 'pre_DCA'; 'pre_DPU';'VST'; 
                  };
             
 %tstar = 25;
             
for ii = 1:length (scan)
 for jj = 1:length (LR)
  for kk = 1:length (region_names)
    region_file = strcat (scan{ii}, pid, LR{jj}, region_names{kk}, '_roi.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (scan{ii}, pid, '_', 'cereb_no_vermis_roi.txt'));
    [BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    %BP = lga (Ct, Cr, outdir, tstar);
%    fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
             strcat (scan{ii}, pid, LR{jj}, region_names{kk}), BP);
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end
 end
end
fclose (fid);


