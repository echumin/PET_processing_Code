%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/Raclopride_2008/PPE/1Aim_mCT/NTS_mCTAim1/PPE076_MR/Subcort_structural_ROI/ROIs/Meredith_ROI/Martinez_TACs/';
outdir = '/portal01/kkyoder/Raclopride_2008/PPE/1Aim_mCT/NTS_mCTAim1/PPE076_MR/Subcort_structural_ROI/ROIs/Meredith_ROI/Martinez_TACs/';
pid = 'PPE';
thalf = 20.4;
%%%%%%%%%%%%%%%%%%%%%%%%%

%modifying for batch Logan kky
%20090107 modifying for batch MRTM and roi names

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T TOUCH STUFF BELOW % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [] = batch_roi (pid, datadir, outdir, thalf);

excel_file = fullfile (outdir, strcat (pid, '_batch-spreadsheet_mrtm.xls'));
fid = fopen (excel_file, 'w');
fprintf (fid, '%s\t%s\t%s\t%s\t%s\t%s', ...
         'region', 'BP', 'R1', 'k2 (1/min)', 'k2a (1/min)', 'k2r (1/min)');
% fclose (fid);
LR = {'L_'; 'R_'};
region_names = { 'Post_DCA'; 'Post_DPU'; 'Pre_DCA'; 'Pre_DPU'; 'VST'; 
                  };
             
 %tstar = 25;
             
  for ii = 1:length (LR)
  for jj = 1:length (region_names)
    region_file = strcat (pid, '-', LR{ii}, region_names{jj}, '_roi.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (pid, '-', 'whole_cerebellum_non-smoking_controls_NPE_n-9_roi.txt'));
    [BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    %BP = lga (Ct, Cr, outdir, tstar);
%    fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
             strcat (pid, LR{ii}, region_names{jj}), BP);
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end
end
fclose (fid);
