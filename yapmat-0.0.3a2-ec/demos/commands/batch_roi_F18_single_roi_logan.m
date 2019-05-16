%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/FAL_human/CSN/CON/CSN005_MR/Individual_shen_roi/ROIs/TACs/';
outdir = '/portal01/kkyoder/FAL_human/CSN/CON/CSN005_MR/Individual_shen_roi/ROIs/TACs/Logan/';
pid = '1_CSN005-';
pid2 = 'CSN005';
thalf = 109.77;
%%%%%%%%%%%%%%%%%%%%%%%%%

%modifying for batch Logan kky
%modifying for fallypride rodent rois kky
%20160315 modifying for Shen human rois kky 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T TOUCH STUFF BELOW % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [] = batch_roi (pid, datadir, outdir, thalf);

excel_file = fullfile (outdir, strcat (pid, '_batch-spreadsheet_logan.xls'));
fid = fopen (excel_file, 'w');
fprintf (fid, '%s\t%s\t%s\t%s\t%s\t%s', ...
       'region', 'BP', 'R1', 'k2 (1/min)', 'k2a (1/min)', 'k2r (1/min)');
% fclose (fid);
%LR = {'L'; 'R'};
region_names = { 'L_ACC_152';'L_ACC_161';'L_ACC_173';'L_ACC_176';'L_accumbens';'L_caudate';'L_caudate_tail';'L_dlPFC_140';
    'L_dlPFC_164';'L_dPFC_227';'L_INS_183';'L_INS_207';'L_INS_216';'L_mPFC_182';'L_mPFC_225';
    'L_OFC_179';'L_OFC_209';'L_PFC_200';'L_PFC_267';'L_putamen';'L_THA_155';'L_vmPFC_170';'L_vmPFC_188';
    'L_vmPFC_271';'L_vmPFC_274';'R_ACC_27';'R_ACC_54';'R_ACC_92';'R_accumbens';'R_caudate';'R_caudate_tail';'R_dlPFC_20';
    'R_dlPFC_88';'R_dlPFC_101';'R_dPFC';'R_INS_15';'R_INS_82';'R_INS_98';'R_mPFC_69';'R_mPFC_81';'R_OFC_129';'R_OFC_134';
    'R_PFC_37';'R_PFC_131';'R_putamen';'R_THA_11';'R_vmPFC_3';'R_vmPFC_22';'R_vmPFC_65';'R_vmPFC_119'};
             
 tstar = 25;
             
%for ii = 1:length (LR)
  for jj = 1:length (region_names)
%region_file = strcat (pid, '-', LR{ii}, region_names{jj}, '_roi.txt');
    region_file = strcat (pid, '-', region_names{jj}, '.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (pid, '-', pid2, '-', 'cereb_roi.txt'))
    %Cr = fullfile (datadir, strcat (pid, '-', LR{ii}, 'CEREB_roi.txt'));
    %[BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    BP = lga (Ct, Cr, outdir, tstar);
%    fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
    strcat (pid, region_names{jj}), BP);
    %strcat (LR{ii}, region_names{jj}), BP)
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end
fclose (fid);
