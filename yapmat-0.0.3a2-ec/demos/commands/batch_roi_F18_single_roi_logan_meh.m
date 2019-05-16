%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/FAL_human/CSN/CON/CSN001_MR/Individual_shen_roi/ROIs/TACs/';
outdir = '/portal01/kkyoder/FAL_human/CSN/CON/CSN001_MR/Individual_shen_roi/ROIs/TACs/Logan/';
pid = 'CSN001-';
pid2 = 'CSN001';
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

scan = {'1_';'2_'}; 
LR = {'L'; 'R'};
region_names = { 'ACC_152';'ACC_161';'ACC_173';'ACC_176';'accumbens';'caudate';'caudate_tail';'dlPFC_140';'dlPFC_164';'dPFC_227';'INS_183';
    'INS_207';'INS_216';'mPFC_182';'mPFC_225';'OFC_179';'OFC_209';'PFC_200';'PFC_267';'putamen';'THA_155';'vmPFC_170';
    'vmPFC_188';'vmPFC_271';'vmPFC_274';'ACC_27';'ACC_54';'ACC_92';'dlPFC_20';'dlPFC_88';'dlPFC_101';'dPFC';'INS_15';
    'INS_82';'INS_98';'mPFC_69';'mPFC_81';'OFC_129';'OFC_134';'PFC_37';'THA_11';'vmPFC_3';'vmPFC_22';'vmPFC_65';'vmPFC_119';};
             
 tstar = 25;
             
for ii = 1:length (scan)
  for jj = 1:length (LR)  
    for kk = 1:length (region_names)
    region_file = strcat scann{ii}, pid, '-', LR{jj}, region_names{kk}, '_roi.txt');
   %region_file = strcat (pid, '-', region_names{jj}, '.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (scanpid, '-', 'CER.txt'))
    %Cr = fullfile (datadir, strcat (pid, '-', LR{ii}, 'CEREB_roi.txt'));
    %[BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    BP = lga (Ct, Cr, outdir, tstar);
%    fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
    strcat (region_names{jj}), BP);
    %strcat (LR{ii}, region_names{jj}), BP)
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end
fclose (fid);
