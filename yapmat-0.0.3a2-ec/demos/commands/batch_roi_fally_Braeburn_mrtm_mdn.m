%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/FAL_human/Braeburn_ATI9242/9242-14/9242-14_MR/9242-14_ROIs/';
%datadir = '~/kky/MDN/';
outdir = '/portal01/kkyoder/FAL_human/Braeburn_ATI9242/9242-14/9242-14_MR/MRTM/';
%outdir = '~/kky/MDN/';
pid = 'DRUG-';
%pid = 'BL-9242-06-';
thalf = 109.77;
%%%%%%%%%%%%%%%%%%%%%%%%%

%modifying for batch Logan kky
%20090107 modifying for batch MRTM and roi names
%20090226 modifying for fallypride rodent rois
%20160722 modifying for Braeburn

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T TOUCH STUFF BELOW % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [] = batch_roi (pid, datadir, outdir, thalf);

excel_file = fullfile (outdir, strcat (pid, 'batch-spreadsheet_mrtm.xls'));
fid = fopen (excel_file, 'w');
fprintf (fid, '%s\t%s\t%s\t%s\t%s\t%s', ...
         'region', 'BP', 'R1', 'k2 (1/min)', 'k2a (1/min)', 'k2r (1/min)');
% fclose (fid);
%LR = {'L'; 'R'};
region_names = {'9242-14_L_Caudate_edit_roi'; '9242-14_R_Caudate_edit_roi'; '9242-14_L_Putamen_edit_roi'; '9242-14_R_Putamen_edit_roi'; '9424-14_L_Accumbens_edit_roi'; '9242-14_R_Accumbens_edit_roi';};
             
%tstar = 25;
             
%for ii = 1:length (LR)
  for jj = 1:length (region_names)
    region_file = strcat (pid, region_names{jj}, '.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (pid, '9242-14_Cerebellum_bin_filled_imcalc_edit_roi', '.txt'));
    %Cr = fullfile (datadir, strcat (pid, '-', LR{ii}, 'CEREB_roi.txt'));
    [BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    %BP = lga (Ct, Cr, outdir, tstar);
    %fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
         strcat (region_names{jj}), BP, R1, k2, k2a, k2r);
         %strcat (LR{ii}, region_names{jj}), BP);
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end

fclose (fid);
