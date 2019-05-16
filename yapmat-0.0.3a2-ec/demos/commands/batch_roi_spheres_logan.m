%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/Raclopride_2008/ABMRF_ARC/RSA_BrainPET_2009_data/RSA_BrainPET_TACs/ABM001/5480';
outdir = '/portal01/kkyoder/Raclopride_2008/ABMRF_ARC/RSA_BrainPET_2009_data/RSA_BrainPET_TACs/ABM001/5480/logan_5480';
pid = '11_5480_rest1';
thalf = 20.4;
%%%%%%%%%%%%%%%%%%%%%%%%%

%modifying for batch Logan kky

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T TOUCH STUFF BELOW % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [] = batch_roi (pid, datadir, outdir, thalf);

excel_file = fullfile (outdir, strcat (pid, '_batch-spreadsheet_logan.xls'));
fid = fopen (excel_file, 'w');
fprintf (fid, '%s\t%s\t%s\t%s\t%s\t%s', ...
         'region', 'BP', 'R1', 'k2 (1/min)', 'k2a (1/min)', 'k2r (1/min)');
% fclose (fid);
LR = {'L'; 'R'};
region_names = { '_DCA_PET_corrected'; '_DPU_PET_corrected'; '_VST_PET_corrected'; 
                 'ANACC1'; 'ANACC2'; 'ANACC3'; ...
                 'PNACC1'; 'PNACC2'; 'PNACC3'; ...
                 'APUT1';  'APUT2';  'PPUT';   ...
                 'CAUD1';  'CAUD2';  'CAUD3' };
             
 tstar = 25;
             
  for ii = 1:length (LR)
  for jj = 1:length (region_names)
    region_file = strcat (pid, '-', LR{ii}, region_names{jj}, '_roi.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (pid, '-', LR{ii}, 'CEREB_roi.txt'));
    %[BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    BP = lga (Ct, Cr, outdir, tstar);
%    fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
             strcat (LR{ii}, region_names{jj}), BP);
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end
end
fclose (fid);
