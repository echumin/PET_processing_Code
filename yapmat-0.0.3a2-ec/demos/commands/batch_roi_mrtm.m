%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/Amit/ROI_BP_20100303/cerebellar_boxes_and_motion_correction_TACs';
outdir = '/portal01/kkyoder/Amit/ROI_BP_20100303/BPs';
pid = '4645';
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
LR = {'Lt'; 'Rt'};
region_names = { '_AverageMR_DCA'; '_AverageMR_VST'; '_DPU_trimmed'; 
                  };
             
 %tstar = 25;
             
  for ii = 1:length (LR)
  for jj = 1:length (region_names)
    region_file = strcat (pid, '-', LR{ii}, region_names{jj}, '_roi.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (pid, '-', 'cereb_roi.txt'));
    [BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, outdir, thalf, 'conventional');
    %BP = lga (Ct, Cr, outdir, tstar);
%    fid = fopen (excel_file, 'a');
    fprintf (fid, '\n%s\t%f\t%f\t%f\t%f\t%f', ...
             strcat (LR{ii}, region_names{jj}), BP);
%    fclose (fid);
%    fprintf ('\nFinished analyzing data file %s with input %s.', Ct, Cr);
  end
end
fclose (fid);
