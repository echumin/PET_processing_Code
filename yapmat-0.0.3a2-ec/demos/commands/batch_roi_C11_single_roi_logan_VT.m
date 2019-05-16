%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFY THESE VALUES  %
%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '/portal01/kkyoder/PBR-Yale/';
outdir = '/portal01/kkyoder/PBR-Yale/logan_Vt/';
pid = 'SB647RAT3_IV';
thalf = 20.4;
%%%%%%%%%%%%%%%%%%%%%%%%%

%modifying for batch Logan kky
%modifying for fallypride rodent rois kky
%20120416 modifying for PBR28 Vt

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T TOUCH STUFF BELOW % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [] = batch_roi (pid, datadir, outdir, thalf);

excel_file = fullfile (outdir, strcat (pid, '_batch-spreadsheet_logan.xls'));
fid = fopen (excel_file, 'w');
%fprintf (fid, '%s\t%s\t%s\t%s\t%s\t%s', ...
 %        'region', 'BP', 'R1', 'k2 (1/min)', 'k2a (1/min)', 'k2r (1/min)');

 fprintf (fid, '%s\t%s\', ...
         'region', 'Vt');
 
% fclose (fid);
%LR = {'L'; 'R'};
region_names = { 'BST', 'CS', 'CAU', 'CER', 'CBnV', 'CNG', 'FRN', 'GP', 'INS', 'NACC', 'OCC', 'PNS', 'PUT', 'THL';};

%20120416 will have to ascertain tstar for PBR28             
 tstar = 25;
             
%for ii = 1:length (LR)
  for jj = 1:length (region_names)
%region_file = strcat (pid, '-', LR{ii}, region_names{jj}, '_roi.txt');
    region_file = strcat (pid, '_', region_names{jj}, '.txt');
    Ct = fullfile (datadir, region_file);
    Cr = fullfile (datadir, strcat (pid, '_', '.txt'))
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
