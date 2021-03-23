% racPET_spm12_5_deltaBP.m
%
% Computed delta BP (baseline-challenge)/baseline using SPM imcalc.
%
% This code follows racPET_spm12_4_MRTM2 utilizing the T1-space MRTM2 BP
% images that are written out by it.
%
% Additional Notes:
%   This code was written with spm12 and fsl5.0.10 (eddy patched)
%   Required MRIread and MRIwrite functions available from the
%   toolbox_matlab_nifti.
%
%Contributons:
%  Evgeny Chumin, Indiana University School of Medicine, 2019
%                 Indiana University, Bloomington, 2020
%
%-------------------------------------------------------------------------%
%% Set path to the SPM directory.
addpath('/usr/local/spm12') % set path to spm12
%-------------------------------------------------------------------------%
%% Set location of the subject directories.
dataDIR='/projects/pet_processing/datadir';
outDIR='/projects/pet_processing/mrtm2_dBP_native';
%-------------------------------------------------------------------------%
%% Set location of the scan order file (if empty as is order is assumed).
% This should be a 3 column tab delimited text file with subjectID, index 
%    of 1st scan, and index of 2nd scan.
orderTXT='/projects/pet_processing/order_list_test.txt';
%-------------------------------------------------------------------------%
%% Subject list selection.
% Run all subjects:
    %subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
% Run a single or set of subjects:
    subjDIRS=dir([dataDIR '/*01']);
%-------------------------------------------------------------------------%
%% End of user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if ~exist(outDIR,'dir')
    mkdir(outDIR)
end
if exist(orderTXT,'file')
    ord=readtable(orderTXT);
end
%% Looping across subjects
for i=1:length(subjDIRS)  
     % set PET subdirectory names
    dircont=dir(fullfile(subjDIRS(i).folder,subjDIRS(i).name)); dircont(1:2)=[];
    petList=struct.empty;
    for p=1:length(dircont)
        if dircont(p).isdir==1 && ~isempty(strfind(dircont(p).name,'PET'))
            petList(end+1).name=dircont(p).name;
        end
    end
%% 
if length(petList) > 1
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    if exist('ord','var')
        idx=find(strcmp(subjDIRS(i).name,ord.Var1));
        pet1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(ord.Var2(idx)).name);
        pet2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(ord.Var3(idx)).name);
    else  
        pet1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name);
        pet2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name);
    end
    
    blPET=dir(fullfile(pet1DIR,'mrtm2_images/*_MRTM2_BP.nii'));
    chPET=dir(fullfile(pet2DIR,'mrtm2_images/*_MRTM2_BP.nii'));
    
    matlabbatch{1}.spm.util.imcalc.input{1,1} = fullfile(blPET.folder,blPET.name);
    matlabbatch{1}.spm.util.imcalc.input{2,1} = fullfile(chPET.folder,chPET.name);                            
    matlabbatch{1}.spm.util.imcalc.output = ['dBP_' subjDIRS(i).name '_MRTM2.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir{1,1} = outDIR;
    matlabbatch{1}.spm.util.imcalc.expression = '(i1-i2)./i1';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
    spm_jobman('run',matlabbatch);
    clear matlabbatch   
else
    warning('One or less PET directories found. Cannot calculate deltaBP. Exiting...')
    return
end
end