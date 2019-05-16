% racPET_spm12_5_deltaBP_untested.m
%
%Contributons:
%  Evgeny Chumin, Indiana University School of Medicine, 2019
%
% Additional Notes:
%   This code was written with spm12 and fsl5.0.10 (eddy patched)
%   Required MRIread and MRIwrite functions available from the
%   toolbox_matlab_nifti.
%
%-------------------------------------------------------------------------%
%% Set path to your SPM directory.
addpath(genpath('/usr/local/spm12')) % set path to spm12
%-------------------------------------------------------------------------%
%% Set location of the subject directories.
dataDIR='/datay2/chumin-F31/data/CNT/SKYRA';
outDIR='/datay2/chumin-F31/mrtm2_images';
%-------------------------------------------------------------------------%
%%
% Looping across subjects
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
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
        
    pet1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name);
    pet2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name);
    
    blPET=dir(fullfile(pet1DIR,'mrtm2_images/*_MRTM2_BP.nii'));
    chPET=dir(fullfile(pet2DIR,'mrtm2_images/*_MRTM2_BP.nii'));
    
    matlabbatch{1}.spm.util.imcalc.input{1,1} = blPET;
    matlabbatch{1}.spm.util.imcalc.input{2,1} = blPET;                            
    matlabbatch{1}.spm.util.imcalc.output = ['dBP_' subjDIRS(i).name '_MRTM2.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir{1,1} = outDIR;
    matlabbatch{1}.spm.util.imcalc.expression = '(i1-i2)/i1';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    
else
    warning('One or less PET directories found. Cannot calculate deltaBP. Exiting...')
    return
end
end



