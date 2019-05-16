% racPET_spm12_1a_stepCoreg.m
%
% Fix bad coregistration of PET timeframes (particularly late frames) by
% coregistering to the nearest properly coregistered frame.
%
% Data should follow the connectvitity pipeline layout, e.g.
%   a PET directory should be created in the subject directory that
%   contains a link (e.g. datadir) to the dicom PET data
%   This code expects that mCT2nii has been ran on the data, which means
%   Nii-FBP_RACd (dynamic data scans) exists and contains nifti PET frames.
%
%Contributons:
%  Evgeny Chumin, Indiana University School of Medicine, 2019
%
% Additional Notes:
%   This code was written with spm12 and fsl5.0.10 (eddy patched)
%   Required MRIread and MRIwrite functions available from the 
%   toolbox_matlab_nifti.
%-------------------------------------------------------------------------%
addpath(genpath('/usr/local/spm12')) % set path to spm12
%-------------------------------------------------------------------------%
% location of subject directories
dataDIR='/datay2/chumin-F31/data/CNT/PRISMA'; 
%-------------------------------------------------------------------------%
% subject ID to be ran (one at a time)
subj='PPE20125';
scan='PET1'; 
%-------------------------------------------------------------------------%
% frames to step coregister
fstart=52; fstop=55;
%-------------------------------------------------------------------------%

% Finding the subject
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
idx=find(strcmp({subjDIRS.name},subj)==1);
subjDIRS=subjDIRS(idx);

for i=1:length(subjDIRS)
    % set subdirectory names
    petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,scan);
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    % if dynamic data exists
    if ~isempty(dir(fullfile(petDIR,'Nii-*FBP_RACd*')))
       niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
       frames=dir(fullfile(niiDIR,'FBP*.nii'));
       if exist(niiDIR,'dir') 
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        for k=fstart:fstop
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {sprintf('%s/r%s,1',niiDIR,frames(k-1).name)};
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sprintf('%s/%s,1',niiDIR,frames(k).name)};
            spm_jobman('run',matlabbatch);
            spm_print(fullfile(niiDIR,'reg2mean_qc.ps'));
            fprintf('%d/%d - %s completed!\n',k,length(frames),frames(k).name)
            clear matlabbatch{1}.spm.spatial.coreg.estwrite.source
            clear matlabbatch{1}.spm.spatial.coreg.estwrite.ref
        end 
        clear matlabbatch
       else
           warning('Step coreg is to be used to fix frames where coregistration to mean failed.')
       end
    else
        disp(petDIR)
        warning('check that mct2nii was ran.')
    end
end