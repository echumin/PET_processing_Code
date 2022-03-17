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
% IDEAL USE CASE:
%  -You have steps 1, 2, and 3 of processing and discovered that there are
%   frames that poorly align with the MRTM fit.
%  -If this apears only in one or some of the frames of the scans, this is likely an
%   issue of poor frame coregistration to the mean PET. To fix it, one
%   needs to coregister the 'bad' frame to an adjacent 'good' frame.
%  -Unfortunately, this mean you must start processing from the beginning.
%   Run 1_preproc again and check the registration of all frames to the mean.
%   If poor realignment is the culprit, input them to this script,
%   coregister bad frames to preceeding good frames (this will update the
%   header transformation matrix), then proceed to step 2 and beyond.
%
%Contributons:
%  Evgeny Chumin, Indiana University School of Medicine, 2019
%
%-------------------------------------------------------------------------%
addpath(genpath('/geode2/soft/hps/rhel7/spm/12/')) % set path to spm12
%-------------------------------------------------------------------------%
% location of subject directories
dataDIR='/N/project/HCPaging/yoderBP_project/BP_DTI_jenya_raw'; 
%-------------------------------------------------------------------------%
% subject ID to be ran (one at a time)
subj='PPE20112';
scan='PET1'; 
%-------------------------------------------------------------------------%
% frames to step coregister
%   must be a vector of integers; must list in order smallest to largest
frames2reg = [50 51 52 53 54 55]; 
%-------------------------------------------------------------------------%
Gframe=49;

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
       if exist(niiDIR,'dir')
        frames=dir(fullfile(niiDIR,'rFBP*.nii'));
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        for k=1:length(frames2reg)
            Bframe=frames2reg(k);
           % Gframe=frames2reg(k)-1;
            matlabbatch{1}.spm.spatial.coreg.estimate.ref = {sprintf('%s/%s,1',niiDIR,frames(Gframe).name)};
            matlabbatch{1}.spm.spatial.coreg.estimate.source = {sprintf('%s/%s,1',niiDIR,frames(Bframe).name)};
            spm_jobman('run',matlabbatch);
            fprintf('%d/%d - %s completed!\n',Bframe,length(frames),frames(Bframe).name)
            clear matlabbatch{1}.spm.spatial.coreg.estwrite.source
            clear matlabbatch{1}.spm.spatial.coreg.estwrite.ref
            clear Bframe %Gframe
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