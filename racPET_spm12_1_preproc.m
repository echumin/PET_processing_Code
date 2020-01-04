% racPET_spm12_1_preproc.m
%
% Raclopride PET image processing batch script.
%
% Data should follow the connectvitity pipeline layout, e.g.
%   a PET directory should be created in the subject directory that
%   contains a link (e.g. datadir) to the dicom PET data
%   This code expects that mCT2nii has been ran on the data, which means
%   Nii-FBP_RACd (dynamic data scans) exists and contains nifti PET frames.
%
% T1 A and B preprocessing through the IUSM-connectivity-pipeline should be
% completed prior to starting these scripts.
%
% For single scan design there should only be a single 'PET' directory.
% For 2 scan designs, where scan2 needs to be aligned to scan1, there
% should be a 'PET1' and 'PET2' directories where:
%   PET1 - a baseline scan / first scan / to which the other scan is
%          registered.
%   PET2 - a challenge scan / second scan / the one that is moved to match
%          PET1
%
%Contributons:
%  Evgeny Chumin, Indiana University School of Medicine, 2019
%  Mario Dzemidzic, Indiana University School of Medicine, 2019
%
% Additional Notes:
%   This code was written with spm12 and fsl5.0.10 (eddy patched)
%   Required MRIread and MRIwrite functions available from the
%   toolbox_matlab_nifti.
%
%-------------------------------------------------------------------------%
%% Set path to your SPM directory.
%addpath(genpath('/usr/local/spm12')) % set path to spm12
addpath(genpath('/home/echumin/Documents/MATLAB/spm12'))
%-------------------------------------------------------------------------%
%% Set location of the subject directories.
dataDIR='/projects/pet_processing/datadir'; 
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
%% Create an early mean PET image
if ~isempty(petList)
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
if exist(fullfile(t1DIR,'T1_fov_denoised.nii'),'file')
    for p=1:length(petList) % loop over PET scans
        petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(p).name);
        % if dynamic data exists
        if ~isempty(dir(fullfile(petDIR,'Nii-*FBP_RACd*')))
            niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
            if ~exist(niiDIR,'dir')
                mkdir(niiDIR)
            end
           % Create an early mean image from frames 2-15 (RACLOPRIDE DATA SPECIFIC)
           copyfile(fullfile(petDIR,'Nii-*FBP_RACd*/*.nii'), niiDIR) 
           frames=dir(fullfile(niiDIR,'FBP*.nii'));
           cd(petDIR)
           fprintf('Generating an early-Mean PET for %s\n',petList(p).name)
           for j=2:15
               copyfile(fullfile(niiDIR,frames(j).name),[niiDIR '/tmp_' frames(j).name])
               matlabbatch{1}.spm.spatial.realign.estwrite.data{1,1}{j-1,1}=sprintf('%s/tmp_%s,1',niiDIR,frames(j).name);
           end
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 7;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            spm_jobman('run',matlabbatch);
            clear matlabbatch
            delete(fullfile(niiDIR,'tmp_*'))
            tmp=dir(fullfile(niiDIR,'*tmp*'));
            movefile([niiDIR '/' tmp(1).name],[niiDIR '/mean' extractAfter(tmp(1).name,'tmp_')])
            movefile([niiDIR '/' tmp(2).name],[niiDIR '/rp_' extractAfter(tmp(2).name,'tmp_')])
        else
            disp(petDIR)
            warning('check that mct2nii was ran.')
        end
        % populate a list of mean images
        mean=dir(fullfile(niiDIR,'mean*.nii'));
        means{p,1}=[fullfile(mean(1).folder,mean(1).name) ',1'];
        clear mean
%% Motion correct each frame to the early-Mean PET (Coregistration)
        disp('Coregistering individual frames to respective mean PET')
        matlabbatch{1}.spm.spatial.coreg.estimate.ref{1,1} = means{p,1};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        for k=2:length(frames)
            if k==2
                matlabbatch{1}.spm.spatial.coreg.estimate.other = {sprintf('%s/%s,1',niiDIR,frames(1).name)};
            else
                matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
            end
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {sprintf('%s/%s,1',niiDIR,frames(k).name)};
        spm_jobman('run',matlabbatch);
        fprintf('%d/%d - %s completed!\n',k,length(frames),frames(k).name)
        clear matlabbatch{1}.spm.spatial.coreg.estwrite.source
        clear matlabbatch{1}.spm.spatial.coreg.estimate.other
        end  
        clear matlabbatch
%% Final Pass Motion Correction
        disp('Realigning individual frames to respective mean PET')
        matlabbatch{1}.spm.spatial.realign.estimate.data{1,1}{1,1} = means{p,1};
        for j=1:length(frames)
            matlabbatch{1}.spm.spatial.realign.estimate.data{1,1}{j+1,1} = sprintf('%s/%s,1',niiDIR,frames(j).name);
        end
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 1;
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 7;
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
        spm_jobman('run',matlabbatch);
        clear matlabbatch
        clear frames
    end
else
    warning('No T1_fov_denoised found. Run IUSM-connectivity-pipeline T1_A and B.')
    return
end
else
    warning('No PET directories found in %s directory. Exiting...',subjDIRS(i).name)
    return
end   
%% Coregister the scan2 mean PET images and frames to scan1
cd(fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name))
if length(petList)>1
    disp('Realigning PET2 to PET1')
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{1,1} = means{1,1};
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1,1} = means{2,1};
    nii2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name,'nii_dynamic_preproc');
    frames=dir(fullfile(nii2DIR,'FBP*.nii'));
    for j=1:length(frames)
        other2{j,1}=[nii2DIR '/' frames(j).name ',1'];
    end
    clear frames
    matlabbatch{1}.spm.spatial.coreg.estimate.other = other2;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end
%% Coregister all scans to T1_fov_denoised
    T1in=fullfile(t1DIR,'T1_fov_denoised.nii');
    T1out=fullfile(t1DIR,'T1_2mm_fov_denoised');
    sentence=sprintf('flirt -in %s -out %s -interp spline -applyisoxfm 2 -ref %s',T1in,T1out,T1in);
    [~,result]=system(sentence);
    if ~isempty(result)
        disp(result)
    end
    sentence=sprintf('gunzip %s.nii.gz',T1out);
    [~,result]=system(sentence);
    if ~isempty(result)
        disp(result)
    end

    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {sprintf('%s.nii,1',T1out)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source{1,1} = means{1,1};
    nii1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'nii_dynamic_preproc');
    frames=dir(fullfile(nii1DIR,'FBP*.nii'));
    for j=1:length(frames)
        other1{j,1}=[nii1DIR '/' frames(j).name ',1'];
    end
    clear frames
    if length(petList)>1
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = vertcat(other1,means{2:end,1},other2);
    else
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = other1; 
    end   
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r2mm_';
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end
