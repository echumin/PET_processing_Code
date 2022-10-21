% racPET_spm12_1_preproc.m
% 
% Raclopride PET image pre-processing batch script.
% 
% Data should follow the IUSM-connectvitity-pipeline layout, e.g.
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
% Contributons:
%  Evgeny Chumin, Indiana University School of Medicine, 2019
%                 Indiana University, Bloomington, 2020 
%  Mario Dzemidzic, Indiana University School of Medicine, 2019
% 
% Additional Notes:
%   This code was written with spm12 (matlab19b) and fsl5.0.10 
%   (eddy patched) later tested with fsl6.0.1.
%   Required MRIread and MRIwrite functions available from the
%   toolbox_matlab_nifti.
% 
%-------------------------------------------------------------------------%
%% Set path to your SPM and PET_Processing_Code directories.
addpath('/usr/local/spm12') % set path to spm12
%addpath('/geode2/soft/hps/rhel7/spm/12/') % set path to spm12
%addpath(genpath('/N/project/HCPaging/yoderBP_project/PET_processing_Code'))
addpath(genpath('/data01/W2D/w2d_proc_mar22/PET_processing_Code'))
%-------------------------------------------------------------------------%
%% Set location of the subject directories.
dataDIR='/data01/W2D/w2d_proc_mar22/datadir_5'; 
%-------------------------------------------------------------------------%
%% Preprocessing is divided into two sections: preprocA and preprocB.
%   - Set the flags to 1 to perform their respective processing.
%
%   - preprocA - generate mean PET, coregistration to mean, and final 
%                realignment to mean.
%   - preprocB - brain mask PET (optional), coregister PET2 to PET1, 
%                coregister all to T1.
preprocA = 1;
preprocB = 0;
%-------------------------------------------------------------------------%
%% Brain mask options. 
% FOR REGISTRATION PURPOCES ONLY
    SkullStripMeanPET = 0; % 1 to use T1-derived brain mask, 0 for no masking
% APPLY MASKING TO ALL PET DATA
    % Works only if the above option to mask mean PET is on.
    maskALLpet = 0; % 1 to apply brain mask to all frames, 0 for no masking
%-------------------------------------------------------------------------%
%% Field-of-view editing options.
    %--% WARNING %--%
% Be very careful when using this parameter to edit field of view as you
% can crop the cerebellum if you are not carefull. This is meant to be used
% when rerunning specific subjects, where there is high noise in the 
% inferior of the image that is believed to be affecting registration.
    Edit_inferiorZslices = 0; % dangerous flag
% to remove inferior slices enter "-N" (+N to add)
%-------------------------------------------------------------------------%
%% Subject list selection.
% Run all subjects:
    subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
% Run a single or set of subjects:
  % subjDIRS=dir([dataDIR '/*2085']);
% For the above specified subjects, run PET 1, 2, or all
    pRUN = []; % options =1, =2, or =[] to run all PET scans.
%-------------------------------------------------------------------------%
%% Advanced flags.
intdpbg = 0; % intermediate debug : writes a 4D image for each step
             % should help identify potential issues.
skipcoreg = 0; % when rerunning for debugging, skip coregistration of frames
               % to mean and use existing r*nii images
%% End of user input
% Looping across subjects
for i=1:length(subjDIRS)
    % set PET subdirectory names
    dircont=dir(fullfile(subjDIRS(i).folder,subjDIRS(i).name)); dircont(1:2)=[];
    petList=struct.empty;
    for p=1:length(dircont)
        if dircont(p).isdir==1 && ~isempty(strfind(dircont(p).name,'PET'))
            petList(end+1).name=dircont(p).name;
        end
    end
    if ~isempty(pRUN)
        petList = petList(pRUN);
    end
%% Create an early mean PET image
if ~isempty(petList)
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
if exist(fullfile(t1DIR,'T1_fov_denoised.nii'),'file')
    fprintf('Processing subject: %s\n',subjDIRS(i).name)
if preprocA == 1
    for p=1:length(petList) % loop over PET scans
        petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(p).name);
        % if dynamic data exists
        if ~isempty(dir(fullfile(petDIR,'Nii-*FBP_RACd*')))
            fprintf('%s: Processing %s\n',subjDIRS(i).name,petList(p).name)
            niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
            if ~exist(niiDIR,'dir')
                mkdir(niiDIR)
            end
%%%%%%% Create an early mean image from frames 2-15 (RACLOPRIDE DATA SPECIFIC)
          % copyfile(fullfile(petDIR,'Nii-*FBP_RACd*/*.nii'), niiDIR) 
          fprintf('Copy nifti frames to processing directory:\n')
           system(sprintf('cp %s/Nii-*FBP_RACd*/*.nii %s',petDIR,niiDIR))
           frames=dir(fullfile(niiDIR,'FBP*.nii'));
           if Edit_inferiorZslices ~= 0
               fprintf('Editing Field of View of Raw Data\n')
               f_edit_Z_slices(frames,Edit_inferiorZslices)
           end
           cd(petDIR)
           fprintf('Generating an early-Mean PET for %s\n',petList(p).name)
           spm_figure('GetWin','Graphics');
           for j=2:15
               copyfile(fullfile(niiDIR,frames(j).name),[niiDIR '/tmp_' frames(j).name])
               matlabbatch{1}.spm.spatial.realign.estwrite.data{1,1}{j-1,1}=sprintf('%s/tmp_%s,1',niiDIR,frames(j).name);
           end
           if intdpbg == 1
               if ~exist([petDIR '/intdbg'],'dir')
                   mkdir([petDIR '/intdbg'])
               end
               infiles = [niiDIR '/tmp_*nii'];
               outfile = [petDIR '/intdbg/raw4_earlyMean'];
               [~,result]=system(sprintf('fslmerge -t %s %s',outfile,infiles));
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
            ps_rename('spm_earlyMean.ps')
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
        if intdpbg == 1
            copyfile(fullfile(mean(1).folder,mean(1).name),[petDIR '/intdbg/' mean(1).name])
        end
        clear mean
%% Motion correct each frame to the early-Mean PET (Coregistration)
if skipcoreg == 1
    chk = dir([niiDIR '/r*nii']);
    if ~isempty(chk)
        disp('Bypassing Coregistration: using existing r*nii images')
        clear chk
        rn = 0;
    else
        disp('Bypassing Coregistration not possible: no r*nii images')
        rn = 1;
    end
else
    rn = 1;
end
if rn == 1
        disp('Coregistering individual frames to respective mean PET')
        disp('Operating on temporary copies, original image headers remain unedited')
        spm_figure('GetWin','Graphics');
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref{1,1} = means{p,1};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        for k=2:length(frames)
            % frame 1 is carried along with frame 2, because it is unlikely
            % to have sufficient snr for registration.
            if k==2 
               copyfile(fullfile(niiDIR,frames(1).name),[niiDIR '/tmp_' frames(1).name])
               matlabbatch{1}.spm.spatial.coreg.estwrite.other = {sprintf('%s/tmp_%s,1',niiDIR,frames(1).name)};
            else
                matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
            end    
        copyfile(fullfile(niiDIR,frames(k).name),[niiDIR '/tmp_' frames(k).name])
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sprintf('%s/tmp_%s,1',niiDIR,frames(k).name)};
        spm_jobman('run',matlabbatch);
        fprintf('%d/%d - %s completed!\n',k,length(frames),frames(k).name)
        clear matlabbatch{1}.spm.spatial.coreg.estwrite.source
        clear matlabbatch{1}.spm.spatial.coreg.estwrite.other
        if k==2
            delete([niiDIR '/tmp_' frames(1).name])
            movefile([niiDIR '/rtmp_' frames(1).name],fullfile(niiDIR,['r' frames(1).name]))
        end
        delete([niiDIR '/tmp_' frames(k).name])
        movefile([niiDIR '/rtmp_' frames(k).name],fullfile(niiDIR,['r' frames(k).name]))
        end  
        clear matlabbatch
        ps_rename('spm_coreg2mean.ps')
        if intdpbg == 1
           infiles = [niiDIR '/r*nii'];
           outfile = [petDIR '/intdbg/coregistered_preRealignEst'];
           [~,result]=system(sprintf('fslmerge -t %s %s',outfile,infiles));
        end
end
%% Final Pass Motion Correction
        disp('Realigning individual frames to respective mean PET')
        spm_figure('GetWin','Graphics');
        matlabbatch{1}.spm.spatial.realign.estimate.data{1,1}{1,1} = means{p,1};
        for j=2:length(frames)
            matlabbatch{1}.spm.spatial.realign.estimate.data{1,1}{j,1} = sprintf('%s/%s,1',niiDIR,['r' frames(j).name]);
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
        ps_rename('spm_finalRealign.ps')
    end
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
if preprocB == 1
cd(fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name))
if length(petList)>1
    nii1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'nii_dynamic_preproc');
    nii2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name,'nii_dynamic_preproc');
    if ~exist('means','var')
        mean=dir(fullfile(nii1DIR,'mean*.nii'));
        means{1,1}=[fullfile(mean(1).folder,mean(1).name) ',1'];
        mean=dir(fullfile(nii2DIR,'mean*.nii'));
        means{2,1}=[fullfile(mean(1).folder,mean(1).name) ',1'];
        clear mean
    end
    disp('Coregistering PET2 to PET1')
    if SkullStripMeanPET == 1
        disp('Mean PET images will be skull stripped prior to coregistration')
        [masked_means] = f_maskPET(means,maskALLpet);
        if maskALLpet == 1
            frames=dir(fullfile(nii2DIR,'brFBP*.nii'));
        else
            frames=dir(fullfile(nii2DIR,'rFBP*.nii'));
        end
        for j=1:length(frames)
            other2{j,1}=[nii2DIR '/' frames(j).name ',1'];
        end
        clear frames
        other2{end+1,1} = means{2,1};
        f_spm_coreg(masked_means{1,1},masked_means{2,1},other2)
    else
        frames=dir(fullfile(nii2DIR,'rFBP*.nii'));
        for j=1:length(frames)
            other2{j,1}=[nii2DIR '/' frames(j).name ',1'];
        end
        clear frames
        f_spm_coreg(means{1,1},means{2,1},other2)
    end    
    ps_rename('spm_coregPET2toPET1.ps')
else
    nii1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'nii_dynamic_preproc');
    if ~exist('means','var')   
        mean=dir(fullfile(nii1DIR,'mean*.nii'));
        means{1,1}=[fullfile(mean(1).folder,mean(1).name) ',1'];
        clear mean
    end
    if SkullStripMeanPET == 1
        disp('Mean PET images will be skull stripped prior to coregistration')
        [masked_means] = f_maskPET(means,maskALLpet);
    end
end
%% Coregister all scans to T1_fov_denoised
    T1in=fullfile(t1DIR,'T1_fov_denoised.nii');
    T1out=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'T1_2mm_fov_denoised');
    sentence=sprintf('flirt -in %s -out %s -interp spline -applyisoxfm 2 -ref %s',T1in,T1out,T1in);
    [~,result]=system(sentence);
    if ~isempty(result)
        disp(result)
    end
    sentence=sprintf('gunzip -f %s.nii.gz',T1out);
    [~,result]=system(sentence);
    if ~isempty(result)
        disp(result)
    end
    
    spm_figure('GetWin','Graphics');
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {sprintf('%s.nii,1',T1out)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source{1,1} = means{1,1};
    
    if maskALLpet==1
        frames=dir(fullfile(nii1DIR,'brFBP*.nii'));
    else
        frames=dir(fullfile(nii1DIR,'rFBP*.nii'));
    end
    for j=1:length(frames)
        other1{j,1}=[nii1DIR '/' frames(j).name ',1'];
    end
    clear frames 
    if length(petList)>1 && SkullStripMeanPET == 1
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = vertcat(other1,masked_means,other2); % THIS MAY BE MISSING A SECOND FILE SPECIFICATION
    elseif length(petList)>1 && SkullStripMeanPET ~= 1
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = vertcat(other1,means{2:end,1},other2);
    else
        if maskALLpet==1
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = vertcat(other1,masked_means); 
        else
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = other1; 
        end
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
    ps_rename('spm_coregAll2T1.ps')
end
clear niiDIR nii1DIR nii2DIR means petList
end
