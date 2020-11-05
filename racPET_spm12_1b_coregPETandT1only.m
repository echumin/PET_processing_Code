% racPET_spm12_1b_coregPETandT1only.m
%
% Performs the last steps from racPET_spm12_1_preproc.m, coregistration of
% mean PET images to each other (if more than 1 PET scan is found) and 
% coregistration of PET to T1_fov_denoised anatomical image.
%
% Useful if you are re-running PET to PET registration with different
% options and don't want to start from scratch.
%
%Contributons:
%  Evgeny Chumin, Indiana University, Bloomington, 2020 
%  Mario Dzemidzic, Indiana University School of Medicine, 2019
%
%-------------------------------------------------------------------------%
%% Set path to your SPM directory.
addpath(genpath('/usr/local/spm12')) % set path to spm12
addpath(genpath('/projects/pet_processing/Jenya_temp/PET_processing_Code'))
%addpath(genpath('/home/echumin/Documents/MATLAB/spm12'))
%-------------------------------------------------------------------------%
%% Set location of the subject directories.
dataDIR='/projects/pet_processing/Jenya_temp/datadir';
%-------------------------------------------------------------------------%
SkullStripMeanPET = 1; % 1 to use T1-derived brain mask

%%
% Looping across subjects
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
%subjDIRS=dir([dataDIR '/*95']); this was to run a specific subject
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
        niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
        mean=dir(fullfile(niiDIR,'mean*.nii'));
        means{p,1}=[fullfile(mean(1).folder,mean(1).name) ',1'];
        clear mean
     end
%% Coregister the scan2 mean PET images and frames to scan1
cd(fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name))
if length(petList)>1
    if SkullStripMeanPET == 1
        disp('Mean PET images will be skull stripped prior to coregistration')
        spm_figure('GetWin','Graphics');
        [masked_means] = f_maskPET(means);
        matlabbatch{1}.spm.spatial.coreg.estimate.ref{1,1} = masked_means{1,1};
        matlabbatch{1}.spm.spatial.coreg.estimate.source{1,1} = masked_means{2,1};
    else
        matlabbatch{1}.spm.spatial.coreg.estimate.ref{1,1} = means{1,1};
        matlabbatch{1}.spm.spatial.coreg.estimate.source{1,1} = means{2,1};
    end
    nii2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name,'nii_dynamic_preproc');
    frames=dir(fullfile(nii2DIR,'FBP*.nii'));
    for j=1:length(frames)
        other2{j,1}=[nii2DIR '/' frames(j).name ',1'];
    end
    clear frames
    if SkullStripMeanPET == 1
        other2{end+1,1}=means{2,1};
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.other = other2;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    disp('Realigning PET2 to PET1')
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    ps_rename('spm_coregPET2toPET1.ps')
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
    nii1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'nii_dynamic_preproc');
    frames=dir(fullfile(nii1DIR,'FBP*.nii'));
    for j=1:length(frames)
        other1{j,1}=[nii1DIR '/' frames(j).name ',1'];
    end
    clear frames
    if length(petList)>1 && SkullStripMeanPET == 1
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = vertcat(other1,other2);
    elseif length(petList)>1 && SkullStripMeanPET ~= 1
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
    ps_rename('spm_coregAll2T1.ps')
else
    warning('No T1_fov_denoised found. Run IUSM-connectivity-pipeline T1_A and B.')
    return
end
else
    warning('No PET directories found in %s directory. Exiting...',subjDIRS(i).name)
    return
end   
end