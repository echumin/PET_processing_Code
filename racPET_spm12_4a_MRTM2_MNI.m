%racPET_spm12_4a_MRTM2_MNI.m

% Creates parametric BPnd images for Raclopride PET data IN MNI SPACE.
%
% This code follows racPET_TACandMRTM code and assumes you have done the
% necessary quality assurance on the regional time activity curve data.

% Additional Notes:
%   While not explicitly test, the code should be able to handle truncated
%   data.
% WARNING: If your data are truncated, double check the time-frames in your
% TAC text files are generated correctly. The user is responsible for the
% accuracy of time-frame data and making sure sufficient data is available
% for accurate BP extimation. 
%
% Contributors:
% Evgeny Chumin, Indiana School of Medicine, 2019
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
    % set system specific paths
addpath(genpath('/datay2/PET_processing_Code/matlabscripts/toolbox_matlab_nifti'))
addpath(genpath('/datay2/PET_processing_Code/yapmat-0.0.3a2-ec/src'))
spm12_path='/usr/local/spm12';
addpath(genpath(spm12_path))
%-------------------------------------------------------------------------%
    % set data directory paths
dataDIR='/datay2/chumin-F31/data/CNT/SKYRA';
outDIR='/datay2/chumin-F31/mrtm2_MNI_images';
%-------------------------------------------------------------------------%
% Raclopride half-life
thalf=20.4;

%% Loop accross subjects   
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
    
    %set hardcoded paths
    disp('%---------------------------------%')
    fprintf('Setting paths to %s data .....\n',subjDIRS(i).name)
    pet1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name);
    nii1DIR=fullfile(pet1DIR,'nii_dynamic_preproc');
    cd(pet1DIR)
    % coregister estimate to template
        % SPM12
        matlabbatch{1}.spm.spatial.coreg.estimate.ref{1} = sprintf('%s/canonical/avg152T1.nii,1',spm12_path);
        matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = sprintf('%s/T1_2mm_fov_denoised.nii,1',pet1DIR);
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    % normalize T1 to MNI
        matlabbatch{2}.spm.tools.oldseg.data{1} = sprintf('%s/T1_2mm_fov_denoised.nii,1',pet1DIR);
        matlabbatch{2}.spm.tools.oldseg.output.GM = [0 0 0];
        matlabbatch{2}.spm.tools.oldseg.output.WM = [0 0 0];
        matlabbatch{2}.spm.tools.oldseg.output.CSF = [0 0 0];
        matlabbatch{2}.spm.tools.oldseg.output.biascor = 1;
        matlabbatch{2}.spm.tools.oldseg.output.cleanup = 2;
        matlabbatch{2}.spm.tools.oldseg.opts.tpm = {
                                                    sprintf('%s/toolbox/OldSeg/grey.nii',spm12_path)
                                                    sprintf('%s/toolbox/OldSeg/white.nii',spm12_path)
                                                    sprintf('%s/toolbox/OldSeg/csf.nii',spm12_path)
                                                    };
        matlabbatch{2}.spm.tools.oldseg.opts.ngaus = [2
                                                      2
                                                      2
                                                      4];
        matlabbatch{2}.spm.tools.oldseg.opts.regtype = 'mni';
        matlabbatch{2}.spm.tools.oldseg.opts.warpreg = 1;
        matlabbatch{2}.spm.tools.oldseg.opts.warpco = 25;
        matlabbatch{2}.spm.tools.oldseg.opts.biasreg = 0.0001;
        matlabbatch{2}.spm.tools.oldseg.opts.biasfwhm = 60;
        matlabbatch{2}.spm.tools.oldseg.opts.samp = 3;
        matlabbatch{2}.spm.tools.oldseg.opts.msk = {''};
        matlabbatch{3}.spm.tools.oldnorm.write.subj.matname{1} = sprintf('%s/T1_2mm_fov_denoised_seg_sn.mat',pet1DIR);
        matlabbatch{3}.spm.tools.oldnorm.write.subj.resample{1} = sprintf('%s/T1_2mm_fov_denoised.nii,1',pet1DIR);
        matlabbatch{3}.spm.tools.oldnorm.write.roptions.preserve = 0;
        matlabbatch{3}.spm.tools.oldnorm.write.roptions.bb = [-90 -126 -72
                                                              91 91 109];
        matlabbatch{3}.spm.tools.oldnorm.write.roptions.vox = [2 2 2];
        matlabbatch{3}.spm.tools.oldnorm.write.roptions.interp = 4;
        matlabbatch{3}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
        matlabbatch{3}.spm.tools.oldnorm.write.roptions.prefix = 'w';

         spm_jobman('run',matlabbatch);
            clear matlabbatch
     
    %write out normalized PET frames (MNI)
        matlabbatch{1}.spm.tools.oldnorm.write.subj.matname{1} = sprintf('%s/T1_2mm_fov_denoised_seg_sn.mat',pet1DIR);
        if length(petList) > 1
                pet2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name);
                nii2DIR=fullfile(pet2DIR,'nii_dynamic_preproc');
                if ~isempty(dir(fullfile(nii2DIR,'r2mm_FBP_RACd*.nii')))
                    frames=dir(fullfile(nii2DIR,'r2mm_FBP_RACd*.nii'));
                    for j=1:length(frames)
                        pet2frames{j,1} = sprintf('%s/%s,1',nii2DIR,frames(j).name);
                    end
                    clear frames
                end
        end
        if ~isempty(dir(fullfile(nii1DIR,'r2mm_FBP_RACd*.nii')))
            frames=dir(fullfile(nii1DIR,'r2mm_FBP_RACd*.nii'));
            for j=1:length(frames)
                pet1frames{j,1} = sprintf('%s/%s,1',nii1DIR,frames(j).name);
            end
            clear frames
            if length(petList) > 1
                matlabbatch{1}.spm.tools.oldnorm.write.subj.resample=vertcat(pet1frames,pet2frames);
            else
                matlabbatch{1}.spm.tools.oldnorm.write.subj.resample=pet1frames; 
            end
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.preserve = 0;
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.bb = [-90 -126 -72
                                                                  91 91 109];
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.vox = [2 2 2];
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.interp = 4;
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.prefix = 'w';
         
            spm_jobman('run',matlabbatch);
            clear matlabbatch
        else
            warning('check that MR registered PET frames are present. Exiting...')
            return
        end
    
    % Run MRTM2
    for p=1:length(petList)
        petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(p).name);
        niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
        frames=dir(fullfile(niiDIR,'r2mm_FBP_RACd*.nii'));
        normDynamicIN=fullfile(niiDIR,['w' frames(1).name]);
    cerebellum_tac=fullfile(petDIR,'roi_TAC_mrtm/cerebellum_tac.txt');
    outSubject=fullfile(petDIR,'mrtm2_mni_images');
    if ~exist(outSubject,'dir')
        mkdir(outSubject)
    end
    
    %usage
    [BP, k2r, R1, k2, k2a] = mrtm2 (normDynamicIN, 'dec', cerebellum_tac, outSubject, thalf, 'unweighted');
    
    if ~exist(outDIR,'dir')
        mkdir(outDIR)
    end
    % COPY THE BP IMAME TO GROUP DIRECTORY
    subjBP=dir(fullfile(outSubject,'wr2mm_*MRTM2*BP.nii'));
    copyfile([subjBP(1).folder '/' subjBP(1).name], [outDIR '/' petList(p).name '_BP_' subjDIRS(i).name '_MNI_MRTM2.nii'])
    end
end