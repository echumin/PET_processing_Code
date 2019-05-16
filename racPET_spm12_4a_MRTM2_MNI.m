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
addpath(genpath('/datay2/chumin-F31/matlabscripts/toolbox_matlab_nifti'))
addpath(genpath('/datay2/chumin-F31/PET_processing_Code/yapmat-0.0.3a2-ec/src'))
addpath(genpath('/usr/local/spm12'))
%-------------------------------------------------------------------------%
    % set data directory paths
dataDIR='/datay2/chumin-F31/data/CNT/SKYRA';
outDIR='/datay2/chumin-F31/mrtm2_images';
scan='PET';
%-------------------------------------------------------------------------%
% Raclopride half-life
thalf=20.4;

% Loop accross subjects   
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
for i=1:length(subjDIRS)
    %set hardcoded paths
    disp('%---------------------------------%')
    fprintf('Setting paths to %s data .....\n',subjDIRS(i).name)
    petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'PET');
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    regDIR=fullfile(t1DIR,'registration');
    niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
    cd(petDIR)
    % normalize T1 to MNI
        % SPM12
        matlabbatch{1}.spm.tools.oldseg.data{1} = sprintf('%s/T1_fov_denoised.nii,1',t1DIR);
        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.biascor = 1;
        matlabbatch{1}.spm.tools.oldseg.output.cleanup = 2;
        matlabbatch{1}.spm.tools.oldseg.opts.tpm = {
                                                    '/usr/local/spm12/toolbox/OldSeg/grey.nii'
                                                    '/usr/local/spm12/toolbox/OldSeg/white.nii'
                                                    '/usr/local/spm12/toolbox/OldSeg/csf.nii'
                                                    };
        matlabbatch{1}.spm.tools.oldseg.opts.ngaus = [2
                                                      2
                                                      2
                                                      4];
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'mni';
        matlabbatch{1}.spm.tools.oldseg.opts.warpreg = 1;
        matlabbatch{1}.spm.tools.oldseg.opts.warpco = 25;
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.0001;
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 60;
        matlabbatch{1}.spm.tools.oldseg.opts.samp = 3;
        matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};
        matlabbatch{2}.spm.tools.oldnorm.write.subj.matname{1} = sprintf('%s/T1_fov_denoised_seg_sn.mat',t1DIR);
        matlabbatch{2}.spm.tools.oldnorm.write.subj.resample{1} = sprintf('%s/T1_fov_denoised.nii,1',t1DIR);
        matlabbatch{2}.spm.tools.oldnorm.write.roptions.preserve = 0;
        matlabbatch{2}.spm.tools.oldnorm.write.roptions.bb = [-90 -126 -72
                                                              91 91 109];
        matlabbatch{2}.spm.tools.oldnorm.write.roptions.vox = [2 2 2];
        matlabbatch{2}.spm.tools.oldnorm.write.roptions.interp = 4;
        matlabbatch{2}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
        matlabbatch{2}.spm.tools.oldnorm.write.roptions.prefix = 'w';

         spm_jobman('run',matlabbatch);
            clear matlabbatch
    
    %write out normalized PET frames (MNI)
        matlabbatch{1}.spm.tools.oldnorm.write.subj.matname{1} = sprintf('%s/T1_fov_denoised_seg_sn.mat',t1DIR);
        niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
        if ~isempty(dir(fullfile(niiDIR,'rFBP_RACd*.nii')))
            frames=dir(fullfile(niiDIR,'rFBP_RACd*.nii'));
            for j=1:length(frames)
                matlabbatch{1}.spm.tools.oldnorm.write.subj.resample{j,1} = sprintf('%s/%s,1',niiDIR,frames(j).name);
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
    

end
    
    
    
    
    
    