%racPET_spm12_4a_MRTM2_MNI.m

% Creates parametric BPnd images for Raclopride PET data.
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
    petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,scan);
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    regDIR=fullfile(t1DIR,'registration');
    niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
    cd(petDIR)
    if ~isempty(dir(fullfile(niiDIR,'r2mm_FBP_RACd*.nii')))
        frames=dir(fullfile(niiDIR,'r2mm_FBP_RACd*.nii'));
    else
        warning('No T1 coregistered r2mm_ images found. Exiting...')
        return
    end
    
    % Run MRTM2
    first_frameIN=fullfile(niiDIR,frames(1).name);
    cerebellum_tac=fullfile(petDIR,'roi_TAC_mrtm/cerebellum_tac.txt');
    outSubject=fullfile(petDIR,'mrtm2_images');
    if ~exist(outSubject,'dir')
        mkdir(outSubject)
    end
    %usage
    [BP, k2r, R1, k2, k2a] = mrtm2 (first_frameIN, 'dec', cerebellum_tac, outSubject, thalf, 'unweighted');
    
    if ~exist(outDIR,'dir')
        mkdir(outDIR)
    end
    % Copy the BP iamge to the group directory
    subjBP=dir(fullfile(outSubject,'r2mm_*MRTM2*BP.nii'));
    copyfile([subjBP(1).folder '/' subjBP(1).name], [outDIR '/' scan '_BP_' subjDIRS(i).name '_MRTM2.nii'])
end    