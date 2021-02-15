%racPET_spm12_4_MRTM2.m
%
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
% Evgeny Chumin, Indiana University School of Medicine, 2019
% Mario Dzemidzic, Indiana University School of Medicine, 2019
%-------------------------------------------------------------------------%
    % set system specific paths
addpath(genpath('/usr/local/spm12')) % set path to spm12
addpath(genpath('/projects/pet_processing/PET_processing_Code'))
%-------------------------------------------------------------------------%
    % set data directory paths
dataDIR='/projects/pet_processing/datadir';
%-------------------------------------------------------------------------%
% Raclopride half-life
thalf=20.4;

%% Loop accross subjects   
%subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
subjDIRS=dir([dataDIR '/*95']);% this was to run a specific subject

for i=1:length(subjDIRS)
    % set PET subdirectory names
    dircont=dir(fullfile(subjDIRS(i).folder,subjDIRS(i).name)); dircont(1:2)=[];
    petList=struct.empty;
    for p=1:length(dircont)
        if dircont(p).isdir==1 && ~isempty(strfind(dircont(p).name,'PET'))
            petList(end+1).name=dircont(p).name;
        end
    end
for p=1:length(petList) % loop over PET scans
    %set hardcoded paths
    disp('%---------------------------------%')
    fprintf('Setting paths to %s %s data .....\n',subjDIRS(i).name,petList(p).name)
    petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(p).name);
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
    cd(petDIR)
    if ~isempty(dir(fullfile(niiDIR,'r2mm_*rFBP_RACd*.nii')))
        frames=dir(fullfile(niiDIR,'r2mm_*rFBP_RACd*.nii'));
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
end
end    