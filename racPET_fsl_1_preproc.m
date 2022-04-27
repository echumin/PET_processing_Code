% racPET_fsl_1_preproc.m
% 
% Raclopride PET image pre-processing batch script, FSL based version
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
%  Evgeny Chumin, Indiana University, Bloomington, 2022 
%  Mario Dzemidzic, Indiana University School of Medicine, 2022
% 
% Additional Notes:
%   This code was written with fsl6.0.5.1 
% 
%-------------------------------------------------------------------------%
%% Set path to your FSL and PET_Processing_Code directories.
fslpath='/N/soft/rhel7/fsl/6.0.5'; %DO NOT PUT A SLASH ON THE END
fsldatapath= strcat(fslpath,'/data'); % default on non-Gentoo fsl
%fsldatapath='/usr/share/fsl/data'; % Gentoo fsl data location
%-------------------------------------------------------------------------%
%% Set location of the subject directories.
dataDIR='/N/project/HCPaging/yoderBP_project/datadir_test'; 
%-------------------------------------------------------------------------%
%% Preprocessing is divided into two sections: preprocA and preprocB.
%   - Set the flags to 1 to perform their respective processing.
%
%   - preprocA - generate mean PET, coregistration to mean, and final 
%                realignment to mean.
%   - preprocB - brain mask PET (optional), coregister PET2 to PET1, 
%                coregister all to T1.
preprocA = 0;
preprocB = 1;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%% Subject list selection.
% Run all subjects:
    subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
% Run a single or set of subjects:
 %  subjDIRS=dir([dataDIR '/*2085']);
% For the above specified subjects, run PET 1, 2, or all
    pRUN = []; % options =1, =2, or =[] to run all PET scans.
%-------------------------------------------------------------------------%

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
% Check that T1_A processing was done
if ~isempty(petList)
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
if exist(fullfile(t1DIR,'T1_fov_denoised.nii'),'file')
    fprintf('Processing subject: %s\n',subjDIRS(i).name)
if preprocA == 1
    for p=1:length(petList) % loop over PET scans
        petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(p).name);
        % if dynamic data exists
        rawDIR = dir(fullfile(petDIR,'Nii-*FBP_RACd*'));
        if ~isempty(rawDIR)
            if size(rawDIR,1)>1
                fprintf(2,'More than 1 dynamic data directory exits. Exiting...\n')
                return
            end
            dynDIR = fullfile(rawDIR(1).folder,rawDIR(1).name);
            fprintf('%s: Processing %s\n',subjDIRS(i).name,petList(p).name)
            niiDIR=fullfile(petDIR,'fsl_dynamic_preproc');
            if ~exist(niiDIR,'dir')
                mkdir(niiDIR)
            end
%%%%%%% Create an early mean image from frames 2-15 (RACLOPRIDE DATA SPECIFIC)
        fprintf('Copy nifti frames to processing directory:\n')
        frames=dir(fullfile(dynDIR,'*.nii'));
            for ff=2:15
                if ff==2
                    filestr = fullfile(frames(ff).folder,frames(ff).name);
                else
                    filestr = [filestr ' ' fullfile(frames(ff).folder,frames(ff).name)];
                end
            end
            meanframes = fullfile(niiDIR,'frames2_15_RACd'); 
            system(sprintf('fslmerge -t %s %s',meanframes,filestr))
        
            cd(petDIR)
            fprintf('Generating an early-Mean PET for %s\n',petList(p).name)
            means{p,1} = fullfile(niiDIR,'earlyMean_f2_15_rac.nii.gz');
            system(sprintf('mcflirt -in %s -cost normmi -spline_final -meanvol',meanframes))
            movefile([meanframes '_mcf_mean_reg.nii.gz'],means{p,1})
            
%% Motion correct each frame to the early-Mean PET (Coregistration)
            fr = length(frames);
            mocoData = fullfile(niiDIR,'rFBP_RACd.nii.gz');
            tmpDIR=fullfile(niiDIR,'tmp-moco');
            tmpout = fullfile(tmpDIR,'rFrame_');
            movDIR=fullfile(niiDIR,'affine_mats_flirt');
            if ~exist(movDIR,'dir')
                mkdir(movDIR)
            end
            if ~exist(tmpDIR,'dir')
                mkdir(tmpDIR)
            end
            for ii=1:fr
                filein = fullfile(frames(ii).folder, frames(ii).name);
                filemat = fullfile(movDIR,['movepar_flirt_frame' num2str(ii) '.mat']);
                fileout = [tmpout num2str(ii) '.nii.gz'];
                % register each frame to mean
                system(sprintf('flirt -in %s -ref %s -out %s -omat %s -cost normmi -searchrx -25 25 -searchry -25 25 -searchrz -25 25 -dof 6 -interp spline',...
                    filein,means{p,1},fileout,filemat))
            end
            % concatenate the registered frames
            for ff=1:fr
                if ff==1
                    filestr = [tmpout num2str(ff) '.nii.gz'];
                else
                    filestr = [filestr ' ' tmpout num2str(ff) '.nii.gz'];
                end
            end
            system(sprintf('fslmerge -t %s %s %s',mocoData,filestr))
            if exist(mocoData,'file')        
                system(['rm -fr ' tmpDIR])
            else
                fprintf(2,'MoCo output file not created by fslmerge\n')
                return
            end 
        else
        disp(petDIR)
        warning('check that mct2nii was ran.')
        end
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
    nii1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'fsl_dynamic_preproc');
    nii2DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(2).name,'fsl_dynamic_preproc');
    if ~exist('means','var')
        mean=dir(fullfile(nii1DIR,'earlyMean*.nii.gz'));
        mean1 = mean(1).name;
        means{1,1}=fullfile(mean(1).folder,mean(1).name);
        mean=dir(fullfile(nii2DIR,'earlyMean*.nii.gz'));
        mean2=mean(1).name;
        means{2,1}=fullfile(mean(1).folder,mean(1).name);
        clear mean
    end
%     
    disp('Coregistering PET2 to PET1')
    mean1 = extractAfter(means{1,1},[nii1DIR '/']);
    mean2 = extractAfter(means{2,1},[nii2DIR '/']);
    fileout = fullfile(nii2DIR,['r_' mean2]); 
    filemat = fullfile(nii2DIR,'pet2flirt_pet1.mat');
    system(sprintf('flirt -in %s -ref %s -out %s -omat %s -cost normcorr -searchrx -25 25 -searchry -25 25 -searchrz -25 25 -dof 6 -interp spline',...
        means{2,1},means{1,1},fileout,filemat))
else
    nii1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'fsl_dynamic_preproc');
    if ~exist('means','var')   
        mean=dir(fullfile(nii1DIR,'earlyMean*.nii.gz'));
        mean1 = mean(1).name;
        means{1,1}=fullfile(mean(1).folder,mean(1).name);
        clear mean
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
 
    % register mean1 to T1 -> 
    % 2 step registration with striatal weighting on second pass.
    disp('Coregistering PET1 to T1')
    
    refweight=fullfile(nii1DIR,'striatum_weighting.nii.gz');
    % Generate T1 striatum weighting mask
    [~,result] = system(sprintf('run_first_all -s L_Accu,L_Caud,L_Pall,L_Puta,R_Accu,R_Caud,R_Pall,R_Puta -i %s.nii.gz -o striatum',T1out));
    [~,~] = system('rm striatum*bvars striatum*vtk striatum*com* striatum*org*');
    [~,~] = system('rm -fr striatum.logs');
    [~,result] = system('fslmaths striatum_all_none_firstseg.nii.gz -bin -dilD -fillh striatum_dil.nii.gz');
    [~,~] = system(sprintf('fslmaths striatum_dil.nii.gz -add 1 %s',refweight));
    
    fileref = [T1out '.nii.gz'];
    fileout = fullfile(nii1DIR,['rT1_' mean1]); 
    filemat = fullfile(nii1DIR,'pet1_to_T1.mat');
    % first pass
    system(sprintf('flirt -in %s -ref %s -out %s -omat %s -cost mutualinfo -dof 6 -interp spline -searchrx -45 45 -searchry -45 45 -searchrz -45 45',...
        means{1,1},fileref,fileout,filemat))
    %
    fileout2 = fullfile(nii1DIR,['wrT1_' mean1]); 
    filemat2 = fullfile(nii1DIR,'pet1_to_T1_2.mat'); 
    % second pass
    system(sprintf('flirt -in %s -ref %s -out %s -omat %s -refweight %s -cost mutualinfo -dof 6 -interp spline -searchrx -15 15 -searchry -15 15 -searchrz -15 15',...
        fileout,fileref,fileout2,filemat2,refweight))
    %apply transform to pet 1 data
    filein = fullfile(nii1DIR,'rFBP_RACd.nii.gz');
   
   fileout = fullfile(nii1DIR,'rT1_rFBP_RACd.nii.gz');
   system(sprintf('applyxfm4D %s %s %s %s --singlematrix',filein,fileref,fileout,filemat))
    
    if length(petList)>1
        % concatenate transformations from pet2topet1 and pet1toT1
        matout = fullfile(nii2DIR,'pet2_to_T1.mat');
        mat1 = fullfile(nii2DIR,'pet2flirt_pet1.mat');
        system(sprintf('convert_xfm -omat %s -concat %s %s',matout,mat1,filemat))
        
        % apply transformation to mean 2 
        fileout = fullfile(nii2DIR,['rT1_' mean2]); 
        system(sprintf('flirt -in %s -ref %s -out %s -init %s -applyxfm -interp spline',...
                   means{2,1},fileref,fileout,matout))
        
        %apply transfrom to pet2 data
        filein = fullfile(nii2DIR,'rFBP_RACd.nii.gz');
        fileout = fullfile(nii2DIR,'rT1_rFBP_RACd.nii.gz');
        system(sprintf('applyxfm4D %s %s %s %s --singlematrix',filein,fileref,fileout,matout))
    end   
end
clear niiDIR nii1DIR nii2DIR means petList mean1 mean2
end
