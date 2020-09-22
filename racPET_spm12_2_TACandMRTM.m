% racPET_spm12_2_TACandMRTM.m
%
% This code follows the racPET_spm12_1_preproc code.
%
% -A named input parcellation from the connectivity pipeline is used to
%  extract ROI data for input indices named with the provided labels.
% -Cerebellum mask is generated on T1_fov_denoised with fsl and the vermis is
% removed via a spatial transformation from an atlas.
% -Time Activity Curves are generated for all ROI.
% -MRTM regional binding potentials are estimated via yapmat-0.0.3a2
%
% Additional Notes:
%   While not explicitly tested, the code should be able to handle truncated
%   data.
% WARNING: If your data are truncated, double check the time-frames in your
% TAC text files are generated correctly. The user is responsible for the
% accuracy of time-frame data and making sure sufficient data is available
% for accurate BP extimation. 
%
% Software Requirements:
%   - AFNI and FSL (paths for both must be defined in your bashrc)
%
% Contributors:
% Evgeny Chumin, Indiana University School of Medicine, 2019
%                Indiana University, Bloomington, 2020
% Mario Dzemidzic, Indiana University School of Medicine, 2019
%-------------------------------------------------------------------------%
%% set system specific paths
addpath(genpath('/usr/local/IUSM-connectivity-pipeline/secondary_analyses/PET_processing_Code/toolbox_matlab_nifti'))
addpath(genpath('/usr/local/IUSM-connectivity-pipeline/secondary_analyses/PET_processing_Code/yapmat-0.0.3a2-ec/src'))
%-------------------------------------------------------------------------%
%% set path to fsl for shape models
fslpath='/usr/bin/fsl'; %DO NOT PUT A SLASH ON THE END
% fsldatapath= strcat(fslpath,'/data'); % default on non-Gentoo fsl
fsldatapath='/usr/share/fsl/data'; % Gentoo fsl data location
%-------------------------------------------------------------------------%
%% set path to MNI cerebellar vermis template
vermisMNI='/usr/local/IUSM-connectivity-pipeline/secondary_analyses/PET_processing_Code/mawlawi_roi_code/cerebellum/vermis_bin_dil.nii.gz';
%-------------------------------------------------------------------------%
%% set data directory path
dataDIR='/data01/W2D/datadir_4';
%-------------------------------------------------------------------------%
%% set OUTPUT directory and file name
outDIR='/data01/W2D/datadir_out_4';
outFILE='mrtm_20200922_test';
%-------------------------------------------------------------------------%
%% set ROI IDs and labels
  % name given to parcellation in pipeline (parcs.plabel(#).name)
  parcs_plabel_name = 'tian_subcortical_S3';
  % example setups for roiDATA structure
  %         BILATERAL CASE         |         LIST CASE
  % roiDATA= {'L_', 'R_', 'label'; | roiDATA= {'ID', 'label';
  %           267,   25, 'pre_dca' |           267, 'L_pre_dca'
  %           264,   52, 'post_dca'|           52,  'R_post_dca'
  %           153,   128, 'pre_dpu'|           153, 'L_pre_dpu'
  %           277,   11, 'post_dpu'|           277, 'L_post_dpu'
  %           155,   114, 'vst'};  |           114, 'R_vst'};
  %
  % these labels currently correspond to the Tian subcortical S3
  % parcellation that is available in the IUSM-connectivity-pipeline.
  roiDATA= {'L_', 'R_', 'label';
            42,    15, 'put_rostv'
            43,    16, 'put_rostd'
            46,    19, 'caud_post'
            47,    20, 'caud_body'
            53,    11, 'caud_head'
            52,    25, 'nacc'};
%-------------------------------------------------------------------------%
%% Raclopride half-life
thalf=20.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse out output file name
outFILE=fullfile(outDIR,outFILE);
if ~exist(outDIR,'dir')
    mkdir(outDIR)
end
count=length(dir(strcat(outFILE,'*')));
if count > 0
    outFILE=fullfile(sprintf('%s_run%d.mat',outFILE,count+1));
else
    outFILE=fullfile(sprintf('%s.mat',outFILE));
end
%% Preallocate output data structures
mrtmOUTbp={'BP'}; mrtmOUTr1={'R1'}; mrtmOUTk2={'k2'}; mrtmOUTk2a={'k2a'}; mrtmOUTk2r={'k2r'};
% find number of roi (k) and whether they are bilateral or single list (j)
ub = size(roiDATA,2)-1;
lp = ub+1; % label position
if ub > 2
    fprintf(2,'Too many columns in roiDATA cell structure. Max 3 colums allowed.\n')
    return
elseif ub < 1
    fprintf(2,'roiDATA must be a minimum of 2 columns (IDs and Labels).\n')
    return
end
nr = size(roiDATA,1);
if nr < 2
    fprintf(2,'roiDATA must be a minumum of 2 rows (Column names and at least 1 ROI).\n')
    return
end
switch ub
    case 1 % in list case, ignore IDs column label
        for k=2:nr % each roi is an output row
            mrtmOUTbp{1,end+1}=roiDATA{k,lp}; %#ok<*SAGROW>
            mrtmOUTr1{1,end+1}=roiDATA{k,lp};
            mrtmOUTk2{1,end+1}=roiDATA{k,lp};
            mrtmOUTk2a{1,end+1}=roiDATA{k,lp};
            mrtmOUTk2r{1,end+1}=roiDATA{k,lp};
        end
    case 2
        for j=1:ub % in bilateral case, prepend hemisphere label
            for k=2:nr % each roi is an output row
                mrtmOUTbp{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,lp});
                mrtmOUTr1{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,lp});
                mrtmOUTk2{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,lp});
                mrtmOUTk2a{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,lp});
                mrtmOUTk2r{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,lp});
            end
        end
end

%% Loop accross subjects        
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
% subjDIRS=dir([dataDIR '/*95']); % this was to run a specific subject
for i=1:length(subjDIRS)
    disp('%---------------------------------%')
    fprintf('Setting paths to %s data .....\n',subjDIRS(i).name)
    % set PET subdirectory names
    dircont=dir(fullfile(subjDIRS(i).folder,subjDIRS(i).name)); dircont(1:2)=[];
    petList=struct.empty;
    for p=1:length(dircont)
        if dircont(p).isdir==1 && ~isempty(strfind(dircont(p).name,'PET'))
            petList(end+1).name=dircont(p).name;
        end
    end
    % Perform file checks
    if ~isempty(petList)
        t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    else
        fprintf(2,'No PET directories found. Exiting...\n')
        return
    end
    if ~exist(fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'T1_2mm_fov_denoised.nii'),'file')
        fprintf(2,'T1_2mm_fov_denoised.nii not found. Make sure preproc was ran. Exiting...\n')
        return
    end
    parcName=['T1_GM_parc_' parcs_plabel_name '.nii.gz'];
    parcFile=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,parcName);
    if ~exist(parcFile,'file')
        parc1mm=fullfile(t1DIR,parcName);
        if exist(parc1mm,'file')
            sentence=sprintf('flirt -in %s -out %s -interp nearestneighbour -applyisoxfm 2 -ref %s',parc1mm,parcFile,parc1mm);
            [~,result]=system(sentence);disp(result)
        else
            fprintf(2,'%s not found. Make sure T1_B was ran. Exiting...\n',parcName)
            return
        end
    end
%% Loop across PET sessions
for p=1:length(petList) % loop over PET scans
    petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(p).name);
    disp(petList(p).name)
    % if dynamic data exists
    if ~isempty(dir(fullfile(petDIR,'Nii-*FBP_RACd*')))
        niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
    else
        fprintf(2,'nii_dymanic_preproc not found. Make sure preproc was ran. Exiting...\n')
        return
    end
    %create a 4D pet dataset
    disp('%---------------------------------%')
    fprintf('Merging %s data into a 4D volume ....\n',subjDIRS(i).name)
    disp('%---------------------------------%')
    pet4D=fullfile(petDIR,'r2mm_T1_FBP_RACd.nii.gz');
    sentence=sprintf('fslmerge -t %s %s',pet4D,fullfile(niiDIR,'r2mm_FBP_RACd*.nii'));
    [~,result]=system(sentence); disp(result)
    frames=dir(fullfile(niiDIR,'FBP_RACd*.nii'));
    numTimePoints=length(frames);
    volNAT=MRIread(fullfile(frames(1).folder,frames(1).name)); volNAT=volNAT.vol;
    numSlices=size(volNAT,3);                                               % find number of slices in each volume
    
    % initialize time-frame data
    TimeFrames(1,1)=0;                                                      % initial TAC point starts at 0
    %referring to dicom data, generate time-frames
    dcmString=dir(fullfile(petDIR,'Nii-*FBP_RACd*/*.dcm'));                 % get list of dicoms
    dcmPath=dcmString(1).folder;                                            % get full path to dicoms
    for f=1:numTimePoints                                                   % for every timepoint
        d=1+(numSlices*(f-1));                                              % find corresponding first slice in volume
        dcmFile=fullfile(dcmPath,dcmString(d).name);                        % build dicom filename for that slice
        % get frame duration (millisecond)
        [~,frameLength]=system(sprintf('dicom_hinfo -tag 0018,1242 %s',dcmFile));
        frameLength=str2double(extractAfter(frameLength,'.dcm '))/1000;     % convert to double in seconds
        TimeFrames(f,2)=TimeFrames(f,1)+frameLength;                        % get end time for timepoint
        if f<numTimePoints                                                  % if this is not the last timepoint
            TimeFrames(f+1,1)=TimeFrames(f,2);                              % set start of next as end time of current
        end
    end
        
    if exist(parcFile,'file') && exist(niiDIR,'dir')
        roiDIR=fullfile(petDIR,'roi_TAC_mrtm');
        if ~exist(roiDIR,'dir')
            mkdir(roiDIR)
        end
        % read in PET data
        fprintf('Working on subject: %s\n',subjDIRS(i).name)
        volPET = MRIread(pet4D);  volPET=volPET.vol;
        [sizeX,sizeY,sizeZ,numTimePoints]=size(volPET);
        % Looping across ROI
        for j=1:ub
            for k=2:nr
                switch ub
                    case 1 % list case
                        rl = roiDATA{k,lp}; % roi label
                    case 2 % bilateral case
                        rl = [roiDATA{1,j} roiDATA{k,lp}];
                end
                fprintf('Starting processing for %s; Parcellation ID: %d\n',rl,roiDATA{k,j})
                disp('Extracting ROI from GM_parc')
                % extract ROI from parcellation
                roiOUT=fullfile(roiDIR,[rl '.nii.gz']);
                sentence=sprintf('fslmaths %s -thr %d -uthr %d %s',parcFile,roiDATA{k,j},roiDATA{k,j},roiOUT);
                [~,result]=system(sentence);disp(result)
                fprintf('Extracting TAC for %s\n\n',rl)
                volROI = MRIread(roiOUT); volROI=volROI.vol;
                % for each PET frame, get mean activity from mask
                for timePoint=1:numTimePoints
                    aux=reshape(volPET(:,:,:,timePoint),[sizeX,sizeY,sizeZ]); 
                    TAC(timePoint,1)=mean(aux(volROI>0));
                    clear aux
                end
                % write out time activity text file
                tac_txt=horzcat(TimeFrames,TAC);
                dlmwrite(fullfile(roiDIR,sprintf('%s_tac.txt',rl)),tac_txt,'delimiter','\t','precision','%.6f')
                clear TAC tac_txt rl
            end
        end
        
        % Creating a cerebellum mask
        disp('%---------------------------------%')
        fprintf('Generating a Cerebellar Mask\n')
        disp('%---------------------------------%')
        crblmDIR=fullfile(roiDIR,'cerebellum');
        if ~exist(crblmDIR,'dir')
            mkdir(crblmDIR)
        end
        T1=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'T1_2mm_fov_denoised.nii');
        outBasename=fullfile(crblmDIR,'subj_2_std_subc');
        [~,result]=system(sprintf('first_flirt %s %s -cort',T1,outBasename));disp(result)
        [~,~]=system(sprintf('rm %s/*subc.mat* %s/*subc.nii* %s/*cort.nii*',crblmDIR,crblmDIR,crblmDIR));
        mat=fullfile(crblmDIR,'subj_2_std_subc_cort.mat');
        for c=1:2
            if c==1
                cereb=fullfile(crblmDIR,'R_cerebellum'); side='R';
            elseif c==2
                cereb=fullfile(crblmDIR,'L_cerebellum'); side='L';
            end
            [~,result]=system(sprintf('run_first -i %s -t %s -o %s -n 40 -m %s/first/models_336_bin/intref_puta/%s_Cereb.bmv -intref %s/first/models_336_bin/05mm/%s_Puta_05mm.bmv',...
                T1,mat,cereb,fsldatapath,side,fsldatapath,side));disp(result)
            [~,result]=system(sprintf('first_boundary_corr -s %s.nii.gz -i %s -b fast -o %s_corr',cereb,T1,cereb));disp(result)
            [~,result]=system(sprintf('fslmaths %s_corr.nii.gz -bin %s_corr_bin',cereb,cereb));disp(result)
        end
        
        [~,result]=system(sprintf('fslmaths %s/R_cerebellum_corr_bin.nii.gz -add %s/L_cerebellum_corr_bin.nii.gz %s/cerebellum_bin',crblmDIR,crblmDIR,crblmDIR));disp(result)
        [~,result]=system(sprintf('fslmaths %s/cerebellum_bin.nii.gz -fillh %s/cerebellum_bin_filled',crblmDIR,crblmDIR));disp(result)
        [~,result]=system(sprintf('fslmaths %s/cerebellum_bin_filled.nii.gz -dilD %s/cerebellum_bin_filled_dil',crblmDIR,crblmDIR));disp(result)
        [~,test]=system(sprintf('rm %s/R_cerebellum* %s/L_cerebellum*',crblmDIR,crblmDIR));
        
        % Use transformations from connectivity pipeline to bring vermis mask to T1-space.
        invwarp=fullfile(t1DIR,'registration','MNI2T1_warp.nii.gz');
        invomat12=fullfile(t1DIR,'registration','MNI2T1_dof12.mat');
        invomat6=fullfile(t1DIR,'registration','MNI2T1_dof6.mat');
        invconcmat=fullfile(t1DIR,'registration','MNI2T1_linear.mat');
        
        vermisNATIVE=fullfile(crblmDIR,'vermis');
        [~,result]=system(sprintf('convert_xfm -omat %s -concat %s %s',invconcmat,invomat6,invomat12));disp(result)
        [~,result]=system(sprintf('applywarp --in=%s --out=%s --ref=%s --warp=%s --postmat=%s',...
            vermisMNI,vermisNATIVE,T1,invwarp,invconcmat));disp(result)
        % use inverse vermins mask to remove it from cerebellum mask
        [~,result]=system(sprintf('fslmaths %s.nii.gz -binv %s_binv',vermisNATIVE,vermisNATIVE));disp(result)
        crblmFull=fullfile(crblmDIR,'cerebellum_noVermis_full.nii.gz');
        [~,result]=system(sprintf('fslmaths %s/cerebellum_bin_filled.nii.gz -mas %s_binv.nii.gz %s',crblmDIR,vermisNATIVE,crblmFull));disp(result)
        % Isolate cerebellar Gray Matter
        GM2mm=fullfile(subjDIRS(i).folder,subjDIRS(i).name,petList(1).name,'T1_2mm_GM_mask.nii.gz');
        if ~exist(GM2mm,'file')
            GM1mm=fullfile(t1DIR,'T1_GM_mask.nii.gz');
            if exist(GM1mm,'file')
                sentence=sprintf('flirt -in %s -out %s -interp nearestneighbour -applyisoxfm 2 -ref %s',GM1mm,GM2mm,GM1mm);
                [~,result]=system(sentence);disp(result) 
            end
        end
        if exist(GM2mm,'file')
            crblmFINAL=fullfile(crblmDIR,'cerebellum_noVermis.nii.gz');
            [~,result]=system(sprintf('fslmaths %s -mas %s %s',crblmFull,GM2mm,crblmFINAL));disp(result)

            fprintf('Extracting TAC for cerebellar reference.\n')
            volROI = MRIread(crblmFINAL); volROI=volROI.vol;
            for timePoint=1:numTimePoints
                aux=reshape(volPET(:,:,:,timePoint),[sizeX,sizeY,sizeZ]); 
                TACcrblm(timePoint,1)=nanmean(aux(volROI>0));
                clear aux
            end       
            tac_crblm_txt=horzcat(TimeFrames,TACcrblm);
            dlmwrite(fullfile(roiDIR,'cerebellum_tac.txt'),tac_crblm_txt,'delimiter','\t','precision','%.6f')
            clear TAC tac_crblm_txt 

           % pleallocate MRTM outputs
            subjMRTMout={'ROI','BP','R1','k2','k2a','k2r'};
            mrtmOUTbp{end+1,1}=subjDIRS(i).name;
            mrtmOUTr1{end+1,1}=subjDIRS(i).name;
            mrtmOUTk2{end+1,1}=subjDIRS(i).name;
            mrtmOUTk2a{end+1,1}=subjDIRS(i).name;
            mrtmOUTk2r{end+1,1}=subjDIRS(i).name;
            % Looping over the ROI
            counter=1;
            for j=1:ub
                for k=2:nr
                    switch ub
                        case 1 % list case
                            rl = roiDATA{k,lp}; % roi label
                        case 2 % bilateral case
                            rl = [roiDATA{1,j} roiDATA{k,lp}];
                    end
                    counter=counter+1;
                    subjMRTMout{end+1,1}=sprintf('%s',rl);
                    fprintf('Running MRTM for %s\n',rl)
                    Ct=fullfile(roiDIR,sprintf('%s_tac.txt',rl));
                    Cr=fullfile(roiDIR,'cerebellum_tac.txt');
                    [BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, roiDIR, thalf, 'conventional');
                    subjMRTMout{end,2}=BP;  mrtmOUTbp{end,counter}=BP;
                    subjMRTMout{end,3}=R1;  mrtmOUTr1{end,counter}=R1;
                    subjMRTMout{end,4}=k2;  mrtmOUTk2{end,counter}=k2;
                    subjMRTMout{end,5}=k2a; mrtmOUTk2a{end,counter}=k2a;
                    subjMRTMout{end,6}=k2r; mrtmOUTk2r{end,counter}=k2r;
                    clear BP R1 k2 k2a k2r
                    close all
                end
            end
            clear counter
            % save subject results
            fileOUT=fullfile(roiDIR,sprintf('%s_mrtm_output.mat',subjDIRS(i).name));
            count=length(dir(strcat(fileOUT(1:end-4),'*')));
            if count>0
                fileOUT=fullfile(roiDIR,sprintf('%s_mrtm_output_run%d.mat',subjDIRS(i).name,count+1));
            end
            save(fileOUT,'subjMRTMout')
        else
            disp(G2mm)
            fprintf(2,'No GM mask available. Makes sure a connectivity pipeline T1_GM_mask exists.\n')
        end
    else
        disp(parcFile)
        disp(niiDIR)
        fprintf(2,'GM parcellation image above not found or/n PET nii dir not found.')
    end
end
end
% save results from the whole subjectlist
save(outFILE,'mrtmOUTbp','mrtmOUTk2','mrtmOUTr1','mrtmOUTk2a','mrtmOUTk2r')
close all
clearvars