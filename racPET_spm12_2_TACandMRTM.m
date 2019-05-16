% racPET_TACandMRTM.m
%
% This code follows the racPET_spm12_1_preproc code.
%
% -The shen286 parcellation is used for striatal ROI extraction (the Mawlawi
% 2001 ROI).
% -Cerebellum mask is generated on T1_fov_denoised with fsl and the vermis is
% removed via a spatial transformation from an atlas.
% -Time Activity Curves are generated for all ROI.
% -MRTM regional binding potentials are estimated via yapmat-0.0.3a2
%
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
    % set system specific paths
addpath(genpath('/usr/local/matlabscripts/toolbox_matlab_nifti'))
addpath(genpath('/datay2/chumin-F31/PET_processing_Code/yapmat-0.0.3a2-ec/src'))
%-------------------------------------------------------------------------%
    % set path to MNI cerebellar vermis template
vermisMNI='/datay2/chumin-F31/PET_processing_Code/mawlawi_roi_code/cerebellum/vermis_bin_dil.nii.gz';
%-------------------------------------------------------------------------%
    % set data directory paths
dataDIR='/datay2/chumin-F31/data/CNT/SKYRA';
scan='PET';
outDIR='/datay2/chumin-F31/results';
outFILE='CNT_SKYRA_mrtm_20190516_test';
%-------------------------------------------------------------------------%
roiDATA= {'L_', 'R_', 'label';
          267,   25, 'pre_dca'
          264,   52, 'post_dca'
          153,   128, 'pre_dpu'
          277,   11, 'post_dpu'
          155,   114, 'vst'};
%-------------------------------------------------------------------------%
% Raclopride half-life
thalf=20.4;

% Preallocate output data structures
mrtmOUTbp={'BP'}; mrtmOUTr1={'R1'}; mrtmOUTk2={'k2'}; mrtmOUTk2a={'k2a'}; mrtmOUTk2r={'k2r'};
for j=1:2
    for k=2:6
        mrtmOUTbp{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,3});
        mrtmOUTr1{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,3});
        mrtmOUTk2{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,3});
        mrtmOUTk2a{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,3});
        mrtmOUTk2r{1,end+1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,3});
    end
end
% Loop accross subjects        
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
for i=1:length(subjDIRS)
    %set hardcoded paths
    disp('%---------------------------------%')
    fprintf('Setting paths to %s data .....\n',subjDIRS(i).name)
    petDIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,scan);
    t1DIR=fullfile(subjDIRS(i).folder,subjDIRS(i).name,'T1');
    parcFile=fullfile(t1DIR,'T1_GM_parc_shen_286.nii.gz');
    niiDIR=fullfile(petDIR,'nii_dynamic_preproc');
    %create a 4D pet dataset
    disp('%---------------------------------%')
    fprintf('Mergning %s data into a 4D volume ....\n',subjDIRS(i).name)
    disp('%---------------------------------%')
    pet4D=fullfile(petDIR,'rT1_FBP_RACd.nii.gz');
    sentence=sprintf('fslmerge -t %s %s',pet4D,fullfile(niiDIR,'rFBP_RACd*.nii'));
    [~,result]=system(sentence); disp(result)
    frames=dir(fullfile(niiDIR,'FBP_RACd*.nii'));
    numTimePoints=length(frames);
    volNAT=MRIread(fullfile(frames(1).folder,frames(1).name)); volNAT=volNAT.vol;
    numSlices=size(volNAT,3);
    
    % initialize time-frame data
    TimeFrames(1,1)=0;
    %referring to dicom data, generate time-frames
    dcmString=dir(fullfile(petDIR,'Nii-*FBP_RACd*/*000001*.dcm'));
    dcmPath=dcmString(1).folder;
    dcmString=strsplit(dcmString(1).name,'-');
    dcmString=strcat(dcmString{1},'-',dcmString{2},'-');
    for f=1:numTimePoints
        d=1+(numSlices*(f-1));
        if (d<10) 
            dcmFile=fullfile(dcmPath,sprintf('%s00000%d.dcm',dcmString,d));
        elseif (d>=10) && (d<100)
            dcmFile=fullfile(dcmPath,sprintf('%s0000%d.dcm',dcmString,d));
        elseif (d>=100) && (d<1000)
            dcmFile=fullfile(dcmPath,sprintf('%s000%d.dcm',dcmString,d));
        elseif (d>=1000) && (d<10000)
            dcmFile=fullfile(dcmPath,sprintf('%s00%d.dcm',dcmString,d));
        end
        [~,frameLength]=system(sprintf('dicom_hinfo -tag 0018,1242 %s',dcmFile));
        frameLength=str2double(extractAfter(frameLength,'.dcm '))/1000;
        TimeFrames(f,2)=TimeFrames(f,1)+frameLength;
        if f<numTimePoints
            TimeFrames(f+1,1)=TimeFrames(f,2);
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
        for j=1:2
            for k=2:6
                fprintf('Starting processing for %s%s; shen ID: %d\n',roiDATA{1,j},roiDATA{k,3},roiDATA{k,j})
                disp('Extracting ROI from GM_parc')
                % extract ROI from parcellation
                 roiOUT=fullfile(roiDIR,sprintf('%s%s.nii.gz',roiDATA{1,j},roiDATA{k,3}));
                sentence=sprintf('fslmaths %s -thr %d -uthr %d %s',parcFile,roiDATA{k,j},roiDATA{k,j},roiOUT);
                [~,result]=system(sentence);disp(result)
                fprintf('Extracting TAC for %s%s\n\n',roiDATA{1,j},roiDATA{k,3})
                volROI = MRIread(roiOUT); volROI=volROI.vol;
                % for each PET frame, get mean activity from mask
                for timePoint=1:numTimePoints
                    aux=reshape(volPET(:,:,:,timePoint),[sizeX,sizeY,sizeZ]); 
                    TAC(timePoint,1)=mean(aux(volROI>0));
                    clear aux
                end
                % write out time activity text file
                tac_txt=horzcat(TimeFrames,TAC);
                dlmwrite(fullfile(roiDIR,sprintf('%s%s_tac.txt',roiDATA{1,j},roiDATA{k,3})),tac_txt,'delimiter','\t','precision','%.6f')   
                clear TAC tac_txt 
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
        T1=fullfile(t1DIR,'T1_fov_denoised.nii');
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
            [~,result]=system(sprintf('run_first -i %s -t %s -o %s -n 40 -m /usr/local/fsl/data/first/models_336_bin/intref_puta/%s_Cereb.bmv -intref /usr/local/fsl/data/first/models_336_bin/05mm/%s_Puta_05mm.bmv',...
                T1,mat,cereb,side,side));disp(result)
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
        crblmFINAL=fullfile(crblmDIR,'cerebellum_noVermis.nii.gz');
        [~,result]=system(sprintf('fslmaths %s/cerebellum_bin_filled.nii.gz -mas %s_binv.nii.gz %s',crblmDIR,vermisNATIVE,crblmFINAL));disp(result)
        % Isolate cerebellar Gray Matter
        [~,result]=system(sprintf('fslmaths %s -mas %s/T1_GM_mask.nii.gz %s',crblmFINAL,t1DIR,crblmFINAL));disp(result)
        
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
        for j=1:2
            for k=2:6
                counter=counter+1;
                roiOUT=fullfile(roiDIR,sprintf('%s%s.nii.gz',roiDATA{1,j},roiDATA{k,3}));
                subjMRTMout{end+1,1}=sprintf('%s%s',roiDATA{1,j},roiDATA{k,3});
                fprintf('Running MRTM for %s\n',roiOUT)
                Ct=fullfile(roiDIR,sprintf('%s%s_tac.txt',roiDATA{1,j},roiDATA{k,3}));
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
        save(fileOUT,'subjMRTMout')
    else
        disp(parcFile)
        disp(niiDIR)
        warning('GM parcellation image above not found or/n PET nii dir not found.')
    end
end
% save results from the whole subjectlist
fileOUTall=fullfile(outDIR,sprintf('%s.mat',outFILE));
save(fileOUTall,'mrtmOUTbp','mrtmOUTk2','mrtmOUTr1','mrtmOUTk2a','mrtmOUTk2r')
close all
clearvars