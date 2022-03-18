function racPET_spm12_dropFrame_TACandMRTM(dataDIR,subjectID,scan,indices2drop)
%%
% This function re-runs MRTM regional BP estimation, allowing to drop
% frames from the fit.
% 
% It required that racPET_spm12_2_TACandMRTM.m has been ran and the 
% time-activity curves for all ROI and cerebellum have been saved.
%
% Existing TAC files are read in, time-points are dropped, and the new TAC
% files and their MRTM output are save to a new directory.
%
% INPUTS:
%   dataDIR - path to data directry that contains the subject to be re-run.
%             (same as that specified in racPET_spm12_2_TACandMRTM.m)
%
%   subjectID - a single subject ID (directory name for subject) that
%               contains the scan data.
%
%   scan - 'PET' or 'PET1' for single scan designs; or 'PET1' or 'PET2' for
%           dual scan designs. The name of the subdirectory within the
%           subject directory that contain the PET processing data.
%
%   indices2drop - a numberic row vector of frames to be excluded from the
%                  fit.
%
% Example Usage:
%   racPET_spm12_dropFrame_TACandMRTM('/projects/pet_processing/datadir','W2D0001','PET1',[45 47 49])   
%
% Evgeny Chumin, IU BLoomington, 2020
%% 
% set system specific paths
addpath(genpath('/N/project/HCPaging/yoderBP_project/PET_processing_Code/toolbox_matlab_nifti'))
addpath(genpath('/N/project/HCPaging/yoderBP_project/PET_processing_Code/yapmat-0.0.3a2-ec/src'))
% Raclopride half-life
thalf=20.4;

  roiDATA= {'L_', 'R_', 'label';
            25,   9, 'NAc-shell'
            26,   10, 'NAc-core'
            29,   13, 'aPUT'
            30,   14, 'pPUT'
            31,   15, 'aCAU'
            32,   16, 'pCAU'};
        
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
%%
% set directory paths
dropString=sprintf('-%d',indices2drop);
petDIR=fullfile(dataDIR,subjectID,scan);  
roiDIR_all=fullfile(petDIR,'roi_TAC_mrtm');

if exist(roiDIR_all,'dir') && exist(fullfile(roiDIR_all,'cerebellum_tac.txt'),'file')
    roiDIR_drop=fullfile(petDIR,['roi_TAC_mrtm_drop' dropString]);
    if ~exist(roiDIR_drop,'dir')
        mkdir(roiDIR_drop)
    end

    % get list of tac files to edit
    tacList=dir([roiDIR_all '/*tac.txt']);
    indices2drop=sort(indices2drop,'descend');
    
    % read in tac, remove rows corresponding to dropped frames, write new tac
    for ii=1:length(tacList)
        tacDATA=dlmread(fullfile(roiDIR_all,tacList(ii).name));
        for jj=1:length(indices2drop)
            tacDATA=vertcat(tacDATA(1:indices2drop(jj)-1,:),tacDATA(indices2drop(jj)+1:end,:));
        end
        dlmwrite(fullfile(roiDIR_drop,tacList(ii).name),tacDATA,'delimiter','\t','precision','%.6f')
        clear tacDATA
    end
    
    % reference region
    Cr=fullfile(roiDIR_drop,'cerebellum_tac.txt');
    
    % pleallocate MRTM outputs
    subjMRTMout={'ROI','BP','R1','k2','k2a','k2r'};
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
            Ct=fullfile(roiDIR_drop,sprintf('%s_tac.txt',rl));
            [BP, R1, k2, k2a, k2r] = mrtm (Ct, Cr, roiDIR_drop, thalf, 'conventional');
            subjMRTMout{end,2}=BP;  
            subjMRTMout{end,3}=R1;  
            subjMRTMout{end,4}=k2;  
            subjMRTMout{end,5}=k2a; 
            subjMRTMout{end,6}=k2r; 
            clear BP R1 k2 k2a k2r
            close all
        end
    end
    clear counter

    % save subject results
    fileOUT=fullfile(roiDIR_drop,sprintf('%s_mrtm_output.mat',subjectID));
    count=length(dir(strcat(fileOUT(1:end-4),'*')));
    if count>0
        fileOUT=fullfile(roiDIR_drop,sprintf('%s_mrtm_output_run%d.mat',subjectID,count+1));
    end
    save(fileOUT,'subjMRTMout')
else 
    fprintf(2,'MRTM_roi time activity curve directory or\n')
    fprintf(2,'contained cerebellum text file do not exist!\n')
end

%% --------------------------------------------------------------------- %%
% get figure names
modelfit=dir(fullfile(roiDIR_drop,'*modelfit.fig'));
wresid=dir(fullfile(roiDIR_drop,'*wresid.fig'));

% load in the figures and get handle of axes and all children
for i=1:length(modelfit)
    h{i,1}=openfig(fullfile(modelfit(i).folder,modelfit(i).name),'reuse');
    ax{i,1}=gca; pause(.5)
    fig{i,1}=get(ax{i,1},'children');
    fig{i,1}(2).MarkerSize=4; fig{i,1}(3).MarkerSize=4;
    h{i,2}=openfig(fullfile(wresid(i).folder,wresid(i).name),'reuse');
    ax{i,2}=gca;pause(.5)
    fig{i,2}=get(ax{i,2},'children');
    fig{i,2}(2).MarkerSize=4;
end

o1=figure('Position',[1 1 675 900]);

m=0;r=0;
for i=1:length(h)
    if mod(i,2)==1
        m=m+1;
        s1{m,1}=subplot(length(h)/2,2,i);
        title(strrep(extractBefore(modelfit(m).name,'_tac_MRTM_modelfit'),'_','-'))
    elseif mod(i,2)==0
        r=r+1;
        s1{r,2}=subplot(length(h)/2,2,i);
    end
end
clear m r
sgtitle(o1,[subjectID,' ',scan,' Page 1/2']);
%p1=mtit(o1,strcat(subjectID,' ',scan,' Page 1/2'));

for i=1:length(h)/2
    copyobj(fig{i,1},s1{i,1})
    if i==3
        ylabel(s1{i,1},'Tracer Concentration')
    end
    if i==5
        xlabel(s1{i,1},'Time (min)')
    end
    copyobj(fig{i,2},s1{i,2})
    if i==3
        ylabel(s1{i,2},'Weighted Residuals')
    end
    if i==5
        xlabel(s1{i,2},'Time (min)')
    end
end

print(fullfile(roiDIR_drop,sprintf('%s-MRTMfits_%s_p1of2.pdf',subjectID,scan)),'-dpdf','-fillpage')

o2=figure('Position',[1 1 675 900]);

m=0; r=0;
for i=1:length(h)
    if mod(i,2)==1
        m=m+1;
        t=m+length(h)/2;
        s2{m,1}=subplot(length(h)/2,2,i);
        title(strrep(extractBefore(modelfit(t).name,'_tac_MRTM_modelfit'),'_','-'))
    elseif mod(i,2)==0
        r=r+1;
        s2{r,2}=subplot(length(h)/2,2,i);
    end
end
clear m r
sgtitle(o2,[subjectID,' ',scan,' Page 2/2']);

for i=1:length(h)/2
    k=i+length(h)/2;
    copyobj(fig{k,1},s2{i,1})
    if i==3
        ylabel(s2{i,1},'Tracer Concentration')
    end
    if i==5
        xlabel(s2{i,1},'Time (min)')
    end
    copyobj(fig{k,2},s2{i,2})
    if i==3
        ylabel(s2{i,2},'Weighted Residuals')
    end
    if i==5
        xlabel(s2{i,2},'Time (min)')
    end
end

print(fullfile(roiDIR_drop,sprintf('%s-MRTMfits_%s_p2of2.pdf',subjectID,scan)),'-dpdf','-fillpage')

close all
end