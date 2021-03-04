% racPET_spm12_3_mrtmFIT_QCsheets.m
%
% This script follows the 2_TACand MRTM, consolidating all figures into pdf
% file summary sheets of MRTM fits for all striatal regions.
%
% Contributors:
% Evgeny Chumin, Indiana University School of Medicine, 2019
%                Indiana University, Bloomington, 2020
% Mario Dzemidzic, Indiana University School of Medicine, 2019
%-------------------------------------------------------------------------%
%% set data directory paths
dataDIR='/projects/pet_processing/datadir';
outDIR='/projects/pet_processing/datadir_out';
%-------------------------------------------------------------------------%
%% Subject list selection.
% Run all subjects:
    %subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];
% Run a single or set of subjects:
    subjDIRS=dir([dataDIR '/*01']);
%-------------------------------------------------------------------------%   
%% End of user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(outDIR,'dir')
    mkdir(outDIR)
end
%% Looping accross subjects        
for s=1:length(subjDIRS)
    disp(subjDIRS(s).name)
    dircont=dir(fullfile(subjDIRS(s).folder,subjDIRS(s).name)); dircont(1:2)=[];
    petList=struct.empty;
    for p=1:length(dircont)
        if dircont(p).isdir==1 && ~isempty(strfind(dircont(p).name,'PET'))
            petList(end+1).name=dircont(p).name;
        end
    end
%% Loop across PET sessions
for p=1:length(petList) % loop over PET scans
    disp(petList(p).name)
    datadir=fullfile(subjDIRS(s).folder,subjDIRS(s).name,petList(p).name,'roi_TAC_mrtm');
% get figure names
modelfit=dir(fullfile(datadir,'*modelfit.fig'));
wresid=dir(fullfile(datadir,'*wresid.fig'));
nr = size(modelfit,1); %number of roi

% load in the figures and get handle of axes and all children
for i=1:nr
    h{i,1}=openfig(fullfile(modelfit(i).folder,modelfit(i).name),'reuse');
    ax{i,1}=gca; pause(.5)
    fig{i,1}=get(ax{i,1},'children'); pause(.5)
    fig{i,1}(2).MarkerSize=4; fig{i,1}(3).MarkerSize=4;
    h{i,2}=openfig(fullfile(wresid(i).folder,wresid(i).name),'reuse');
    ax{i,2}=gca; pause(.5)
    fig{i,2}=get(ax{i,2},'children'); pause(.5)
    fig{i,2}(2).MarkerSize=4;
end

nump = ceil(nr/5);
for pg=1:nump
    switch pg
        case nump % last page
            lft = nr-(5*(pg-1)); % number of remaining regions
            ofig{pg}=figure('Position',[1 1 675 900]);
            m=0;r=0;
            for i=1:lft*2
                if mod(i,2)==1
                    m=m+1;
                    t=m+(5*(pg-1));
                    sfig{m,1}=subplot(5,2,i);
                    title(strrep(extractBefore(modelfit(t).name,'_tac_MRTM_modelfit'),'_','-'))
                elseif mod(i,2)==0
                    r=r+1;
                    sfig{r,2}=subplot(5,2,i);
                end
            end
            clear m r
            sgtitle(ofig{pg},sprintf('%s %s Page %d/%d',subjDIRS(s).name,petList(p).name,pg,nump));

            for i=1:lft
                k=i+5*(pg-1);
                copyobj(fig{k,1},sfig{i,1})
                copyobj(fig{k,2},sfig{i,2})
                if i==ceil(lft/2)
                    ylabel(sfig{i,1},'Tracer Concentration')
                    ylabel(sfig{i,2},'Weighted Residuals')
                end
                if i==lft
                    xlabel(sfig{i,1},'Time (min)')
                    xlabel(sfig{i,2},'Time (min)')
                end
            end
            fileout=fullfile(datadir,sprintf('%s-MRTMfits_%s_p%dof%d.pdf',subjDIRS(s).name,petList(p).name,pg,nump));
            print(ofig{pg},'-dpdf',fileout,'-fillpage')
            copyfile(fileout,outDIR)
            clear sfig fileout lft

        otherwise % all preceeding pages
            ofig{pg}=figure('Position',[1 1 675 900]);
            m=0;r=0;
            for i=1:10
                if mod(i,2)==1
                    m=m+1;
                    t=m+5*(pg-1);
                    sfig{m,1}=subplot(5,2,i);
                    title(strrep(extractBefore(modelfit(t).name,'_tac_MRTM_modelfit'),'_','-'))
                elseif mod(i,2)==0
                    r=r+1;
                    sfig{r,2}=subplot(5,2,i);
                end
            end
            clear m r
            sgtitle(ofig{pg},sprintf('%s %s Page %d/%d',subjDIRS(s).name,petList(p).name,pg,nump));

            for i=1:5
                k=i+5*(pg-1);
                copyobj(fig{k,1},sfig{i,1})
                copyobj(fig{k,2},sfig{i,2})
                if i==3
                    ylabel(sfig{i,1},'Tracer Concentration')
                    ylabel(sfig{i,2},'Weighted Residuals')
                end
                if i==5
                    xlabel(sfig{i,1},'Time (min)')
                    xlabel(sfig{i,2},'Time (min)')
                end
            end
            fileout=fullfile(datadir,sprintf('%s-MRTMfits_%s_p%dof%d.pdf',subjDIRS(s).name,petList(p).name,pg,nump));
            print(ofig{pg},'-dpdf',fileout,'-fillpage')
            copyfile(fileout,outDIR)
            clear sfig fileout
    end
end
close all
end
end