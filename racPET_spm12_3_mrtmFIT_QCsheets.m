 % set data directory paths
dataDIR='/datay2/chumin-F31/data/CNT/SKYRA';
outDIR='/datay2/chumin-F31/results/mrtm_qc_sheets/CNT-SKYRA';
scan='PET';

if ~exist(outDIR,'dir')
    mkdir(outDIR)
end

% Loop accross subjects        
subjDIRS=dir(dataDIR);subjDIRS(1:2)=[];

for s=1:length(subjDIRS)
    datadir=fullfile(subjDIRS(s).folder,subjDIRS(s).name,scan,'roi_TAC_mrtm');
% get figure names
modelfit=dir(fullfile(datadir,'*modelfit.fig'));
wresid=dir(fullfile(datadir,'*wresid.fig'));

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
p1=mtit(o1,strcat(subjDIRS(s).name,' 1/2'));

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
    
print(fullfile(datadir,sprintf('%s-MRTMfitsV2_1of2.pdf',subjDIRS(s).name)),'-dpdf','-fillpage')
copyfile(fullfile(datadir,sprintf('%s-MRTMfitsV2_1of2.pdf',subjDIRS(s).name)),outDIR)

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
p2=mtit(o2,strcat(subjDIRS(s).name,' 2/2'));

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

print(fullfile(datadir,sprintf('%s-MRTMfitsV2_2of2.pdf',subjDIRS(s).name)),'-dpdf','-fillpage')
copyfile(fullfile(datadir,sprintf('%s-MRTMfitsV2_2of2.pdf',subjDIRS(s).name)),outDIR)

close all
end