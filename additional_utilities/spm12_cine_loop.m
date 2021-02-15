function spm12_cine_loop(path2data,fileroot)

% If spm12 is not in matlab path prior to execution, add it here.
addpath(genpath('/usr/local/spm12/'))
spm_figure('GetWin','Graphics');
fileList=dir(path2data);
fileList(1:2)=[];
 
for i=1:length(fileList)
    indx = strfind(fileList(i).name,fileroot);
    if ~isempty(indx) && ~exist('fileFrames','var')==1
        if ~isempty(strfind(fileList(i).name,'nii'))
            if indx==1
                fileFrames(1,:)=[path2data '/' fileList(i).name];
            end
        else
            warning('Fileroot match found, but it is not a nifti image.')
            disp(fileList(i).name)
        end
    elseif ~isempty(indx) && exist('fileFrames','var')==1
        if ~isempty(strfind(fileList(i).name,'nii'))
            if indx==1
                fileFrames(end+1,:)=[path2data '/' fileList(i).name];
            end
        else
            warning('Fileroot match found, but it is not a nifti image.')
            disp(fileList(i).name)
        end
    end
    
end
if ~exist('fileFrames','var')==1 || isempty(fileFrames)
    warning('No files matching provided fileroot found');return
end
fprintf('%d images met fileroot criteria.\n',size(fileFrames,1))
        
spm_check_registration(fileFrames(:,:));