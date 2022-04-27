function f_edit_Z_slices(frames,numslices)
% remove or pad slices on the inferior edge of an image.
% REQUIRED: AFNI installed and in your bashrc
% INPUTS:
%   frames - structure of pet image files from dir command 
%   numslices - number of slices to be cut (negative values) or added
%               positive values.
% OUTPUT:
%   Overwrites existing raw data with new field of view

for f = 1:length(frames)
    prefix = fullfile(frames(f).folder,['c' frames(f).name]);
    fileIN = fullfile(frames(f).folder,frames(f).name);
    [~,result1]=system(sprintf('3dZeropad -I %d -prefix %s %s',numslices,prefix,fileIN));
    [~,result2]=system(sprintf('mv %s %s',prefix,fileIN));
    if f == 1
        disp(result1)
        disp(result2)
    end
end