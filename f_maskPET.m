function masked_means = f_maskPET(means)
% Use T1_fov_denoised and its coresponding brain mask to brain mask the PET
% data

for m=1:2
    tmp = extractAfter(means{m},'nii_dynamic_preproc');
    PETm = extractBefore(means{m},tmp); clear tmp
    T1 = [PETm(1:end-24) 'T1'];
    
    maskDIR = fullfile(PETm,'T1mask');
    if ~exist(maskDIR,'dir')
        mkdir(maskDIR)
    end
    cd(maskDIR)
    
    % Prepare Inputs
    % T1_fov_denosied
    T1fov = '/T1_fov_denoised.nii';
    % T1_brain_mask
    Mask1 ='/T1_brain_mask_filled.nii.gz';
    [~,result]=system(sprintf('cp %s%s %s; cp %s%s %s; gunzip -f %s%s',T1,T1fov,maskDIR,T1,Mask1,maskDIR,maskDIR,Mask1));
    disp(result)
    T1mask = [maskDIR Mask1(1:end-3) ',1'];
    T1fov = [maskDIR '/T1_fov_denoised.nii,1'];
    
    fprintf('Coregistering T1 to PET%d\n',m)
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref{1,1} = means{m,1};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source{1,1} = T1fov;
    matlabbatch{1}.spm.spatial.coreg.estwrite.other{1,1} = T1mask;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    
    fprintf('Masking PET%d\n',m)
    output = sprintf('masked_meanFBP_PET%d.nii,1',m);
    masked_means{m,1}= [maskDIR '/' output];
    matlabbatch{1}.spm.util.imcalc.input{1,1} = means{m,1};
    matlabbatch{1}.spm.util.imcalc.input{2,1} = [maskDIR '/r' Mask1(2:end-3) ',1'];                                  
    matlabbatch{1}.spm.util.imcalc.output = output;
    matlabbatch{1}.spm.util.imcalc.outdir{1,1} = maskDIR;
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*(i2>0)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = -1;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end
end