function f_spm_coreg(ref,source,other)

spm_figure('GetWin','Graphics');
        
matlabbatch{1}.spm.spatial.coreg.estimate.ref{1,1} = ref;
matlabbatch{1}.spm.spatial.coreg.estimate.source{1,1} = source;
matlabbatch{1}.spm.spatial.coreg.estimate.other = other;

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

spm_jobman('run',matlabbatch);
clear matlabbatch

end
