function make_hist (pid, cond, num);
BP = load_nii (sprintf ('%s_%s_%d_010100_MRTM_BP.nii', pid, cond, num));
R1 = load_nii (sprintf ('%s_%s_%d_010100_MRTM_R1.nii', pid, cond, num));
k2 = load_nii (sprintf ('%s_%s_%d_010100_MRTM_k2.nii', pid, cond, num));
k2a = load_nii (sprintf ('%s_%s_%d_010100_MRTM_k2a.nii', pid, cond, num));
nii=BP;
BP=BP.img;
R1=R1.img;
k2=k2.img;
k2a=k2a.img;
nvox = prod (size (BP));
BPreshape = reshape (BP, nvox, 1);
R1reshape = reshape (R1, nvox, 1);
k2reshape = reshape (k2, nvox, 1);
k2areshape = reshape (k2a, nvox, 1);
keyboard
idx = find (BPreshape>1 & R1reshape>0.5 & k2reshape>0.05 & k2areshape>0.01);
nmask = length (idx)
k2r = median (k2reshape(idx) ./ R1reshape(idx))
figure; hold on;
hist (k2reshape(idx) ./ R1reshape(idx), 50, 1); colormap(white)
plot (k2r, 0, 'rx', 'MarkerSize', 12);
grid on;
xlabel ('k2ref (1/min)');
ylabel ('Frequency of occurrance');
title (sprintf ('Histogram of k2ref values for the %d voxels passing thresholds', nmask));
print ('-dpng', sprintf ('k2r-%s-%s.png', pid, cond));
mask = zeros (size (BPreshape));
mask(idx) = 1;
nii.img = reshape (mask, size (BP));
save_nii (nii, sprintf ('mask-%s-%s.nii', pid, cond));
