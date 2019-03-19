clear all
close all
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/'))


%% generate fake fMRI data for test and save it in test_glm.nii
mask=load_nii('MNI152_T1_2mm_brain_mask.nii');
inmask=find(mask.img>0);

T=200;
img=cumsum(randn(91,109,91,T));
for t=1:T
   img(:,:,:,t)= img(:,:,:,t).*double(mask.img);
end
filename='test_glm.nii';
nii=make_nii(img);
save_nii(nii,filename);

%% run the non parametric GLM
N=2; % number of regressors in the model
cfg=[];
cfg.infile='test_glm.nii';
cfg.regressor = cumsum(randn(T,N)); % put your model here
cfg.cdtP = 0.05;
%   cft.cdtR = r value of cluster defining threshold (overrides cdtP)
cfg.alpha = 0.05; % % alpha level of the permutation test
cfg.clusterstatistic = 'maxsum'; % 'maxsize' or 'maxsum'
cfg.NPERM = 5000;
cfg.seed = 0; % seed for the random
cfg.toi = find(squeeze(cfg.regressor)); % time points of interest
cfg.CR=[0 1]; % range to estimate correlations

cfg=bramila_glm_np(cfg);

% store the results into a new nifti file

for r=1:N
    out=cfg.vol(:,:,:,r).*sign(cfg.cmask(:,:,:,r));
    save_nii(make_nii(out),['results_regressor_' num2str(r) '.nii']);
end
