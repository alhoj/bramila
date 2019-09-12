function nii=bramila_fixOriginator(targetimg,refimg)
% Usage:
%   targetimg='/m/nbe/scratch/braindata/something.nii'
%   refimg='/m/nbe/scratch/braindata/shared/toolboxes/HarvardOxford/MNI152_T1_2mm_brain.nii';
%   nii=bramila_fixOriginator(targetimg,refimg);
%   save_nii(nii,filename);
%

if ischar(refimg)
    ref=load_nii(refimg);
else
    ref=refimg;
end

target=load_nii(targetimg);


target.hdr.dime.pixdim = ref.hdr.dime.pixdim;
target.hdr.dime.scl_slope = ref.hdr.dime.scl_slope;
target.hdr.dime.xyzt_units = ref.hdr.dime.xyzt_units;
target.hdr.hist = ref.hdr.hist;


target.original.hdr.dime.pixdim = ref.original.hdr.dime.pixdim;
target.original.hdr.dime.scl_slope = ref.original.hdr.dime.scl_slope;
target.original.hdr.dime.xyzt_units = ref.original.hdr.dime.xyzt_units;
target.original.hdr.hist = ref.original.hdr.hist;

%fileout = [target.fileprefix '_W_ORIGINATOR.nii'];

%save_nii(target,fileout);

nii=target;
