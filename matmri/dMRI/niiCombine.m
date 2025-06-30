function niiCombine(nii1,nii2,inds1,inds2,outname,bvec1,bvec2,bval1,bval2)
% Combines nii files nii1 and nii2, where the indices the volumes
% occupy in the output are given by the files inds1 and inds2, respectively.
%
% Intended use: first use niiExtractLTE to separate LTE and STE (savemode
% 3), then apply fsl eddy on LTE and STE separately, then recombine with
% niiCombine.
% - the indices in the files inds1 and inds2 should use 0 for first index
% - if indices are in both inds1 and inds2, only the ones from nii1 and
% inds1 will be used.
% - filenames must include extensions! (e.g., .nii or .nii.gz)
% optional: supply bvec and bval files to be combined.
%
% (C) Corey Baron, 2024

% Load data
im_info = niftiinfo(nii1);
im1 = niftiread(nii1);
im2 = niftiread(nii2);
inds1_a = load(inds1);
inds2_a = load(inds2);

% Build output
sz = size(im1);
sz(4) = max([inds1_a(:); inds2_a(:)]) + 1;
imout = zeros(sz,'like',im1);
imout(:,:,:,inds2_a+1) = im2;
imout(:,:,:,inds1_a+1) = im1;

% Save the combined nii file
if size(im1,4) == 1 
    % Input was 3D, so there is more to change to make it 4D
    if length(im_info.PixelDimensions) < 4
        warning('using arbitrary temporal spacing in nifti')
        im_info.PixelDimensions = [im_info.PixelDimensions, 5];
    end
    im_info.ImageSize = sz;
    im_info.raw.dim(1) = 4;
    im_info.raw.dim(5) = sz(4);
    im_info.raw.pixdim(5) = im_info.PixelDimensions(4);
else
    im_info.ImageSize(4) = sz(4);
    im_info.raw.dim(5) = sz(4);
end
if isequal(outname(end-2:end),'.gz')
    outname = outname(1:end-3);
    niftiwrite(imout, outname, im_info, 'Compressed', true);
else
    niftiwrite(imout, outname, im_info);
end

% Save combined bvec files
if (nargin > 5) && ~isempty(bvec1)
    bvec1_a = load(bvec1);
    bvec2_a = load(bvec2);
    bvecOut = zeros(3, sz(4));
    bvecOut(:,inds2_a+1) = bvec2_a;
    bvecOut(:,inds1_a+1) = bvec1_a;
    fid = fopen([outname(1:end-4), '.bvec'], 'w');
    for n=1:3
        fprintf(fid, '%.6f ', bvecOut(n,1:end-1));
        fprintf(fid, '%.6f\n', bvecOut(n,end));
    end
    fclose(fid);
end

% Save combined bval files
if (nargin > 7) && ~isempty(bval1)
    bval1_a = load(bval1);
    bval2_a = load(bval2);
    bvalOut = zeros(3, sz(4));
    bvalOut(:,inds2_a+1) = bval2_a;
    bvalOut(:,inds1_a+1) = bval1_a;
    fid = fopen([outname(1:end-4), '.bval'], 'w');
    for n=1:3
        fprintf(fid, '%.6f ', bvalOut(n,1:end-1));
        fprintf(fid, '%.6f\n', bvalOut(n,end));
    end
    fclose(fid);
end


