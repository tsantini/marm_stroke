function niiExtractLTE(nii_in,bmatFile, bvalFile, bvecFile, saveMode, bthresh)
% Function to extract LTE scans from a b-tensor encoding scan nifti file
%    Alternate usage: provide empty array for nii_in, and files are saved
%    listing the indices related to the chosen savemode. Indices use base 0
%    (i.e., the first volume is index 0)

if nargin < 3
    bvalFile = [];
end
if nargin < 4
    bvecFile = [];
end
if (nargin < 5) || isempty(saveMode)
    % 0: save LTE + b0. bthresh is ignored for this case, and b0 is
    %    included based on code in bmatRank.m
    % 1: save LTE + b0 in one file, and STE in another file. bthresh is
    %    ignored for this case, and b0 is included based on code in bmatRank.m
    % 2: save b0, LTE, and STE in three separate files. 
    % 3: save LTE + b0 in one file and STE + b0 in another.  
    %    CAUTION: duplicates b0 images. Be sure not to use duplicated b0
    %    images in later analysis.
    saveMode = 0;
end
if (nargin < 6) || isempty(bthresh)
    % Threshold for considering scan to be b ~ 0
    bthresh = 20;
end

% Load data
info = [];
im = [];
bval = [];
bvec = [];
if ~isempty(nii_in)
    info = niftiinfo(nii_in);
    im = niftiread(nii_in);
end
[brank, bvecFromMat] = bmatRank(load(bmatFile));
if ~isempty(bvalFile)
    bval = load(bvalFile);
end
if ~isempty(bvecFile)
    bvec = load(bvecFile);
    % Replace with bvec from bmat
    sgn = sign(bvec);
    sgn(sgn==0) = 1;
    bvec = sgn.*abs(bvecFromMat);
end

% Get filename
fname = nii_in;
if isempty(fname)
    fname = bmatFile;
end
[fpath,fname] = fileparts(fname);
while contains(fname,'.')
    [~,fname] = fileparts(fname);
end

% Extract and save subsets
if saveMode < 2
    inds = brank < 1.5;
    save_nii(fpath,inds,im,info,[fname,'_LTE'],bval,bvec);
    if (saveMode == 1) && any(brank > 2.5)
        inds = brank > 2.5;
        save_nii(fpath,inds,im,info,[fname,'_STE'],bval,bvec);
        if any(and(brank>=1.5,brank<=2.5))
            inds = and(brank>=1.5,brank<=2.5);
            save_nii(fpath,inds,im,info,[fname,'_PTE'],bval,bvec);
        end
    end
elseif saveMode == 2
    inds = bval < bthresh;
    save_nii(fpath,inds,im,info,[fname,'_b0'],bval,bvec);
    %
    inds = and(bval>=bthresh,brank<1.5);
    save_nii(fpath,inds,im,info,[fname,'_LTE'],bval,bvec);
    %
    if any(brank > 2.5)
        inds = and(bval>=bthresh,brank>2.5);
        save_nii(fpath,inds,im,info,[fname,'_STE'],bval,bvec);
    end
    %
    if any(and(brank>=1.5,brank<=2.5))
        inds = and(bval>=bthresh,and(brank>=1.5,brank<=2.5));
        save_nii(fpath,inds,im,info,[fname,'_PTE'],bval,bvec);
    end
elseif saveMode == 3
    inds_b0 = bval < bthresh;
    %
    inds = or(inds_b0,brank<1.5);
    save_nii(fpath,inds,im,info,[fname,'_b0_LTE'],bval,bvec);
    %
    if any(brank > 2.5)
        inds = or(inds_b0,brank>2.5);
        save_nii(fpath,inds,im,info,[fname,'_b0_STE'],bval,bvec);
    end
    %
    if any(and(brank>=1.5,brank<=2.5))
        inds = or(inds_b0,and(brank>=1.5,brank<=2.5));
        save_nii(fpath,inds,im,info,[fname,'_b0_PTE'],bval,bvec);
    end
else
    error('unknown saveMode')
end

end

function save_nii(fpath,inds,im,info,fname,bval,bvec)
    if ~isempty(fpath)
        fpath = [fpath,filesep];
    end
    if ~isempty(im)
        im_a = im(:,:,:,inds);
        info.ImageSize(4) = size(im_a,4);
        info.raw.dim(5) = size(im_a,4);
        niftiwrite(im_a,[fpath,fname, '.nii'],info,'Compressed',true);
    else
        [~, inds_a] = find(inds);
        inds_a = inds_a - 1; % Use base 0
        fid = fopen([fpath,fname, '_inds.txt'], 'w');
        if length(inds_a)>1
            fprintf(fid, '%d,', inds_a(1:end-1));
        end
        fprintf(fid, '%d', inds_a(end));
        fclose(fid);
    end

    if ~isempty(bval)
        bval = bval(inds);
        fid = fopen([fpath,fname, '.bval'], 'w');
        fprintf(fid, '%d ', round(bval(1:end-1)));
        fprintf(fid, '%d', round(bval(end)));
        fclose(fid);
    end

    if ~isempty(bvec)
        bvec = bvec(:,inds);
        fid = fopen([fpath,fname, '.bvec'], 'w');
        for n=1:3
            fprintf(fid, '%.6f ', bvec(n,1:end-1));
            fprintf(fid, '%.6f\n', bvec(n,end));
        end
        fclose(fid);
    end
end