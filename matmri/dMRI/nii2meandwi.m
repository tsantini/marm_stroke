function dwi = nii2meandwi(name,name_bval,savename,bthresh)
% Loads in a dwi scan, and outputs a mean dwi at each b-value. Include file
% extensions in names.
%    name_bval: name of bval file or bmat file. If bmat file is provided,
%       mean dwi will be computed separately for different b-tensor shapes.
%    bthresh: variation allowed for each b-shell (s/mm2)
%

if nargin<3
    savename = [];
end
if nargin<4
    bthresh = 50;
end

% Create savename for output
if isempty(savename)
    [filepath,savename,ext] = fileparts(name);
    if strcmp(ext,'.gz')
        doCompress = 1;
        [~,savename,~] = fileparts(savename); 
    else
        doCompress = 0;
    end
    savename(savename=='.') = '_'; % niftiwrite does not like periods in name
    savename = [filepath,filesep,savename];
else
    doCompress = 1;
end

% Read in imaging data
im = niftiread(name);

% Read in b-value data
bval = load(name_bval);
brankList = [];
if min(size(bval)) > 1
    % A b-matrix has been supplied
    [brank, ~, bval] = bmatRank(bval);
    for n=1:3
        if any(brank==n)
            brankList = cat(1,brankList,n);
        end
    end
else
    brankList = 0;
    brank = zeros(1,length(bval));
end

% Convert to single
im = single(im);

% find all b-values
blist = [];
for n=1:length(bval)
    if any(abs(blist-bval(n)) < bthresh)
        continue
    else
        blist = cat(1,blist,bval(n));
    end
end
% disp(blist);

%  Create dwi's
im_info = niftiinfo(name);
im_info.PixelDimensions = im_info.PixelDimensions(1:3);
im_info.ImageSize = im_info.ImageSize(1:3);
im_info.raw.dim(1) = 3;
im_info.raw.dim(5) = 1;
im_info.raw.pixdim(5) = 0;
im_info.raw.dim_info = ' ';
%
im_info.Datatype = 'single';
im_info.BitsPerPixel = 32;
im_info.raw.datatype = 16;
im_info.raw.bitpix = 32;
for m=1:length(brankList)
    switch brankList(m)
        case 0
            txt = 'dwi';
        case 1 
            txt = 'lte';
        case 2
            txt = 'pte';
        case 3
            txt = 'ste';
    end
    for n=1:length(blist)
        inds = and(brank==brankList(m), abs(bval-blist(n))<bthresh);
        if any(inds)
            dwi = mean(im(:,:,:,inds),4);
            niftiwrite(dwi, sprintf('%s_%s_b%d', savename, txt, round(blist(n))), im_info, 'Compressed', doCompress);
        end
    end
end


