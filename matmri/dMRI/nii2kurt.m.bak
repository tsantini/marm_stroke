function [Dmean,Dpar,Dperp,Wmean,Kpar,Kperp,FA,fshells,FAvec,Dpowder,Wpowder] = nii2kurt(file, bfiles, freqs, maskfile, opt, savename)
%
%   Estimate kurtosis using axially symmetric model. See Hansen B, Shemesh N, Jespersen SN. Fast imaging of mean, axial and radial diffusion kurtosis. Neuroimage. 2016 Nov 15;142:381-393. doi: 10.1016/j.neuroimage.2016.08.022. Epub 2016 Aug 15. PMID: 27539807; PMCID: PMC5159238.
%   Supports multifreq OGSE.
%   Includes spatial regularization (optional)
%   Requires matMRI toolbox: https://gitlab.com/cfmm/matlab/matmri
%
%   Two step algorithm:
%   1. a. Solves argmin ||Ax-b||^2_2 + opt.TVregVec||Dx||^2_2
%       A = basic DTI model
%       x = diffusion tensor (one tensor is fit over all OGSE freqs)
%       b = log-transformed signal data
%       D = 3D spatial gradient
%      b. determines eigenvector for domininant eigenvalue, which is used
%      to determing the voxel-wise axis of symmetry
%
%   2. Solves argmin ||Ax-b||^2_2 + opt.TVreg||Dx||^2_2 + opt.TVreg*opt.tvFreq||DGx||^2_2
%       A = diffusion kurtosis model with axial symmetry
%       x = diffusion kurtosis parameter maps
%       b = log-transformed signal data
%       D = 3D spatial gradient
%       G = 1D spatial gradient along OGSE freq dimension, applied circularly
%
%   INPUTS
%   file:       path to NIFTI file containing raw image data
%   bfiles:     path to .bval and .bvec files (in same format as FSL expects 
%                   and dcm2niix outputs). Include name, but not extension
%               Units: s/mm2 (important for spatial regularization so that
%               D and K are scaled properly relative to each other).
%   freqs:      path to file in the same format as .bval file, but instead specifying 
%                   OGSE frequencies (use 0 for PGSE). Only used for
%                   binning to get estimates at each freq. Units of Hz.
%                   Alternate Usage: the user may also input a matlab vector
%   maskfile:   (optional) path to NIFTI file containing a binary mask
%   opt:        (optional) structure that contains options for algorithm.
%                   See comments where defaults are set in code for info.
%   savename:   (optional) filename for output nifti. Exclude extension.
%                It can be the fullpath+filename for a different output directory.
%
%   OUTPUTS
%   Dmean,Dpar,Dperp,Wmean,Kpar,Kperp,FA,fshells,FAvec,Dpowder,Wpowder
%      See Hansen et al for definitions of Dmean,Dpar,Dperp,Wmean,Kpar,Kperp
%   FA is the usual definition of fractional anisotropy
%   fshells: the ogse frequency at each shell (not an image!)
%   FAvec:    mean FA over all OGSE freqs scaled by dominant eigenvector direction
%   Dpowder, Wpowder: apparent mean diffusivity and kurtosis obtained from
%                     powder average fitting
%
%   (c) 2023, Corey Baron
%

% Start timer
tic_a = tic;

%% Set option defaults
if nargin < 4
    maskfile = [];
end

% Options 
if nargin<6
    savename = [];
end
if nargin<5 || ~isfield(opt,'bthresh') || isempty(opt.bthresh)
    % Threshold for difference between b-shells
    opt.bthresh = 50;
end
if nargin<5 || ~isfield(opt,'fthresh') || isempty(opt.fthresh)
    % Threshold for difference between OGSE frequency  shells
    opt.fthresh = 5;
end
if ~isfield(opt,'noGPU') || isempty(opt.noGPU)
    % Disable automatic usage of GPU
    opt.noGPU = 0;
end
if ~isfield(opt,'saveNifti') || isempty(opt.saveNifti)
    % Save nifti outputs. 
    opt.saveNifti = true;
end
if ~isfield(opt,'verbose') || isempty(opt.verbose)
    % Print out information
    opt.verbose = 1;
end
if ~isfield(opt,'nIter') || isempty(opt.nIter)
    % Number of iterations 
    opt.nIter = 100;
end
if ~isfield(opt,'tvRegVec') || isempty(opt.tvRegVec)
    % Whether to do spatial regularization for diffusion tensor calc for
    % finding main eigenvector
    opt.tvRegVec = 0.5;
end
if ~isfield(opt,'tvReg') || isempty(opt.tvReg)
    % Whether to do spatial regularization for main calculations
    opt.tvReg = 0.1;
end
if ~isfield(opt,'tvFreq') || isempty(opt.tvFreq)
    % What to scale TV^2 of differences betweem ogse frequencies for during TV_^2 regularization 
    % Set to 0 to not penalize contrast differences across freqs
    % IN DEVELOPMENT - NOT RECOMMENDED
    opt.tvFreq = 0;
end

fprintf('nii2kurt processing...\n')
tic1 = tic;

% Check for gpu support
useGPU = 0;
if gpuDeviceCount >= 1 && ~opt.noGPU
    if opt.verbose
        fprintf('  Using GPU...\n')
    end
    G = gpuDevice;
    useGPU = 1;
end

%% Read in files
im = niftiread(file);
im = single(im); % convert to single. Log function doesnt work with integer types.
im(im<10E-15)=10E-15; %We can get negative values after preprocessing
szIm = size(im);
bval = load([bfiles,'.bval']); 
bvec = load([bfiles,'.bvec']);
bval = bval(:);
bvec = bvec.';
bvec = bvec./sum(bvec.^2,2); % normalize
bvec(isnan(bvec))=0; % avoid nan if bvec=[0 0 0] (b0)

if ischar(freqs)
    % Text file, in the same format as .bval file, which indicates which
    % acquisitions are STE (1) and LTE (0). If this is not a string, we
    % presume the user inputted the list directly
    freqs = load(freqs);
elseif isempty(freqs)
    freqs = zeros(size(bval),'like',bval);
end
freqs = freqs(:);
if isempty(maskfile)
    mask = true(szIm(1:3));
else
    mask = niftiread(maskfile);
    mask = logical(mask);
end

% Divide b-values by 1000 to make them close to 1 (better properties for
% numerical calculations)
bval = bval/1000;
opt.bthresh = opt.bthresh/1000;

% Perform some checks
if max(bval) > 3
    warning('b-value greater than 3000 s/mm2 detected')
end
if opt.tvReg > 0 || opt.tvRegVec > 0
    fprintf('  INFO: using spatial regularization; tvReg = %g; tvRegVec = %g\n', opt.tvReg, opt.tvRegVec)
end
if opt.tvFreq > 0
    warning('Using OGSE freq regularization - not well tested and not recommended.')
end

% Crop mask to exclude voxels that are 0 in all acquisitions (can happen
% from eddy current correction)
mask = and(mask, sum(abs(im),4) > 2*eps);

% Pixel values of 0 can lead to NaN that propagate throughout the image due
% to the spatial regularization.
im(im==0) = eps;

% Move to GPU
if ~opt.noGPU
    bval = gpuArray(bval);
    bvec = gpuArray(bvec);
end

% Find and sort b-value and frequency shells
bshells = shellSort(bval,opt.bthresh);
fshells = shellSort(freqs,5);

% Estimation of diffusion tensor for goal of getting main eigenvector.
% Consider all freqs together, since eigenvector is not expected to depend
% on freq. Ditto for b-values
inds = true(size(bval)); % to select subset, can use something like inds = bval < (1 + opt.bthresh); 
A = [-ones(sum(inds),1), bval(inds).*(bvec(inds,1).^2), bval(inds).*(bvec(inds,2).^2), ...
    bval(inds).*(bvec(inds,3).^2), 2*bval(inds).*bvec(inds,1).*bvec(inds,2),...
    2*bval(inds).*bvec(inds,1).*bvec(inds,3), 2*bval(inds).*bvec(inds,2).*bvec(inds,3)];
S = -log(im(:,:,:,inds));
S = permute(S, [4,1,2,3]);
DT = zeros([size(A,2),szIm(1:3)],'like', A);
DT(:,mask) = A\S(:,mask);
if opt.tvRegVec > 0
    gamRelVec = [0 1 1 1 2 2 2]'; % [s0 Dxx Dyy Dzz Dxy Dxz Dyz]
    NitMax = 500;
    szin = size(im);
    szin = szin(1:3);
    bin = 0*S(:,:);
    bin(:,mask) = S(:,mask);
    x0 = DT;
    x0 = reshape(x0,[size(x0,1), szin]);
    optC.resThresh = 1e-5;
    % Determine ratio between data consistency and regularizer
    Rin = @(x,transp) spDiff(x,transp,[0 1 1 1]);
    x0diff = sqrt(gamRelVec).*Rin(x0,'notransp');
    resid = 0*bin;
    resid(:,mask) = A*x0(:,mask) - bin(:,mask);
    norm2diff  = sum(x0diff(:).^2);
    norm2resid = sum(resid(:).^2);
    gam = opt.tvRegVec*gamRelVec*norm2resid/norm2diff;
    [DT, resSqAll, ~,~,~,~,status] = cgne(@(x,transp) AsubDTI(x,transp,A,szin,mask),...
        bin,x0,NitMax,optC,gam,Rin);
end
clear S x0
DT = cat(1, DT(2,:), DT(5,:), DT(6,:), DT(5,:), DT(3,:), DT(7,:), DT(6,:), DT(7,:), DT(4,:));
DT = reshape(DT,[3 3 size(DT,2)]);

% Compute eigenvalues and eigenvectors. Using svd because pagefun does not
% support eig. Also, svd conveniently sorts the singular values
if ~isMATLABReleaseOlderThan('R2022b') && ~opt.noGPU
    vec = zeros(size(DT),'like',DT);
    eigvals = zeros(size(DT),'like',DT);
    [vec(:,:,mask), eigvals(:,:,mask)] = pagefun(@svd,DT(:,:,mask)); 
    vec = reshape(vec(:,1,:), [3, size(DT,3)]);
    eigvals = cat(1,eigvals(1,1,:),eigvals(2,2,:),eigvals(3,3,:));
else
    warning('nii2kurt is MUCH faster with MATLAB R2022b or higher with gpu')
    vec = zeros(3,size(DT,3),'like',DT);
    eigvals = zeros(3,size(DT,3),'like',DT);
    for n=1:size(DT,3)
        if mask(n)
            [tmp, tmp2] = svd(DT(:,:,n));
            vec(:,n) = tmp(:,1);
            eigvals(:,n) = diag(tmp2);
        end
    end
end
clear DT
vec = vec./sum(vec.^2,1); % Normalize (shouldn't be necessary)
apparentMD = mean(eigvals,1);
apparentFA = sqrt(3/2*sum((eigvals - apparentMD).^2,1)./...
    sum(eigvals.^2,1)); % NB: this is not accurate if all b-values were used.

% Compute full model. The cosine of the polar angle theta used in Eq. 15 of
% Hansen et al (full citation above) is given by the dot product of the
% dominant eigenvector and the diffusion encoding direction (which are both
% unit vectors)
cosThet = sum(reshape(vec, [1 3 size(vec,2)]).*bvec, 2);
vec = gather(vec);
cos2Thet = 2*cosThet.^2 - 1;
cos4Thet = 2*cos2Thet.^2 - 1;

% Create system matrix. Columns correspond to:
% S0, Derp, Dpar, Wperp, Wpar, Wmean
A = cat( 2, ones(size(cosThet),'like',cosThet), ...
    -bval.*(1-cosThet.^2), -bval.*cosThet.^2, ...
    bval.^2/6.*(10/16*cos4Thet - 8/16*cos2Thet -2/16), ...
    bval.^2/6.*(5/16*cos4Thet + 8/16*cos2Thet+3/16), ...
    bval.^2/6.*(-15/16*cos4Thet + 15/16) );
clear cosThet cos2Thet cos4Thet

% Solve system for all freqs
Dpowder = zeros([szIm(1:3), length(fshells)],'like',im);
Wpowder = zeros([szIm(1:3), length(fshells)],'like',im);
S = log(im);
S = permute(S, [4 5 1 2 3]);
% Average all the b0's so that all freqs use the same b0. We also have to
% give some b0 to each freq, so that looping through freqs works in
% AsubAllf.
freqs0 = freqs;
bval0 = bval;
S(bval<opt.bthresh,:,:,:,:) = repmat(mean(S(bval<opt.bthresh,:,:,:,:),1), ...
    [sum(bval<opt.bthresh) 1 1 1 1]);
indsb0 = bval < opt.bthresh;
Sb0 = S(bval<opt.bthresh,:,:,:,:);
Ab0 = A(bval<opt.bthresh,:,:);
freqs(indsb0) = fshells(1);
b0 = mean(bval(bval<opt.bthresh));
bval(bval<opt.bthresh) = b0;
for nf = 2:length(fshells)
    S = cat(1,S,Sb0);
    freqs = cat(1,freqs(:),fshells(nf)*ones(sum(indsb0),1));
    bval = cat(1,bval(:),b0*ones(sum(indsb0),1));
    A = cat(1,A,Ab0);
end
% Get initial guess
y = AsubAllf(S,'pinv',A,szIm(1:3),freqs,fshells,opt.fthresh,mask); 
if opt.tvReg > 0
    gamRel = [0 1 1 1 1 1]'; % [s0 Dperp Dpar Wperp Wpar Wmean]
    NitMax = 500;
    opt.resThresh = 1e-5;
    bin = 0*S(:,:,:);
    bin(:,:,mask) = S(:,:,mask);
    % Find residual ratio to help scale gamma. We do not include
    % differences across ogse freqs for the scaling so that adding in ogse
    % freq diffs predictably increases net regularization
    if opt.tvFreq > 0
        Rin = @(x,transp) RsubFreq(x, transp, opt.tvFreq);
    else
        Rin = @(x,transp) spDiff(x,transp,[0 0 1 1 1]);
    end
    x0diff = sqrt(gamRel).*spDiff(y,'notransp',[0 0 1 1 1]);
    resid = AsubAllf(y,'notransp',A,szIm(1:3),freqs,fshells,opt.fthresh,mask) - bin(:,:,:);
    norm2diff  = sum(x0diff(:).^2);
    norm2resid = sum(resid(:).^2);
    clear x0diff resid
    gam = opt.tvReg*gamRel*norm2resid/norm2diff;
    [y, resSqAll2, ~,~,~,~,status2] = cgne(@(x,transp) AsubAllf(x,transp,A,szIm(1:3),freqs,fshells,opt.fthresh,mask),...
        bin,y,NitMax,opt,gam,Rin);
end
clear bin S A* 
y = gather(y);
y = permute(y,[1 3 4 5 2]);
y = y.*reshape(mask, [1 szIm(1:3)]);
Dperp = reshape(y(2,:), [szIm(1:3), length(fshells)]);
Dpar  = reshape(y(3,:), [szIm(1:3), length(fshells)]);
Kperp = reshape(y(4,:), [szIm(1:3), length(fshells)]);
Kpar  = reshape(y(5,:), [szIm(1:3), length(fshells)]);
Wmean = reshape(y(6,:), [szIm(1:3), length(fshells)]);
clear y
bval = bval0;
freqs = freqs0;

% Compute powder average kurtosis
for nf = 1:length(fshells)    
    S_fp = mean(im(:,:,:,bval<opt.bthresh),4);
    for nb = 2:length(bshells)
        inds_b = and( abs(freqs-fshells(nf))<opt.fthresh, ...
            abs(bval-bshells(nb))<opt.bthresh );
        S_fp = cat(4, S_fp, mean(im(:,:,:,inds_b),4));
    end
    S_fp = permute(S_fp,[4 1:3]);
    S_fp = log(S_fp);
    Ap = [ones(length(bshells),1,'like',bshells), -bshells(:), bshells(:).^2/6];
    y = zeros([size(Ap,2), szIm(1:3)],'like',im);
    y(:,mask) = Ap\S_fp(:,mask);
    Dpowder(:,:,:,nf) = gather(reshape(y(2,:), szIm(1:3)));
    Wpowder(:,:,:,nf) = gather(reshape(y(3,:), szIm(1:3)));
end
clear y S_fp

% Final computation. So far all kurtosis parameters have included
% diffusivity, so here we need to divide that out.
Dmean = 1/3*Dpar + 2/3*Dperp;
Kperp = Kperp./(Dperp.^2); 
Kpar = Kpar./(Dpar.^2);
Wmean = Wmean./(Dmean.^2);
Wpowder = Wpowder./(Dpowder.^2);
vec = reshape(permute(vec,[2 1]), [szIm(1:3) 3]);
FA = sqrt(3/2* ((Dpar-Dmean).^2+2*(Dperp-Dmean).^2) ./ (Dpar.^2+2*Dperp.^2) );
FAvec = mean(FA,4).*vec;

% Create NIFTI files. Inherit header information from input NIFTI
if opt.saveNifti
    im_info = niftiinfo(file);
    if length(fshells) == 1 % Case for only one frq, so final images are 3D
        im_info.PixelDimensions = im_info.PixelDimensions(1:3);
        im_info.ImageSize = im_info.ImageSize(1:3);
        im_info.raw.dim(1) = 3;
        im_info.raw.dim(5) = 1;
        im_info.raw.pixdim(5) = 0;
        im_info.raw.dim_info = ' ';
    else
        im_info.ImageSize(4) = length(fshells);
        im_info.raw.dim(5) = length(fshells);
        im_info.raw.dim_info = ' ';
    end

    %
    im_info.Datatype = 'single';
    im_info.BitsPerPixel = 32;
    im_info.raw.datatype = 16;
    im_info.raw.bitpix = 32;
    if isempty(savename)
        [filepath,savename,ext] = fileparts(file);
        if strcmp(ext,'.gz')
            [~,savename,~] = fileparts(savename); 
        end
        savename(savename=='.') = '_'; % niftiwrite does not like periods in name
    end
    
    % Bringing the calculated values to a resonable range
    Dmean(Dmean>10 | Dmean<0) = 0;
    Dpar(Dpar>10 | Dpar<0) = 0;
    Dperp(Dpar>10 | Dperp<0) = 0;
    Wmean(Wmean>10 | Wmean<0) = 0;
    Kpar(Kpar>10 | Kpar<0) = 0;
    Kperp(Kperp>10 | Kperp<0) = 0;
    FA(FA>1.25 | FA<0) = 0;
    Dpowder(Dpowder>10 | Dpowder<0) = 0;
    Wpowder(Wpowder>10 | Wpowder<0) = 0;

    niftiwrite(single(Dmean), sprintf('%s_Dmean', savename), im_info, 'Compressed', true);
    niftiwrite(single(Dpar), sprintf('%s_Dpar', savename), im_info, 'Compressed', true);
    niftiwrite(single(Dperp), sprintf('%s_Dperp', savename), im_info, 'Compressed', true);
    niftiwrite(single(Wmean), sprintf('%s_Wmean', savename), im_info, 'Compressed', true);
    niftiwrite(single(Kpar), sprintf('%s_Kpar', savename), im_info, 'Compressed', true);
    niftiwrite(single(Kperp), sprintf('%s_Kperp', savename), im_info, 'Compressed', true);
    niftiwrite(single(FA), sprintf('%s_FA', savename), im_info, 'Compressed', true);
    niftiwrite(single(Dpowder), sprintf('%s_Dpowder', savename), im_info, 'Compressed', true);
    niftiwrite(single(Wpowder), sprintf('%s_Wpowder', savename), im_info, 'Compressed', true);

    % update 4D nifti info data for FA vec file if multiple freq avialable
    if length(fshells) == 1 % Case for only one frq, so final images are 3D
        im_info = niftiinfo(file);
        im_info.ImageSize(4) = 3;

        im_info.Datatype = 'single';
        im_info.BitsPerPixel = 32;
        im_info.raw.datatype = 16;
        im_info.raw.bitpix = 32;
    else
        im_info.ImageSize(4) = 3;
        im_info.raw.dim(5) = 3;
    end
    
    niftiwrite(single(abs(FAvec)), sprintf('%s_FAvec', savename), im_info, 'Compressed', true);
    dlmwrite(sprintf('%s.fshells', savename),fshells)
end

fprintf('took %d sec\n', round(toc(tic1)));

end

function shellVals = shellSort(listIn, thresh)
% Rough estimate of shells
shellVals = [];
for n=1:length(listIn)
    if isempty(shellVals)
        shellVals = listIn(n);
    elseif ~any(abs(listIn(n)-shellVals) < thresh)
        shellVals = cat(1,shellVals,listIn(n));
    end
end

% Use mean b-value for each shell, then sort the lists
for n=1:length(shellVals)
    inds = abs(listIn-shellVals(n)) < thresh;
    shellVals(n) = mean(listIn(inds));
end
shellVals = sort(shellVals,'ascend');
end


function y = AsubDTI(x,transp,A,szin,mask)
if strcmp(transp,'notransp')
    y = zeros([size(A,1) szin], 'like', A);
    y(:,mask) = A*x(:,mask);
else
    y = zeros([size(A,2) szin], 'like', A);
    y(:,mask) = (A')*x(:,mask);
    y = reshape(y, [size(y,1) szin(1:3)]);
end
end

function y = AsubAllf(x,transp,A,szin,freqs,fshells,fthresh,mask)
if strcmp(transp,'notransp')
    y = zeros([size(A,1) 1 prod(szin)],'like',A);
else
    y = zeros([size(A,2) length(fshells) szin],'like',A);
end
szy = size(y);
for nf = 1:length(fshells)
    inds_f = abs(freqs-fshells(nf))<fthresh ;
    if strcmp(transp,'notransp')
        if isgpuarray(A)
            y(inds_f,:,mask) = pagefun(@mtimes, A(inds_f,:,mask), x(:,nf,mask));
        else
            for n=1:prod(szy(3:end))
                if mask(n)
                    y(inds_f,:,n) = A(inds_f,:,n)*x(:,nf,n);
                end
            end
        end
    elseif strcmp(transp,'transp')
        if isgpuarray(A)
            y_a = pagefun(@ctranspose, A(inds_f,:,mask));
            y(:,nf,mask) = pagefun(@mtimes, y_a, x(inds_f,:,mask));
        else
            for n=1:prod(szy(3:end))
                if mask(n)
                    y(:,nf,n) = A(inds_f,:,n)'*x(inds_f,:,n);
                end
            end
        end
    elseif strcmp(transp,'pinv')
        if isgpuarray(A)
            y(:,nf,mask) = pagefun(@mldivide, A(inds_f,:,mask), x(inds_f,:,mask));
        else
            for n=1:prod(szy(3:end))
                if mask(n)
                    y(:,nf,n) = A(inds_f,:,n)\x(inds_f,:,n);
                end
            end
        end
    end
end
end

function y = RsubFreq(x,transp,fScl)
    if strcmp(transp,'notransp')
        y = cat(7, fScl*spDiff(spDiff(x, transp,[0 1 0 0 0]),transp,[0 0 1 1 1]), ...
            spDiff(x, transp,[0 0 1 1 1]));
    else
        y = fScl*spDiff(spDiff(x(:,:,:,:,:,:,1), transp,[0 0 1 1 1]),transp,[0 1 0 0 0]) + ...
            spDiff(x(:,:,:,:,:,:,2), transp,[0 0 1 1 1]);
    end
end

