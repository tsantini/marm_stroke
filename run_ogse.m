clearvars
addpath matmri
setPath

files = dir('E:\250425-processing_marm_stroke\2022_marmStroke\all_OGSE_uFA\BM4D_PCAMoeda_20230116*OGSEa.nii');

for i=1:length(files)
    tic
    cd (files(i).folder)
    opt.noGPU = 0;
    %opt.tvRegVec = 0;
    %opt.tvReg = 0;
    % img = niftiread([files(i).folder filesep erase(files(i).name,'.nii')]);
    % open_nii_TS(img(:,:,:,1))
	bval=load([files(i).folder filesep erase(files(i).name,'.nii') '.bval']);
	freqs = repmat([0 57 33 12],1, length(bval)/4);
    [Dmean,Dpar,Dperp,Wmean,Kpar,Kperp,FA,fshells,FAvec,Dpowder,Wpowder] = nii2kurt([files(i).folder filesep files(i).name], [files(i).folder filesep erase(files(i).name,'.nii')], freqs, [], opt, [files(i).folder filesep erase(erase(files(i).name,'.gz'),'.nii')]);
    toc
end



