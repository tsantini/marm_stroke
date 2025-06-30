clc; clearvars

addpath matmri
setPath

files = dir('E:\250425-processing_marm_stroke\2022_marmStroke\all_OGSE_uFA\BM4D_PCAMoeda_*uFAa.nii');

for i = 1:length(files)
    bmat=load([files(i).folder filesep erase(files(i).name,'.nii') '.bmat']);
    bval=load([files(i).folder filesep erase(files(i).name,'.nii') '.bval']);
    isiso=std([bmat(1,:); bmat(5,:); bmat(9,:)],[],1)./mean([bmat(1,:); bmat(5,:); bmat(9,:)],1)<0.1 & ([bmat(1,:)+ bmat(5,:)+ bmat(9,:)] > 500);

    info = niftiinfo([files(i).folder filesep files(i).name]);

    repetition=info.ImageSize(4)./size(bmat,2);

	[uFA,Kiso,Klin,D,sf,uA2, uFA_noFWE,Klin_noFWE,Kiso_noFWE,D_noFWE,sSTE,sLTE] = nii2uFA_fwe([files(i).folder filesep files(i).name], [files(i).folder filesep  [erase(files(i).name,'.nii') '.bval']], ...
	 	repmat(isiso,1, repetition),[],[],[],[files(i).folder filesep erase(erase(files(i).name,'.gz'),'.nii')]);
end
