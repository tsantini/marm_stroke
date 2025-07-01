data=/mnt/e/250425-processing_marm_stroke/2022_marmStroke/all_OGSE_uFA

for files in $data/*_LTE_all.nii.gz; do
	file_name=`basename $files`
	echo processing $file_name
	docker run -v $data:/data -e file_name=$file_name -w /data nyudiffusionmri/designer2:main bash -c 'tmi /data/${file_name},/data/${file_name%_LTE_all.nii.gz}_STE_all.nii.gz -bshape 1,0 -DTI /data/${file_name}_designer_DTI'
	docker run -v $data:/data -e file_name=$file_name -w /data nyudiffusionmri/designer2:main bash -c 'tmi /data/${file_name},/data/${file_name%_LTE_all.nii.gz}_STE_all.nii.gz -bshape 1,0 -DKI -WDKI -akc_outliers -fit_smoothing 25 /data/${file_name}_designer_DKI'
	docker run -v $data:/data -e file_name=$file_name -w /data nyudiffusionmri/designer2:main bash -c 'tmi /data/${file_name},/data/${file_name%_LTE_all.nii.gz}_STE_all.nii.gz -bshape 1,0 -SMI -compartments EAS,IAS,FW  /data/${file_name}_designer_SMI'
done



