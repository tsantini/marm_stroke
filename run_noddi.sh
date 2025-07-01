## Working in progress
## TODO: process DTI, DKI, SM, NODDI

data=/mnt/e/250425-processing_marm_stroke/2022_marmStroke/all_OGSE_uFA

for files in $data/*_LTE_all.nii.gz; do
	file_name=`basename $files`
	echo processing $file_name
	docker run -v $data:/data -e file_name=$file_name -w /data tsantini/dwipreproc -w /data/${file_name} -ndky
done
