## Working in progress
## TODO: process DTI, DKI, SM, NODDI

data=/mnt/e/250425-processing_marm_stroke/2022_marmStroke/all_OGSE_uFA
docker run -v $data:/data -w /data nyudiffusionmri/designer2:main bash -c 'designer'
docker run -v $data:/data -w /data nyudiffusionmri/designer2:main bash -c 'tmi'