antsRegistrationSyN.sh -d 3 -m ~/fsl_605/data/atlases/JHU/JHU-ICBM-FA-2mm.nii.gz -f mprage.nii.gz -n 20 -o reg &&

antsApplyTransforms -d 3 -i /home/local/VANDERBILT/kanakap/fsl_605/data/atlases/JHU/JHU-ICBM-labels-2mm.nii -r mprage.nii.gz -n NearestNeighbor -t reg1Warp.nii.gz -t [reg0GenericAffine.mat] -o T1atlas2subj.nii.gz
