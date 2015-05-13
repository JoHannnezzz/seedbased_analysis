#!/bin/bash
# to run: 
# ./x.grouplevel_analysis.sh -2 5 23 [/path/to/results/dir] [/path/to/group/level/directory]
# where -2 5 23 are the coordinates of the seed
# and glm.con and glm.mat are the outputs of the fsl group-level model from the fsl gui with 1's in the first column of the EVs and de-meaned covariates of interest in the second column of the EVs
# the contrasts should be set as:
	# 1 0
	# -1 0 
	# 0 1
	# 0 -1
results_dir=${4}
group_dir=${5}

fslmerge -t ${group_dir}/merged_file_roi_${1}.${2}.${3}.nii.gz ${results_dir}/mean_conn_maps/_roi_${1}.${2}.${3}/_subject_id_*/corr_map_smoothed_merged_mean.nii.gz 

fslmaths ${group_dir}/merged_file_roi_${1}.${2}.${3}.nii.gz -abs -Tmin -bin ${group_dir}/mask_roi_${1}.${2}.${3}.nii.gz

randomise -i ${group_dir}/merged_file_roi_${1}.${2}.${3}.nii.gz -o ${group_dir}/results_roi_${1}.${2}.${3} \
	-m ${group_dir}/mask_roi_${1}.${2}.${3}.nii.gz \
  -d ${group_dir}/glm.mat -t ${group_dir}/glm.con -T

fslmaths ${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tfce_corrp_tstat3.nii.gz \
	-thr 0.95 -bin -mul \
  ${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tstat3.nii.gz \
  ${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tfce_corrp_tstat3_threshed.nii.gz

fslmaths ${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tfce_corrp_tstat4.nii.gz \
	-thr 0.95 -bin -mul \
  ${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tstat4.nii.gz \
  ${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tfce_corrp_tstat4_threshed.nii.gz

cluster --in=${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tfce_corrp_tstat3_threshed.nii.gz --thresh=0.0001 --oindex=${group_dir}/results_roi_${1}.${2}.${3}/cluster_positive/results_roi_${1}.${2}.${3}_tfce_corrp_tstat3_threshed_clusterIndex --olmax=${group_dir}/results_roi_${1}.${2}.${3}/cluster_positive/results_roi_${1}.${2}.${3}_tfce_corrp_tstat3_threshed_clusterIndexLmax.txt --osize=${group_dir}/results_roi_${1}.${2}.${3}/cluster_positive/results_roi_${1}.${2}.${3}_tfce_corrp_tstat3_threshed_size -mm > ${group_dir}/results_roi_${1}.${2}.${3}/cluster_positive/cluster_info_tstst3.txt

cluster --in=${group_dir}/results_roi_${1}.${2}.${3}/results_roi_${1}.${2}.${3}_tfce_corrp_tstat4_threshed.nii.gz --thresh=0.0001 --oindex=${group_dir}/results_roi_${1}.${2}.${3}/cluster_negative/results_roi_${1}.${2}.${3}_tfce_corrp_tstat4_threshed_clusterIndex --olmax=${group_dir}/results_roi_${1}.${2}.${3}/cluster_negative/results_roi_${1}.${2}.${3}_tfce_corrp_tstat4_threshed_clusterIndexLmax.txt --osize=${group_dir}/results_roi_${1}.${2}.${3}/cluster_negative/results_roi_${1}.${2}.${3}_tfce_corrp_tstat4_threshed_size -mm > ${group_dir}/results_roi_${1}.${2}.${3}/cluster_negative/cluster_info_tstst4.txt

