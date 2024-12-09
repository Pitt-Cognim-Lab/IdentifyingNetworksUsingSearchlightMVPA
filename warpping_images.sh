#!/bin/bash

## Code to run the wrapping for all the subjects
#Usage ./wrapping_images $subj_num
# subj=$1
subjects=("s103" "s105" "s107" "s108" "s109" "s110" "s112" "s113" "s114" "s115" "s116" "s117" "s118" "s119" "s120" "s121" "s122" "s123" "s126" "s128")

for subj in "${subjects[@]}"
do
	echo $subj
	cd $subj

	# First wrapping the anatomical image to the standard space. This will retain the good spatial resolution.
	auto_warp.py -base MNI152_T1_2009c+tlrc -input 1_MPRAGE_ns+orig

	# Resampling the functional images to the anatomical resolution so that I can use the matrices that I got from above command

	3dresample -prefix pb05.$subj.r02.empty_re+orig -input pb05.$subj.r02.empty+orig -master 1_MPRAGE_ns+orig

	3dresample -prefix pb05.$subj.r04.empty_re+orig -input pb05.$subj.r04.empty+orig -master 1_MPRAGE_ns+orig

	3dresample -prefix pb05.$subj.r06.empty_re+orig -input pb05.$subj.r06.empty+orig -master 1_MPRAGE_ns+orig

	3dresample -prefix pb05.$subj.r08.empty_re+orig -input pb05.$subj.r08.empty+orig -master 1_MPRAGE_ns+orig


	## Wrapping these functional images to the standard space- this will shift the space as well to 3.125 mm space as well. 
	3dNwarpApply -master ../MNI152_T1_2009c+tlrc -dxyz 3.125 -source pb05.$subj.r02.empty_re+orig pb05.$subj.r04.empty_re+orig pb05.$subj.r06.empty_re+orig pb05.$subj.r08.empty_re+orig -nwarp 'awpy/anat.un.aff.qw_WARP.nii awpy/anat.un.aff.Xat.1D' -suffix _stan

	cd ..
done
