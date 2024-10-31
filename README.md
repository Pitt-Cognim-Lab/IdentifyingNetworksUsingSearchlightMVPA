# IdentifyingNetworksUsingSearchlightMVPA
This repo contains the final code for the project "Identifying networks within an fMRI multivariate searchlight analysis"


## Part 1: Searchlight classification.

1. Warpping the functioncal data to the standard space
   Template used: MNI152_T1_2009c+tlrc
   Here is the link to the bash script used for warpping the images: https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/wrapping_images.sh
   First, warped the anatomical image of each subject to the standard image. This will give us the transformation matrix. The warped functional images to the anatomical images. After that used the transformation matrices to warp the functional images to the standard space.
     
2. Running the searchlight classification on these standardized images
   Here is the code for the searchlight classification: https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/searchlight_classification.m
   In this script I am reading the standardized images of all the subjects and running a 4-fold 4-way searchlight classification. This script takes about a day to run, since searchlight classification for each subject takes about an hour. I am filling in the labels between the trs where one label was showed. There was an instance where the voxels in few of the searchglights was empty. There is a piece of code before the classification to find remove those voxels.
   I am saving out the mean accuracies (across the folds) and predicted labels of each searchlight. I am also saving the whole workspace after running the searchlight. There are few variables of interest in this saved workspace. Here are few variables that I used most frequently.
   **individual_results** - This is the binary vector of the accuracies for each searchlight.
   **mean_searchlight_accuracies** - This file saves the mean accuracies for each searchlight. I use this to save it to afni and get significant accuracies.
   **individual_labels** - This a vector of the predicted labels for each of the time points.

   <!---side note: I didnt save them before and wrote another script to get them which is [here](https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/GetPredictedLabels.m)--->
   
4. Getting the significant accuracies from the above classification
   After this I used the mean_searchlight_accuracies_WB+tlrc for each subject and ran t-test to get significant accuracies. The mean searchlight accuracies were basically the mean of accuracies for each fold for the searchlight centered at that particular voxel. 

   The ttest command that I ran:

   3dttest++ -prefix ttest.results/ttest_results -mask ../Functional_data/Whole_brain_mask_stan_re+tlrc -singletonA 0.25 -setB group_analysis 01 mean_searchlight_results_WBs103+tlrc 02 mean_searchlight_results_WBs105+tlrc 03 mean_searchlight_results_WBs107+tlrc 04 mean_searchlight_results_WBs108+tlrc 05 mean_searchlight_results_WBs109+tlrc 06 mean_searchlight_results_WBs110+tlrc 07 mean_searchlight_results_WBs112+tlrc 08 mean_searchlight_results_WBs113+tlrc 09 mean_searchlight_results_WBs114+tlrc 10 mean_searchlight_results_WBs115+tlrc 11 mean_searchlight_results_WBs116+tlrc 12 mean_searchlight_results_WBs117+tlrc 13 mean_searchlight_results_WBs118+tlrc 14 mean_searchlight_results_WBs119+tlrc 15 mean_searchlight_results_WBs120+tlrc 16 mean_searchlight_results_WBs121+tlrc 17 mean_searchlight_results_WBs122+tlrc 18 mean_searchlight_results_WBs123+tlrc 19 mean_searchlight_results_WBs126+tlrc 20 mean_searchlight_results_WBs128+tlrc
   
   To get the mask for p<0.05 threshold. I ran 3dcalc (2.0930 was the threshold value of the voxel with p<0.05, I think it was a t-value map.) 

   3dcalc -a 'ttest_results+tlrc[1]' -expr 'or(isnegative(a+2.0930) ,ispositive(a-2.0930) )' -prefix thr05_mask

   Then ran the get_significant_accuracies script to get the voxels of significant accuracies and then run the ICA on them.

   t-test result: without cluster correction:
   <img src="https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/Images/ttest_axial-without_bg.png" width="30%" height="30%">
   <img src="https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/Images/ttest_coronal-without_bg.png" width="30%" height="30%">
   <img src="https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/Images/ttest_sag-without_bg.png" width="30%" height="30%">

## Part 2: Secondary Analysis

### 1. ICA
   After getting the searchlights that were significant I put them through ICA classification. Here is the link to the python script for [ICA](https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/ICA.ipynb). In this script I read each subjects significant searchlight accuracies and run ICA on it. Details are as follows:

   - I standardize the data where it makes the mean zero and variance 1. I dont think it changes the covariance of the column. 
      
   - I run PCA to just get the number of components, but not as a data reduction step. I look at the point where the variance of the data stops changing significantly. It looks like after 7 components it reduces very less. The total variance of the top 7 components is less than 20%. It saves a lot less variance in the data. 
      
   - After getting the components I threshold the values to get which voxels were significant part of the component. I did that by getting the mean and standard deviation. The voxels that were included had value that were 3 standard deviation far from mean. [ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5594474/]
      
   - I get the components from here and save them as csv files. The I use [general_write2afni](https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/general_write2afni.m) script to save these components in afni (PS: I make changes in it according to the requirement, so double check before running).
      
   - After that to make sense of the components I calculated the confusion matrices for each component using this [code](https://github.com/Pitt-Cognim-Lab/IdentifyingNetworksUsingSearchlightMVPA/blob/main/Component_cluster_confusion_matrix.m). Here are the results for that. There are 7 components right now.
     <!--The components are present on the server at "smb://data.lrdcfile.pitt.edu/project/Coutanche/Shared/Projects/Searchlight/Component_maps/"-->


-----------------------------------------updated till here-------------------------------------------------------------------


Confusion Matrices for all the clusters in each component:
(Component had no significant voxels)

|        |  Word    |   Face  |  Tool  |  Number  |
|--------|----------|---------|--------|----------|
| Comp1: |  NaN     |   NaN   |  NaN   |   NaN    |
| Comp2: |  0.2693  |  0.2743 | 0.2849 |  0.4441  |
| Comp3: |  0.3013  |  0.4353 | 0.3127 |  0.2960  |
| Comp4: |  0.3003  |  0.3412 | 0.3024 |  0.2827  |
| Comp5: |  0.2464  |  0.2959 | 0.2679 |  0.3315  |
| Comp6: |  0.3264  |  0.5325 | 0.2908 |  0.2812  |
| Comp7: |  0.3680  |  0.3705 | 0.2937 |  0.2822  |

**Component: 2** [255]    0.3167

<img src="https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/Images/7ComponentsC2.png" width="70%" height="70%">

Cluster: 1 [184]  0.3201

Brain areas  center of mass (28.5, -72.2, 48.9): Parietal (superior, intra, posterior), recognition memory, reasoning, rule, tasks, retrieval, spatial, abilities, abstract, abuse. 

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 0.2823  |  0.2628  |  0.2885 |   0.4469 |

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix2_1.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 11.2916  |  9.8269  |  9.5497  |  9.3318 | 
| 10.9095  | 10.5136  |  9.7935  |  8.7834 |
| 9.5576   | 9.9361  | 11.5416   | 8.9647 |
| 7.2970  |  7.3549  |  7.4736  | 17.8745 |

Cluster: 2  [70] 0.3132

Brain areas CoM (-29.3,-73.7, 47.1) : Switching, calculation, parietal (intra, superior), naming, retrieval, angular, learn, items, angular gyrus, memory, abilities, abstract, abuse.

 | Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 0.2357  |  0.3038  |  0.2763  |  0.4372 |

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix2_2.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 9.4286 |  10.0786 |  10.1364  | 10.3564 |
| 9.9129 |  12.1529 |   9.9971  |  7.9371 |
| 9.4750 |   9.7643 |  11.0507  |  9.7100 |
| 7.3936 |   6.9736  |  8.1450  | 17.4879 |

**Component: 3**   [279]     0.3363

<img src="https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/Images/7ComponentsC3.png" width="70%" height="70%">

Cluster: 1   [279]   0.3363

Brain areas CoM (-42.9, -57.3, -19.2): Fusiform, visual word, word form, orthographic, fusiform (gyrus, face), reading, ffa(x3), face, word (x2), written, readers, fusiform, semantic.

 | Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 0.3013  |  0.4353  |  0.3127 |   0.2960 |

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix3_1.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 12.0520 |   8.0514  | 10.1367  |  9.7599 |
| 7.7627  | 17.4100  |  7.4065  |  7.4208 |
| 8.9930  |  7.5873 |  12.5097 |  10.9100 |
| 9.6961  |  7.3047 |  11.1611 |  11.8382 |

**Component: 4**  [274]    0.3065

<img src="https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/Images/7ComponentsC4.png" width="70%" height="70%">

Cluster: 1    [273]   0.3065

Brain areas CoM (3.1, -89.7, -4.2): v1, visual, risk taking, mt, medial pfc, negative emotions, asd,coordination, risky, v5, occipital, abilities (x2), abstract, abuse, acc, accumbens, accurate (x2), acoustic. 

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 0.3002  |  0.3409 |   0.3022 |   0.2827 |

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix4_1.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 12.0073  |  8.7408 |  10.4874  |  8.7645 |
| 8.3963 |  13.6353  |  8.7606  |  9.2077 |
| 9.7678  |  8.2288 |  12.0899  |  9.9136 |
| 8.3985  |  9.8158 |  10.4762 |  11.3095 |

**Component: 5**    [114]    0.2865

(This just shows the largest cluster)

<img src="https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/Images/7ComponentsC5.png" width="80%" height="80%">

Cluster: 1     [51]      0.2850

Brain areas CoM (2.2, 10.4, 49.6): No data available. 

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 0.2221 |   0.2407   | 0.2927  |  0.3844|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix5_1.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 8.8843  | 10.4206 |  10.4422 |  10.2529 |
| 9.7569  |  9.6265 |  10.5931 |  10.0235 |
| 9.2657  | 10.1961 |  11.7088  |  8.8294 |
| 8.2186  |  8.1245 |   8.2794 |  15.3775 |

Cluster: 2     [25]    0.2811

Brain area CoM (-1.8, -61.9, 27.7): Autobiographical, posterior cingulate, precuneus, autobiographical memory, default, precuneus posterior, theory mind, default mode, social, episodic, mental states, memory retrieval, theory mentalizing, cingulate, posterior, mind. 

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
|0.2619  |  0.3290 |   0.2302  |  0.3031|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix5_2.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 10.4760  |  8.5440 |  10.9900  |  9.9900 |
| 8.6200 |  13.1600  |  8.4160  |  9.8040 |
| 10.5300  |  9.3020  |  9.2100 |  10.9580 |
| 9.4600  |  8.6280  |  9.7860 |  12.1260 |

Cluster: 3    [21]    0.2859

Brain area CoM (-45.9, -75.9, 26.7): autobiographical (memory), default mode (x2), episodic, angular, memory, angular gyrus, semantic (memory), retrieval, memory retrieval, retrosplenial, names, dorsal attention, remembering, thinking.  

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
|0.2688  |  0.3275  |  0.2612  |  0.2860|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix5_3.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 10.7524  |  9.3000 |  10.1429  |  9.8048 |
| 9.2238  | 13.1000  |  8.7167  |  8.9595 |
| 9.5881 |  10.2952  | 10.4476  |  9.6690 |
| 9.1786  |  9.3548 |  10.0286  | 11.4381 |

Cluster: 4     [11]         0.2940

Brain area CoM (53.8, -63.7, 19.6): temporoparietal junction (x4), theory mind, mind tom, amnestic, imagine, temoporo parietal, default, junction (x3), posterior cingulate, temoporo, memories

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
|0.2691  |  0.3998  | 0.2566   | 0.2505|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix5_4.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 10.7636  |  8.2727 |  11.0091  |  9.9545 |
| 7.9318 |  15.9909  |  7.2318  |  8.8455 |
| 10.8318 |   8.8227 |  10.2636 |  10.0818 |
| 9.8500  |  9.8636 |  10.2682 |  10.0182 |

**Component: 6**   [343]      0.3577

<img src="https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/Images/7ComponentsC6.png" width="80%" height="80%">

Cluster: 1     [343]       0.3577

Brain area CoM (40.1, -61.8, -17.3): fusiform (x5), face (x2), recognition, objects, occipito-temoporal, recognize, visual, object, occipito, face recognition, letters, characters, occipital (inferior). 

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
|0.3264  |  0.5325  |  0.2908  |  0.2812|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix6_1.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 13.0567  |  6.1529  | 10.5586 |  10.2318 |
| 6.1829 |  21.3006 |  6.2882  |  6.2283 |
| 11.4872  |  6.9423 |  11.6306  |  9.9399 |
| 10.7668  |  7.4854 |  10.5010  | 11.2468 |

**Component: 7**    [286]      0.3284

<img src="https://github.com/Pitt-Cognim-Lab/Searchlight/blob/main/Images/7ComponentsC7.png" width="80%" height="80%">

Cluster: 1     [147]         0.3357

Brain area CoM (+33.6, -88.5, 9.1): visual, occipital (x2), occipitotemporal, v5, encoding, fusiform gyri, ventral visual, perceptual, ffa, occipital temporal, expertise, characters, lateral occipital, face, visual stream, ffa, gain, abilities, ability. 

| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
|0.3754  |  0.3849  |  0.2955  |  0.2870|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix7_1.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 15.0146  |  7.6861  |  9.7173  |  7.5820 |
| 8.2796 |  15.3949  |  7.5354  |  8.7901 |
| 10.2929 |   7.8915  | 11.8218  |  9.9939 |
| 9.1102  |  9.3163 |  10.0918 |  11.4816 |

Cluster: 2    [139]          0.3211

Brain area CoM (-30.3, -94.1, 10.4):  occipital, visual, occipital cortex, fusiform (x2), middle occipital, lateral occipital, rotation, visual stream, visual perception, occipitotemporal, abilities 9x2), abstract, abuse, acc, accumbens, accurate (x2), acoustic. 

 | Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
|0.3603  |  0.3553  |  0.2917  |  0.2770|

<img src="https://github.com/Pitt-Cognim-Lab/Lab_Notebooks/blob/main/To_add_image/Searchlight/ConfusionMatrix7_2.png" width="40%" height="40%">

Averaged over:
| Word | Face  | Tool | Number | 
|-------|----------|--------|---------|
| 14.4119  |  7.9784  |  9.5953 |  8.0144 |
| 9.1241 |  14.2137  |  8.2104  |  8.4518 |
| 10.2687  |  8.0273 |  11.6694 |  10.0345 |
| 8.6255  |  9.8169  | 10.4763 |  11.0813 |


Paper that talks about how to threshold the components : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5594474/

Few resources talking about the negative weights of the ICA components:
   - https://www.brainvoyager.com/ubb/Forum4/HTML/000613.html
   - https://www.researchgate.net/post/What_does_negativity_mean_in_the_unmixed_components_in_Independent_Component_Analysis_ICA
