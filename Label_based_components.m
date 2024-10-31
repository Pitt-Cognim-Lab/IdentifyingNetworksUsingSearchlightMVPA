%% This script is to figure out which components are relevant to which label. 
% % I am taking in the components and getting the corresponding accuracy vectors for each searchlight
% % and calculating that accuracy. 
subject_dirs=["s103","s105","s107","s108","s109","s110","s112","s113","s114","s115","s116","s117","s118","s119","s120","s121","s122","s123","s126","s128"];

% Stay in the main directory
% % first have to read in each component as a mask
% This will tell which searchlight belongs to which component
% comps=readmatrix("../../Components_map_concatenated_7c.csv"); 
% comps=readmatrix("Components_map_7c_2std.csv"); % this will have 7 columns
% comps = readmatrix("Components_map_7c_2std.csv");
comps=readmatrix("Univariate_masks/Components_map_face_8c_1_2std.csv");

% Getting the mask for flags;
subj_mask = init_subj('Flags','ttest_mask');
subj_mask=load_afni_mask(subj_mask,'roi_mask','Component_maps/Whole_brain_mask_stan_re+tlrc');
% Hardcoded the path here
% sig_acc='subjects_data/Mean_searchlight_results/ttest.results/thr05_mask+tlrc';
sig_acc='Univariate_masks/stats_face/ttest.results/thrp005_a05_1_mask+tlrc';
subj_mask = load_afni_pattern(subj_mask,'flag','roi_mask',sig_acc);
flag_mat=get_mat(subj_mask,'pattern','flag');

spatial_map=zeros(length(flag_mat),size(comps,2));
ind=1;
for i=1:length(flag_mat)
    if flag_mat(i)~=0
        spatial_map(i,:)=comps(ind,:);
        ind=ind+1;
    end
end

% Sense check
for c=1:7
    if ~(length(find(spatial_map(:,1)))==length(find(comps(:,1))))
        disp("The number of voxels in comps and spatial map are not matching for component")
        disp(c)
    end
end

cd 'subjects_data/Functional_data/'
average_subj_accuracies=zeros(length(flag_mat),4); % this 4 is for 4 labels. 
% Taking a running average of the accuracies for each label for each
% participant.
for s=1:length(subject_dirs)
    subject=subject_dirs(s);
    cd (subject);
    fname=sprintf('Searchlight_classification_results_%s.mat',subject);
    load(fname,'conditions_matrix','individual_results');
    % I will need regressors and individual accuracies - conditions matrix
    % just have to shrink it.
    labels=getLabels(conditions_matrix);
    labels=labels(find(labels));
    % This loop is for each label
    for i=1:4
        indices=find(labels==i);
        results=individual_results(:,indices);
        av_results=mean(results,2); % average over the label
        % Taking a running average for all the subjects
        average_subj_accuracies(:,i)=(average_subj_accuracies(:,i)*(s-1) + av_results)/s;
    end
    cd ..
end
cd ../../
% %% Lets save what we need and then load it back in
% save("Average_subject_accuracies_7c_2std.mat","average_subj_accuracies","spatial_map");
% clear;

%% Now I have my component maps and I have the accuracies as well. I will run it
% % For each component it is finding the accuracy of each label in each
% searchlight belonging to that component. 
% load("Average_subject_accuracies_7c_2std.mat");
[num_voxels,gg]=size(spatial_map);
% num_comps=max(max(spatial_map));
num_comps=gg;
accuracies_maps=zeros(num_voxels,num_comps,4); 
for c=1:num_comps
    voxel_indices=find(spatial_map(:,c)==c); % this will give me the indices of all the voxels in component c
    accuracies_maps(voxel_indices,c,:)=average_subj_accuracies(voxel_indices,:);
end

% getting average of component accuracy for each label
for c=1:num_comps
    temp_label_acc=zeros(1,4);
    for a=1:4
        comp_acc_map=accuracies_maps(find(accuracies_maps(:,c,a)),c,a);
        % length(comp_acc_map)
        temp_label_acc(1,a)=mean(comp_acc_map);
    end
    disp(temp_label_acc);
end

%% Getting the accuracy for each cluster of each component for each label
subj_mask = init_subj('Components','comp_mask');
subj_mask=load_afni_mask(subj_mask,'roi_mask','Component_maps/Whole_brain_mask_stan_re+tlrc');
% cd Component_maps/Comps_masks_2std/
cd Univariate_masks/facevsall/comps_mask/
for c=1:num_comps
    % cname=sprintf('comp%d_2std_mask+tlrc',c);
    cname=sprintf('comps_mask_facevsall_8c2std_C%d_mask+tlrc',c);
    pattern_name=sprintf('component%d',c);
    subj_mask=load_afni_pattern(subj_mask,pattern_name,'roi_mask',cname);
    flags=get_mat(subj_mask,'pattern',pattern_name);
    numClusters=max(unique(flags));
    fprintf("Component: %d\n",c);
    component_accuracy=[];
    for clus=1:numClusters
        fprintf("Cluster: %d\n",clus);
        clusterLabelsInds=find(flags==clus);
        % Have to get the labels that belong to that cluster
        temp_label_acc=zeros(1,4);
        for a=1:4
            comp_acc_map=accuracies_maps(clusterLabelsInds,c,a);
            % length(comp_acc_map)
            temp_label_acc(1,a)=mean(comp_acc_map);
        end
    disp(temp_label_acc);
    end
end

%Creating labels
function labels_temp = getLabels(labels)
%This is assuming that labels are one-hot encoded. Columns are the
%timepoints and rows are the different label rows.
    [numLabels ,numTimepoints] = size(labels);
    labels_temp=zeros(1,numTimepoints);
    for j=1:numTimepoints
        ind=find(labels(:,j));
        if isempty(ind)
            continue;
        end
        labels_temp(1,j)=ind;
    end
end