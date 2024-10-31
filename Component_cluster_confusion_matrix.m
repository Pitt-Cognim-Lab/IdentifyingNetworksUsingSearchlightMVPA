%% This script is to get confusion matrix for each cluster in each component. 
clear;
subject_dirs=["s103","s105","s107","s108","s109","s110","s112","s113","s114","s115","s116","s117","s118","s119","s120","s121","s122","s123","s126","s128"];

% average_subject=zeros(61910,160); %% Hardcoding the indices here.
% In this label 1 will be saved from 1:40, label 2 in 41:80, label 3 in
% 81:120 and label 4 in 121:160
cd ("subjects_data/Functional_data")
%% First lets get the average accuracy across subjects.
numVoxels=61910; %% Hardcoding it here, cause I dont want to think about it.
sum_confusion_matrices=zeros(numVoxels,4,4); % Hardcoding cause I know that there are 4 classes.
avg_confusion_matrices=zeros(numVoxels,4,4);
for s=1:length(subject_dirs)
    subject=subject_dirs(s);
    cd (subject);
    fname=sprintf('Searchlight_classification_results_%s.mat',subject);
    load(fname,'conditions_matrix');
    load('individual_labels.mat');
    % Getting the appropriate label;
    labels=getLabels(conditions_matrix);
    labels=labels(find(labels)); % True labels
    numVoxels=length(individual_labels);
    confusionMatrices=zeros(numVoxels,4,4);
    for v=1:numVoxels
        predictedLabs=individual_labels(v,:); % Predicted labels
        confusionMatrices(v,:,:)=confusionmat(labels,predictedLabs); % this will give me the confusion matrix. 
        % x-label - predicted class & y - true classes
        % y-axis - predicted labels & z-axis - true labels
    end
    sum_confusion_matrices=sum_confusion_matrices+confusionMatrices;
    % For each searchlight I will get the confusion matrix and then add it
    % all up for the whole cluster.
    % for l=1:4
    %     lab_inds=find(labels==l);
    %     predLabs=predictedLabs(:,lab_inds);
    %     % average_subject(:,(l-1)*40+1:l*40)=(average_subject(:,(l-1)*40+1:l*40)*(s-1)+results_sub)/s; % accuracy score
    % end
    cd ../
end
avg_confusion_matrices=sum_confusion_matrices/20;

% Now, we have the average accuracies of each voxel at each time point for
% each label across all the subjects. 
clear subj_mask;
%% Getting average true positives, true negatives, false positives and false negatives 
% % for each cluster of each component
% Getting the mask for flags;
subj_mask = init_subj('Components','comp_mask');
subj_mask=load_afni_mask(subj_mask,'roi_mask','Whole_brain_mask_stan_re+tlrc');
% cd ../../Component_maps/Individual_voxelValue_components/
cd ../../Component_maps/Comps_masks_2std %% This will give the confusion matrix for the 2 std ICA components
% cd ../../Univariate_masks/facevsall/comps_mask/ %% Masks for facevsall components
num_comps=7;
for c=1:num_comps
    cname=sprintf('comp%d_2std_mask+tlrc',c); %% Have to change the name here as well.
    % cname=sprintf('comps_mask_facevsall_8c2std_C%d_mask+tlrc',c);
    pattern_name=sprintf('component%d',c);
    subj_mask=load_afni_pattern(subj_mask,pattern_name,'roi_mask',cname);
    flags=get_mat(subj_mask,'pattern',pattern_name);
    numClusters=max(unique(flags));
    fprintf("Component: %d\n",c);
    component_accuracy=[];
    for clus=1:numClusters
        fprintf("Cluster: %d\n",clus);
        % Have to get the labels that belong to that cluster
        clusterLabelsInds=find(flags==clus);
        disp("Cluster Size: ");
        disp(length(clusterLabelsInds));
        temp_cm=sum_confusion_matrices(clusterLabelsInds,:,:);
        temp_avg_cm=avg_confusion_matrices(clusterLabelsInds,:,:);
        final_sum_cm=sum(temp_cm);
        final_sum_cm=reshape(final_sum_cm,[4,4]);
        final_avg_cm=mean(temp_avg_cm);
        final_avg_cm=reshape(final_avg_cm,[4,4]);
        cm_percent=final_sum_cm./sum(sum(final_sum_cm));
        cm_percent=cm_percent*100;
        accuracy_final=(final_sum_cm(1,1)+final_sum_cm(2,2)+final_sum_cm(3,3)+final_sum_cm(4,4))/sum(sum(final_sum_cm));
        component_accuracy=[component_accuracy,accuracy_final];
        disp("Accuracy of the cluster: ");
        disp(accuracy_final);
        disp("Summed over: ")
        disp(final_sum_cm)
        % cm=confusionchart(final_sum_cm);
        h=heatmap(cm_percent);
        h.XLabel="Predicted Class";
        h.YLabel="True Class";
        chart_name=sprintf('ConfusionMatrix%%%d_%d.png',c,clus);
        saveas(h,chart_name);
        disp("Averaged over: ")
        disp(final_avg_cm)
    end
    disp("Component accuracy: ");
    mean(component_accuracy)

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