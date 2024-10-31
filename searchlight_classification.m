%% Running through all the subjects
% Initializing global variables
subjects=["s103", "s105", "s107","s108","s109","s110","s112","s113","s114","s115","s116","s117","s118","s119","s120","s121","s122","s123","s126","s128"];
num_runs=4;
% Be present in the Functional data directory

%% Normal classification of a spherical ROI
for itr=1:length(subjects)
    subj=subjects(itr);

    subj_sl = init_subj('classification',subj);
    
    %%% create the mask that will be used when loading in the data
    % subj_sl=load_afni_mask(subj_sl,'roi_mask','VT_Mask_final_re+orig');
    subj_sl=load_afni_mask(subj_sl,'roi_mask','Whole_brain_mask_stan_re+tlrc'); % Common Whole brain mask
    % Then go into the subject dir since mask is outside. 
    cd(subj)
    
    % now, read and set up the actual data. load_AFNI_pattern reads in the
    % EPI data from a BRIK file, keeping only the voxels active in the
    % mask (see above)
    raw_fns={};
    for rr=2:2:num_runs*2
        run=strcat('pb05.',subj,'.r0',int2str(rr),'.empty_re_stan+tlrc'); % changed the name for standard space
        raw_fns=[raw_fns,run];
    end
    raw_fns=cellstr(raw_fns);

    subj_sl = load_afni_pattern(subj_sl,'epi','roi_mask',raw_fns);
    
    subj_sl = init_object(subj_sl,'regressors','conditions');
    regressor_fn=sprintf('../../Regressors/%s/Cantlon_Binary_Regressor_4way_concat',subj);
    load(regressor_fn);
    subj_sl = set_mat(subj_sl,'regressors','conditions', Cantlon_Binary_Regressor_4way_concat);
    condnames = {'word', 'face','tool','number'}; 
    subj_sl = set_objfield(subj_sl,'regressors','conditions','condnames',condnames);
    % Hard coded the length of runs here.
    runs=[ones(1,80) ones(1,80).*2 ones(1,80).*3 ones(1,80).*4]; % we know that the runs are 80 timepoints long. 
    subj_sl=init_object(subj_sl,'selector','runs');
    subj_sl=set_mat(subj_sl,'selector','runs',runs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRE-PROCESSING - z-scoring in time and no-peeking anova
    
    % Shifting the regressors by 2 TR
    subj_sl=shift_regressors(subj_sl,'conditions','runs',2);
    % we want to z-score the EPI data (called 'epi'),
    % individually on each run (using the 'runs' selectors). This is supposed
    % to be within each run. treats each run as a separate timecourse
    subj_sl = zscore_runs(subj_sl,'epi','runs');
    
    subj_sl = create_xvalid_indices(subj_sl,'runs');
    
    conds_mat=get_mat(subj_sl,'regressors','conditions_sh2');
    
    % Filling in the labels here - hard coded here as well
    for i=1:4
        for j=2:319
            if(conds_mat(i,j-1)==1) && (conds_mat(i,j+1)==1)
                conds_mat(i,j)=1;
            end
        end
    end
    
    subj_sl=set_mat(subj_sl,'regressors','conditions_sh2',conds_mat);
    
    sum_cond_mat=sum(conds_mat);
    
    % setting up the run splits - so that I dont take useless information.
    for r=1:4
        runs_xval{r}=get_mat(subj_sl,'selector',sprintf('runs_xval_%d',r));
        runs_xval{r}(sum_cond_mat==0)=0;
        subj_sl=set_mat(subj_sl,'selector',sprintf('runs_xval_%d',r),runs_xval{r});
    end
    
    tic;
    %% Running the classifier
    data=get_mat(subj_sl,'pattern','epi_z');
    accuracies4way=zeros(4,1);
    subj_sl.adj_sphere = create_adj_list(subj_sl,'roi_mask','radius',3);
    flag=0;
    numVoxels=size(data,1);
    Timepoints=size(data,2);
    % Storing the accuracies for each cross validation for each searchlight. 
    numRandomization=10000;
    searchlight_accuracies=zeros(length(subj_sl.adj_sphere),4);
    mean_searchlight_accuracies=zeros(length(subj_sl.adj_sphere),1);
    mean_searchlight_accuracies_rand1=zeros(1,numRandomization);% for 100000 randomization
    flag=zeros(length(subj_sl.adj_sphere),1);
    % Storing the if correct label was predicted for each timepoint for each
    % searchlight (it will not give for each voxel, just each for each
    % timepoint in that searchlight)
    % Getting the number of relevant columns
    conditions_colNum=sum(sum_cond_mat);
    individual_results=zeros(length(subj_sl.adj_sphere),conditions_colNum);
    individual_labels=zeros(length(subj_sl.adj_sphere),conditions_colNum);
    
    % Getting non-zero indices
    % There are voxels which has 0 as its data
    run1_mat=get_mat(subj_sl,'selector','runs_xval_1');
    run2_mat=get_mat(subj_sl,'selector','runs_xval_2');
    run3_mat=get_mat(subj_sl,'selector','runs_xval_3');
    run4_mat=get_mat(subj_sl,'selector','runs_xval_4');
    % Getting the indices of the timepoints in the data mat
    run1_indices=find(run1_mat==2);
    run2_indices=find(run2_mat==2);
    run3_indices=find(run3_mat==2);
    run4_indices=find(run4_mat==2);
    % Getting the actual data to check if it is nonzero
    run1_data=data(:,run1_indices);
    run2_data=data(:,run2_indices);
    run3_data=data(:,run3_indices);
    run4_data=data(:,run4_indices);
    % Finding the nonzeros data and intersecting the common voxels.
    non01=find(sum(run1_data,2));
    non02=find(sum(run2_data,2));
    non03=find(sum(run3_data,2));
    non04=find(sum(run4_data,2));
    non0304=intersect(non03,non04);
    non0102=intersect(non01,non02);
    % Finally the indices of the voxels where the data is not zero.
    % This is only used for training and testing and I am not using it for
    % anything else. 
    NonZero_indices=intersect(non0102,non0304);
    %NonZero_indices=find(sum(data,2));
    for j=1:length(subj_sl.adj_sphere)
        %subj_sl.adj_sphere(j,:)
        inds_cols=find(subj_sl.adj_sphere(j,:)~=0);
        sphere_inds_temp=subj_sl.adj_sphere(j,inds_cols);
        %Let me remove the useless indices from here itself. Voxels which
        %have 0 throughout all the TRs.
        sphere_inds=transpose(intersect(sphere_inds_temp,NonZero_indices));
        accuracies4way_temp=zeros(4,numRandomization); % for 20 randomization
        for i=1:4 % this takes about 0.1 sec in total
            % running the correct classification
            runs_mat{i}=get_mat(subj_sl,'selector',sprintf('runs_xval_%d',i));
            training_indices{i}=find(runs_mat{i}==1);
            testing_indices{i}=find(runs_mat{i}==2);
            training_data=transpose(data(sphere_inds,training_indices{i}));
            testing_data=transpose(data(sphere_inds,testing_indices{i}));
            
            % % Need one label vector
            conditions_matrix=get_mat(subj_sl,'regressors','conditions_sh2');
            train_label_temp=conditions_matrix(:,training_indices{i});
            test_label_temp=conditions_matrix(:,testing_indices{i});
            train_labels=transpose(getLabels(train_label_temp));
            test_labels=transpose(getLabels(test_label_temp));
            
            % % run the classifier
            model=fitcnb(training_data,train_labels,'ClassNames',{'1','2','3','4'});
            %model=fitcnb(training_data,train_labels,'ClassNames',{'word','face','tool','number'});
            % Use the model on the test data to get predictions
            predicted_labels=predict(model,testing_data);
            predicted_labels=str2double(predicted_labels);
            predictions=(test_labels==predicted_labels);
            %cm = confusionchart(test_labels,predicted_labels);
            acc=sum(test_labels==predicted_labels,'all')/length(test_labels);
            accuracies4way(i,1)=acc;
            searchlight_accuracies(j,i)=acc;
            % It is 40 cause the test set has 40 time points - this is hard
            % coded
            individual_results(j,40*(i-1)+1:40*i)=predictions;
            individual_labels(j,40*(i-1)+1:40*i)=predicted_labels;
            
        end
        mean_searchlight_accuracies_rand1=mean(accuracies4way_temp,1);
        mean_sa=mean(accuracies4way);
    end
    mean_searchlight_accuracies=mean(searchlight_accuracies,2);
    timeElapsed1=toc
    
    %% Writing to afni
    subj_sl=duplicate_object(subj_sl,'pattern','epi_z','mean_searchlight_results_WB');
    subj_sl=set_mat(subj_sl,'pattern','mean_searchlight_results_WB',mean_searchlight_accuracies);
    
    %Writing it out to AFNI
    % run is the latest functional image that I read - at the top.
    args.view='+tlrc';
    write_to_afni(subj_sl,'pattern','mean_searchlight_results_WB',run,args);
    
    % Saving the workspace:
    fname=sprintf('Searchlight_classification_results_%s.mat',subj);
    save(fname);
    subj
cd ..
end


%Creating labels
function labels_temp = getLabels(labels)
%This is assuming that labels are one-hot encoded. Columns are the
%timepoints and rows are the different label rows.
    [numLabels ,numTimepoints] = size(labels);
    labels_temp=zeros(1,numTimepoints);
    for j=1:numTimepoints
        ind=find(labels(:,j));
        labels_temp(1,j)=ind;
    end
end

%% -----------------------CODE DUMP---------------------------------------
% % COMMENTED CODE WHICH IS OBSOLETE NOW.

% running the permutations - commenting for now. 
% for k=1:numRandomization % reduced it to 10k
%     rand_test_labels=test_labels(randperm(length(test_labels)));
%     acc=sum(rand_test_labels==predicted_labels,'all')/length(rand_test_labels);
%     accuracies4way_temp(i,k)=acc;
% end

% % Random model - where I have randomized the labels- this takes a lot of
% time
% rand_accuracies4way=zeros(4,10);
% % Storing the accuracies for each cross validation for each searchlight. 
% mean_searchlight_accuracies_rand2=zeros(length(subj_sl.adj_sphere),100);
% % Storing the if correct label was predicted for each timepoint for each
% % searchlight (it will not give for each voxel, just each for each
% % timepoint in that searchlight)
% %rand_individual_results=zeros(length(subj_sl.adj_sphere),160);
% % tic;
% for j=1:length(subj_sl.adj_sphere)
%     %subj_sl.adj_sphere(j,:)
%     inds_cols=find(subj_sl.adj_sphere(j,:)~=0);
%     sphere_inds=subj_sl.adj_sphere(j,inds_cols);
%     %sphere_inds
%     for i=1:4
%         for k=1:100
%             runs_mat{i}=get_mat(subj_sl,'selector',sprintf('runs_xval_%d',i));
%             training_indices{i}=find(runs_mat{i}==1);
%             testing_indices{i}=find(runs_mat{i}==2);
%             training_data=transpose(data(sphere_inds,training_indices{i}));
%             testing_data=transpose(data(sphere_inds,testing_indices{i}));
%             %size(testing_data)
%             % % Need one label vector
%             conditions_matrix=get_mat(subj_sl,'regressors','conditions_sh2');
%             %rand_cond_mat=conditions_matrix(:,randperm(length(conditions_matrix)));-cant
%             %do this since it will give us the 0 labels as well. 
%             % Just randomizing the indices which have labels
%             % cond_indices=find(sum(conditions_matrix));
%             % rand_cond_indices=cond_indices(randperm(length(cond_indices)));
%             % rand_cond_mat=conditions_matrix;
%             % rand_cond_mat(:,cond_indices)=conditions_matrix(:,rand_cond_indices);
%             train_label_temp=conditions_matrix(:,training_indices{i});
%             test_label_temp=conditions_matrix(:,testing_indices{i});
%             train_labels=transpose(getLabels(train_label_temp));
%             train_labels=train_labels(randperm(length(train_labels))); % random
%             test_labels=transpose(getLabels(test_label_temp));
%             test_labels=test_labels(randperm(length(test_labels))); % random
%             %size(train_labels)
%             %size(test_labels)
%             % put the classifier
%             model_rand=fitcnb(training_data,train_labels,'ClassNames',{'1','2','3','4'});
%             %model=fitcnb(training_data,train_labels,'ClassNames',{'word','face','tool','number'});
%             % Use the model on the test data to get predictions
%             rand_predicted_labels=predict(model_rand,testing_data);
%             rand_predicted_labels=str2double(rand_predicted_labels);
%             predictions=test_labels==rand_predicted_labels;
%             %cm = confusionchart(test_labels,rand_predicted_labels);
%             acc=sum(test_labels==rand_predicted_labels,'all')/length(test_labels);
%             rand_accuracies4way(i,k)=acc;
%             % It is 40 cause the test set has 40 time points
%             %rand_individual_results(j,40*(i-1)+1:40*i)=predictions;
%         end
%     end
%     mean_searchlight_accuracies_rand2(j,:)=mean(rand_accuracies4way,1);
%     %mean(accuracies4way)
% end
% 
% timeElapsed2=toc

 % Removing all the non-significant accuracies - commenting this as well
% for ind=1:length(flag)
%     if flag(ind)==0
%         mean_searchlight_accuracies(ind)=0;
%     end
% end

% for i=1:10
%   raw_filenames{i} = sprintf('haxby8_r%i+orig',i);
% end
%raw_fns={'data/sub-128/func/sub-128_task-matchjudgment_run-01_bold+orig','data/sub-128/func/sub-128_task-matchjudgment_run-02_bold+orig','data/sub-128/func/sub-128_task-matchjudgment_run-03_bold+orig','data/sub-128/func/sub-128_task-matchjudgment_run-04_bold+orig'};
%raw_fns={'data/s128/pb05.s128.r02.empty+orig','data/s128/pb05.s128.r04.empty+orig','data/s128/pb05.s128.r06.empty+orig','data/s128/pb05.s128.r08.empty+orig'};
%raw_fns={'pb05.s128.r02.empty+orig','pb05.s128.r04.empty+orig','pb05.s128.r06.empty+orig','pb05.s128.r08.empty+orig'};
%subj_sl = feature_select(subj_sl,'epi_z','conditions_sh2','runs_xval');
    
% class_args.train_funct_name = 'train_gnb';
% class_args.test_funct_name = 'test_gnb';

%% Pearson correlation - commenting from here. 

% R=corrcoef(transpose(individual_results));
% imagesc(R);
% 
% filtered_sa=[];
% filtered_ir=[];
% ind=1;
% % Filter accuracies that are more than 0.275 (more than chance)
% for s=1:length(searchlight_accuracies)
%     if mean(searchlight_accuracies(s,:))>=0.275
%         filtered_sa(ind,:)=searchlight_accuracies(s,:);
%         filtered_ir(ind,:)=individual_results(s,:);
%         ind=ind+1;
%     end
% end
% 
% %% Writing the searchlight results to afni
% mean_searchlight_accuracies=mean(searchlight_accuracies,2);
% subj_sl=duplicate_object(subj_sl,'pattern','epi_z','mean_sa');
% subj_sl=set_mat(subj_sl,'pattern','mean_sa',mean_searchlight_accuracies);
% write_to_afni(subj_sl,'pattern','mean_sa','pb05.s128.r02.empty+orig');
% 

% they have used anova for feature select, but I dont need that. Directly
% going to the classification. 
%[subj_sl results] = cross_validation(subj_sl,'epi_z','conditions_sh2','runs_xval','roi_mask',class_args);
%[subj_sl results_fs] = cross_validation(subj_sl,'epi_z','conditions_sh2','runs_xval','epi_z_thresh0.05',class_args);
