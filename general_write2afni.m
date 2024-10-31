%% Writing to afni 
% flag contains the indices where the searchlight is significant.
% Getting the mask for flags;
subj_mask = init_subj('Flags','ttest_mask');
% subj_mask=load_afni_mask(subj_mask,'roi_mask','../subjects_data/Functional_data/Whole_brain_mask_stan_re+tlrc');
subj_mask=load_afni_mask(subj_mask,'roi_mask','Whole_brain_mask_stan_re+tlrc');
% Hardcoded the path here
sig_acc='../subjects_data/Mean_searchlight_results/ttest.results/thr05_mask+tlrc';
subj_mask = load_afni_pattern(subj_mask,'flag','roi_mask',sig_acc);
flag_mat=get_mat(subj_mask,'pattern','flag');

% Read the corresponding matrix here
comps=readmatrix("../Components_map_7c_2std.csv");
% comps=readmatrix("../Components_map_concatenated2_7c_2std.csv");

num_components=size(comps,2);
%spatial_maps=zeros(length(flag_mat),num_components);
spatial_map=zeros(length(flag_mat),num_components);
ind=1;
for i=1:length(flag_mat)
        if flag_mat(i)==1
            spatial_map(i,:)=comps(ind,:);
            ind=ind+1;
        end
end

%% Writing the maps out to afni
% Inidividual components
for col=1:size(spatial_map,2)
    obj_name=append('Component_maps_7c_2std_',int2str(col)); %% Change name here
    subj_mask=duplicate_object(subj_mask,'pattern', 'flag',obj_name);
    subj_mask=set_mat(subj_mask,'pattern',obj_name, spatial_map(:,col));
    args.view='+tlrc';
    write_to_afni(subj_mask,'pattern',obj_name,'pb05.s128.r02.empty_re_stan+tlrc',args);
end
% one single components
% obj_name='Component_map__voxelValue_9c';
% subj_mask=duplicate_object(subj_mask,'pattern', 'flag',obj_name);
% subj_mask=set_mat(subj_mask,'pattern',obj_name, spatial_map);
% args.view='+tlrc';
% write_to_afni(subj_mask,'pattern',obj_name,'../pb05.s128.r02.empty_re_stan+tlrc',args);

