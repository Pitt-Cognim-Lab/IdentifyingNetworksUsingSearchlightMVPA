subject_dirs=["s103","s105","s107","s108","s109","s110","s112","s113","s114","s115","s116","s117","s118","s119","s120","s121","s122","s123","s126","s128"];

subj_mask = init_subj('Flags','ttest_mask');
% Be present in the Functional data directory.
subj_mask=load_afni_mask(subj_mask,'roi_mask','Whole_brain_mask_stan_re+tlrc');
% Hardcoded the path here
sig_acc='../Mean_searchlight_results/ttest.results/thr05_mask+tlrc';
subj_mask = load_afni_pattern(subj_mask,'flag','roi_mask',sig_acc);
flag_mat=get_mat(subj_mask,'pattern','flag');
for s=1:length(subject_dirs)
    subj=subject_dirs(s);
    subj
    cd (subj);
    fname=sprintf('Searchlight_classification_results_%s.mat',subj);
    load(fname);
    significant_individual_results=individual_results(find(flag_mat),:);
    output_name=sprintf('Significant_searchlight_acc_%s.csv',subj);
    writematrix(significant_individual_results,output_name);
    cd ..
end

    