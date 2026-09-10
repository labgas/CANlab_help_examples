# CANlab_help_examples (LaBGAS fork) — dependency overview

> **Generated file — do not edit.** Regenerate with
> `LaBGAScore_dep_report('/data/master_github_repos/CANlab_help_examples', 'files', <the 19 files documented here>, 'outdir', '/data/master_github_repos/CANlab_help_examples/Second_level_analysis_template_scripts')`
> (see `clean/LaBGAScore_dep_report.m` in LaBGAScore).
> Generated 2026-09-10 by MATLAB 2021a.

This document records which **external** functions each file calls and which
repository those live in. Calls that resolve back into this repository, and
MathWorks' own functions, are omitted.

Direct calls only (depth 1). The transitive closure of a second-level
script runs to thousands of files and reduces to "CanlabCore and SPM";
it is what the provenance tooling uses, for a different purpose.

## How to read this

Resolution is static, and MATLAB makes some of it genuinely undecidable.
Every edge carries a confidence, and the uncertain ones are reported rather
than guessed:

| Confidence | Meaning |
|---|---|
| `resolved` | exactly one definition exists for this name |
| `ambiguous` | several classes define it — e.g. `threshold` exists in `@atlas`, `@glm_map`, `@image_vector` and `@statistic_image`. Deciding which one runs needs type inference these workspace-chained scripts do not support. All candidates are listed. |
| `dotcall` | called as `obj.name(...)`, matched to a class method by name. Could also be a struct field. |
| `dynamic` | the file uses `feval`/`eval`/`str2func`, so its real call set cannot be recovered statically |
| `unparseable` | the file has a syntax error and could not be walked |

## Summary

| Script | Domain | Depends on | Direct calls | Caveats |
|---|---|---|---:|---:|
| `a2_set_default_options` | Second_level_analysis_template_scripts | — | 0 | 0 |
| `a_set_up_paths_always_run_first` | Second_level_analysis_template_scripts | LaBGAScore | 2 | 0 |
| `c2_SVM_contrasts_masked` | Second_level_analysis_template_scripts | CanlabCore, LaBGAScore | 11 | 1 |
| `c2a_second_level_regression` | Second_level_analysis_template_scripts | CanlabCore, LaBGAScore | 17 | 2 |
| `c2f_run_MVPA_regression_single_trial` | Second_level_analysis_template_scripts | CanlabCore, ooFmriDataObjML | 22 | 18 |
| `c2g_run_multivariate_mediation_single_trial` | Second_level_analysis_template_scripts | CanlabCore, MediationToolbox | 11 | 4 |
| `d10_signature_riverplots` | Second_level_analysis_template_scripts | CanlabCore | 3 | 0 |
| `d_signature_responses_generic` | Second_level_analysis_template_scripts | — | 0 | 0 |
| `e1_corr_patterns` | Second_level_analysis_template_scripts | CanlabCore | 10 | 0 |
| `h_signature_responses_group_diff` | Second_level_analysis_template_scripts | CanlabCore | 5 | 0 |
| `prep_1_set_conditions_contrasts_colors` | Second_level_analysis_template_scripts | CanlabCore | 2 | 0 |
| `prep_1b_prep_behavioral_data` | Second_level_analysis_template_scripts | — | 0 | 0 |
| `prep_2_load_image_data_and_save` | Second_level_analysis_template_scripts | CanlabCore, canlab_single_trials, spm12 | 9 | 1 |
| `prep_3_calc_univariate_contrast_maps_and_save` | Second_level_analysis_template_scripts | CanlabCore | 4 | 0 |
| `prep_3a_run_second_level_regression_and_save` | Second_level_analysis_template_scripts | CanlabCore, LaBGAScore | 31 | 4 |
| `prep_3c_run_SVMs_on_contrasts_masked` | Second_level_analysis_template_scripts | CanlabCore, LaBGAScore, ooFmriDataObjML | 22 | 5 |
| `prep_3f_create_fmri_data_single_trial_object` | Second_level_analysis_template_scripts | CanlabCore, canlab_single_trials | 3 | 1 |
| `prep_3g_create_fmri_data_runwise_contrast_object` | Second_level_analysis_template_scripts | CanlabCore, canlab_single_trials | 2 | 1 |
| `prep_4_apply_signatures_and_save` | Second_level_analysis_template_scripts | CanlabCore, MasksPrivate, Neuroimaging_Pattern_Masks | 3 | 0 |

## Dependencies by repository

| Repository | Call edges | Distinct functions |
|---|---:|---:|
| CanlabCore | 159 | 59 |
| ooFmriDataObjML | 17 | 12 |
| LaBGAScore | 11 | 8 |
| canlab_single_trials | 4 | 2 |
| spm12 | 3 | 2 |
| MasksPrivate | 1 | 1 |
| MediationToolbox | 1 | 1 |
| Neuroimaging_Pattern_Masks | 1 | 1 |

## Per-script detail

### `a2_set_default_options`

`Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/a2_set_default_options.m`

No external dependencies.

### `a_set_up_paths_always_run_first`

`Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/a_set_up_paths_always_run_first.m`

**LaBGAScore**

- `LaBGAScore_firstlevel_s1_options_dsgn_struct`
- `LaBGAScore_prep_s0_define_directories`

### `c2_SVM_contrasts_masked`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2_SVM_contrasts_masked.m`

**CanlabCore**

- `barplot_columns`
- `colormap_tor`
- `create_figure`
- `load_atlas`
- `montage` *(@region)* — `ambiguous`
- `print_matrix`
- `region` *(@region)*
- `roc_plot`
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*

**LaBGAScore**

- `thresholded_fmri_data_from_statistic_image`

### `c2a_second_level_regression`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2a_second_level_regression.m`

**CanlabCore**

- `addblobs` *(@fmridisplay)*
- `apply_mask` *(@image_vector)*
- `atlas2region` *(@atlas)*
- `canlab_results_fmridisplay`
- `downsample_parcellation` *(@atlas)* — `dotcall`
- `fmri_mask_image` *(@fmri_mask_image)*
- `load_atlas`
- `montage` *(@region)* — `ambiguous`
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `table_of_atlas_regions_covered` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*

**LaBGAScore**

- `LaBGAScore_region_table`
- `LaBGAScore_region_table_safe`
- `LaBGAScore_smart_parallel_pool_setup`
- `tfce_fwe_from_null`

### `c2f_run_MVPA_regression_single_trial`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2f_run_MVPA_regression_single_trial.m`

**CanlabCore**

- `addblobs` *(@fmridisplay)*
- `apply_mask` *(@image_vector)* — `dotcall`
- `canlab_results_fmridisplay`
- `fit` *(@glm_map)* — `dotcall`, 2 candidates
- `fmri_mask_image` *(@fmri_mask_image)*
- `line_plot_multisubject`
- `pipeline` *(@pipeline)* — `ambiguous`
- `plot` *(@fmri_data)* — `dotcall`, 11 candidates
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `rescale` *(@fmri_data)* — `dotcall`
- `test` *(@algorithm)* — `dotcall`
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*

**ooFmriDataObjML**

- `bayesOptCV`
- `crossValScore`
- `cvpartition2`
- `fmri2VxlFeatTransformer`
- `get_mse`
- `mlpcrRegressor`
- `pcrRegressor`
- `pipeline` — `ambiguous`
- `plsRegressor`

### `c2g_run_multivariate_mediation_single_trial`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2g_run_multivariate_mediation_single_trial.m`

**CanlabCore**

- `apply_mask` *(@image_vector)* — `dotcall`
- `autolabel_regions_using_atlas` *(@region)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `history` *(@image_vector)* — `dotcall`
- `load_atlas`
- `montage` *(@region)* — `ambiguous`
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `rescale` *(@fmri_data)* — `dotcall`
- `title_montage` *(@fmridisplay)*

**MediationToolbox**

- `multivariateMediation`

### `d10_signature_riverplots`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/d10_signature_riverplots.m`

**CanlabCore**

- `load_image_set`
- `riverplot` *(@fmri_data)* — `ambiguous_within_repo`, 2 candidates
- `seaborn_colors`

### `d_signature_responses_generic`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/d_signature_responses_generic.m`

No external dependencies.

### `e1_corr_patterns`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/e1_corr_patterns.m`

**CanlabCore**

- `addblobs` *(@fmridisplay)*
- `apply_mask` *(@image_vector)*
- `canlab_results_fmridisplay`
- `fmri_mask_image` *(@fmri_mask_image)*
- `load_atlas`
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `searchlight_correlation`
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*

### `h_signature_responses_group_diff`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/h_signature_responses_group_diff.m`

**CanlabCore**

- `barplot_columns`
- `create_figure`
- `mediansplit`
- `seaborn_colors`
- `ttest2_printout`

### `prep_1_set_conditions_contrasts_colors`

`Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/prep_1_set_conditions_contrasts_colors.m`

**CanlabCore**

- `colorcube_colors`
- `format_strings_for_legend`

### `prep_1b_prep_behavioral_data`

`Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/prep_1b_prep_behavioral_data.m`

No external dependencies.

### `prep_2_load_image_data_and_save`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_2_load_image_data_and_save.m`

**CanlabCore**

- `check_valid_imagename`
- `create_figure`
- `enforce_variable_types` *(@image_vector)*
- `fmri_data` *(@fmri_data)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `qc_metrics_second_level` *(@image_vector)*

**canlab_single_trials**

- `fmri_data_st` *(@fmri_data_st)*

**spm12**

- `conditions` *(@meeg)* — `dotcall`
- `spm_select` — `ambiguous_within_repo`, 2 candidates

### `prep_3_calc_univariate_contrast_maps_and_save`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3_calc_univariate_contrast_maps_and_save.m`

**CanlabCore**

- `create_figure`
- `enforce_variable_types` *(@image_vector)*
- `qc_metrics_second_level` *(@image_vector)*
- `replace_empty` *(@image_vector)*

### `prep_3a_run_second_level_regression_and_save`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3a_run_second_level_regression_and_save.m`

**CanlabCore**

- `addblobs` *(@fmridisplay)*
- `apply_mask` *(@image_vector)*
- `apply_parcellation` *(@image_vector)*
- `atlas2region` *(@atlas)*
- `barplot_columns`
- `canlab_results_fmridisplay`
- `create_figure`
- `downsample_parcellation` *(@atlas)* — `dotcall`
- `estimateBayesFactor` *(@statistic_image)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `get_wh_image` *(@image_vector)* — `dotcall`
- `getvif`
- `hansen_neurotransmitter_maps` *(@image_vector)*
- `load_atlas`
- `merge_atlases` *(@atlas)*
- `montage` *(@region)* — `ambiguous`
- `plot_correlation_matrix`
- `plot_vertical_line`
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `robfit_parcelwise` *(@fmri_data)*
- `run_diagnostics` *(@glm_map)*
- `scn_standard_colors`
- `ste`
- `test` *(@algorithm)* — `dotcall`
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*
- `validate_object` *(@fmri_data)* — `ambiguous_within_repo`, 2 candidates

**LaBGAScore**

- `LaBGAScore_smart_parallel_pool_setup`
- `group_tfce_from_subject_maps`
- `thresholded_fmri_data_from_statistic_image`

### `prep_3c_run_SVMs_on_contrasts_masked`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3c_run_SVMs_on_contrasts_masked.m`

**CanlabCore**

- `apply_mask` *(@image_vector)*
- `enforce_variable_types` *(@image_vector)*
- `fit` *(@glm_map)* — `dotcall`, 2 candidates
- `fmri_mask_image` *(@fmri_mask_image)*
- `group` *(@group)*
- `montage` *(@region)* — `ambiguous`
- `pipeline` *(@pipeline)* — `ambiguous`
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `searchlight_disti_Lukas`
- `sec2hms`
- `stability_selection` *(@predictive_model)*
- `title_montage` *(@fmridisplay)*
- `trim_mask` *(@image_vector)*

**LaBGAScore**

- `LaBGAScore_smart_parallel_pool_setup`

**ooFmriDataObjML**

- `bayesOptCV`
- `crossValScore`
- `cvpartition2`
- `fmri2VxlFeatTransformer`
- `get_f1_macro`
- `get_hinge_loss`
- `linearSvmClf`
- `pipeline` — `ambiguous`

### `prep_3f_create_fmri_data_single_trial_object`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3f_create_fmri_data_single_trial_object.m`

**CanlabCore**

- `get_wh_image` *(@image_vector)* — `dotcall`
- `remove_empty` *(@image_vector)*

**canlab_single_trials**

- `fmri_data_st` *(@fmri_data_st)*
- `get_wh_image` *(@fmri_data_st)*

### `prep_3g_create_fmri_data_runwise_contrast_object`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3g_create_fmri_data_runwise_contrast_object.m`

**CanlabCore**

- `image_math` *(@image_vector)*

**canlab_single_trials**

- `fmri_data_st` *(@fmri_data_st)*

### `prep_4_apply_signatures_and_save`

`Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_4_apply_signatures_and_save.m`

**CanlabCore**

- `ste`

**MasksPrivate**

- `apply_nps`

**Neuroimaging_Pattern_Masks**

- `apply_all_signatures`

