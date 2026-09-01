# Second-Level Analysis Template Scripts — LaBGAS Usage Guide

This document describes how **LaBGAS** (Laboratory for Brain-Gut Axis Studies, KU Leuven) uses the second-level (group) fMRI analysis scripts in this directory. It covers only the scripts LaBGAS actively maintains and runs — not the full ~80-script generic CANlab toolbox this directory also contains. For that broader background, see `a0_begin_here_readme.m` and `list_of_scripts_and_workflow.m`.

If you are new to this repo, read this file top to bottom once, then use the [script reference table](#script-reference) as a lookup when running or editing a specific script.

## Contents

- [What this is](#what-this-is)
- [Directory structure and naming convention](#directory-structure-and-naming-convention)
- [The two script groups](#the-two-script-groups)
- [Per-study workflow](#per-study-workflow)
- [Shared data model](#shared-data-model)
- [`a2_set_default_options.m` and the capitalization convention](#a2_set_default_optionsm-and-the-capitalization-convention)
- [Script reference](#script-reference)
- [The `LaBGAScore` dependency](#the-labgascore-dependency)
- [Out of scope](#out-of-scope)

## What this is

`CANlab_help_examples/Second_level_analysis_template_scripts` is a template-script framework (not a software package — no build system, linter, or test suite) for running group-level analyses on first-level beta/contrast images: univariate GLMs, cross-validated SVMs, CANlab "signature" pattern responses (NPS, SIIPS1, etc.), searchlight correlations, and single-trial/runwise MVPA and mediation analyses. This working copy is a **LaBGAS lab fork** of the generic CANlab template. LaBGAS uses a curated, actively-maintained subset of the scripts (listed below); the rest of the toolbox is generic CANlab machinery that LaBGAS does not currently use in this form.

Every in-scope script follows one shape: a MATLAB comment header (USAGE/OPTIONS/NOTES), a call to `a_set_up_paths_always_run_first` (which sets up paths and options), optional custom-option overrides, a guarded reload of any `.mat` files it needs, the actual computation, then save/publish. Scripts are meant to be run as MATLAB `.m` scripts (cell-by-cell or with `publish()`), not called as functions — and, per the convention below, never run directly from this repo.

## Directory structure and naming convention

Every LaBGAS project is organized as a **DataLad superdataset** `proj_xxx`, with subdatasets `sourcedata`, `BIDS`, `derivatives`, `code`, `firstlevel`, and `secondlevel`. This framework's scripts mainly touch `firstlevel`, `secondlevel`, `code`, and `BIDS` (for phenotype data).

For any given second-level model, there is a matching pair of folders: `code/secondlevel/model_x/` (the scripts) and `secondlevel/model_x/` (the scripts' outputs — results, masks, figures, published HTML). `a_set_up_paths_always_run_first.m` builds the `secondlevel/model_x/` side of this (creating `masks/`, `results/`, `results/figures/`, `results/notes/`, `results/html/` if missing).

**None of the scripts documented here are ever run directly from this GitHub repo.** Every script — both the 4 Group 1 scripts and all Group 2 scripts — is copied into the study's `code/secondlevel/model_x/` folder and renamed:

```
projname_mM_sN_scriptname.m
```

where `M` is the model number and `N` is the script's sequential order within that model — e.g. `myproj_m1_s5_prep_3a_run_second_level_regression_and_save.m`. The `core_scripts_to_run_without_modifying/` directory name refers to not changing the *analysis logic* of Group 2 scripts, not to running them in place from the repo.

- **Group 1 copies** are substantively edited per study (paths, conditions, contrasts, behavioral data).
- **Group 2 copies** are only *lightly* edited — typically just the first section, updating the call to `a_set_up_paths_always_run_first` to point at the study's own renamed copy of that script.
- **Versioning within a model:** if the same analysis needs to run multiple times with different options (e.g. voxelwise vs. parcelwise second-level regression via `prep_3a`/`c2a`), a single `a2_set_default_options.m` is maintained per model, holding the defaults for the primary variant. Additional variants are separate copies of just the varying script, with a letter appended to the **sequence number** rather than the script's own name — e.g. `..._s5_prep_3a...` (variant 1), `..._s5a_prep_3...` (variant 2), `..._s5b_prep_3a...` (variant 3). Each variant overrides only the options that differ from `a2`'s defaults, in its own "SET CUSTOM OPTIONS" section near the top of the script.

## The two script groups

**Group 1 — user-edited entry points**, in `b_copy_to_local_scripts_dir_and_modify/`:

- `a_set_up_paths_always_run_first.m`
- `a2_set_default_options.m`
- `prep_1_set_conditions_contrasts_colors.m`
- `prep_1b_prep_behavioral_data.m`

**Group 2 — core analysis scripts**, in `core_scripts_to_run_without_modifying/` (15 scripts, listed in full in the [script reference](#script-reference) below).

Both directories also contain many other scripts not covered here — see [Out of scope](#out-of-scope).

Not every script gets copied into every study's model folder: the 4 Group 1 scripts plus `prep_2_load_image_data_and_save.m` and `prep_3_calc_univariate_contrast_maps_and_save.m` are always needed (core setup and image/contrast loading, a prerequisite for everything else). The remaining 12 Group 2 scripts are added only as a given model's specific analyses require them.

## Per-study workflow

1. Copy the always-needed scripts (4 Group 1 + `prep_2` + `prep_3`) into `code/secondlevel/model_x/`, renamed per the convention above.
2. Edit the Group 1 copies substantively: `a_set_up_paths_always_run_first.m` (paths, calls `LaBGAScore_prep_s0_define_directories` + `LaBGAScore_firstlevel_s1_options_dsgn_struct`) and `a2_set_default_options.m` (all per-script defaults, in one place, capitalized-header convention explained [below](#a2_set_default_optionsm-and-the-capitalization-convention)).
3. Run `prep_1_set_conditions_contrasts_colors.m` (defines `DAT`), optionally `prep_1b_prep_behavioral_data.m` (attaches `DAT.BEHAVIOR`/`DAT.BETWEENPERSON` — required for group definitions). Unlike `a_set_up_paths_always_run_first.m` and `a2_set_default_options.m`, these two scripts are where the study's conditions, contrasts, colors, behavioral variables, and between-person groups are actually defined, so they typically need **extensive**, not light, study-specific editing — the checked-in versions are worked examples for one study's design, not templates to run with minor tweaks. Many more worked examples are available in LaBGAS's project folders on the KU Leuven server and in LaBGAS's GIN/GitHub repos.
4. Lightly edit the `prep_2`/`prep_3` copies' first section (the `a_set_up_paths_always_run_first` call).
5. Run `prep_2_load_image_data_and_save.m` → `prep_3_calc_univariate_contrast_maps_and_save.m`.
6. Copy in and run whichever further Group 2 scripts the model needs, lightly edited as above, in roughly this dependency order:
   - **Univariate GLM:** `prep_3a_run_second_level_regression_and_save.m` → `c2a_second_level_regression.m`
   - **SVM:** `prep_3c_run_SVMs_on_contrasts_masked.m` → `c2_SVM_contrasts_masked.m`
   - **Signatures:** `prep_4_apply_signatures_and_save.m` → `d_signature_responses_generic.m` / `d10_signature_riverplots.m` / `h_signature_responses_group_diff.m`
   - **Searchlight correlation:** `e1_corr_patterns.m` (no `prep_4` dependency)
   - **Single-trial / runwise MVPA & mediation:** `prep_3f_create_fmri_data_single_trial_object.m` or `prep_3g_create_fmri_data_runwise_contrast_object.m` → `c2f_run_MVPA_regression_single_trial.m` / `c2g_run_multivariate_mediation_single_trial.m`. These additionally require single-trial con images already produced by `LaBGAScore_firstlevel_s2_fit_model.m` (first-level, not part of this repo).
7. For a script that must run multiple times within the model with different options, add a lettered sequence-number copy (`..._s5a_...`, `..._s5b_...`) rather than duplicating `a2_set_default_options.m`; override only the differing options in that copy's own "SET CUSTOM OPTIONS" section.

You do **not** need to explicitly reload `.mat` files before Group 2 scripts. From `prep_3a` onward (and in every lettered on-demand script), each script contains its own guarded reload — `if ~exist('DAT','var') ... load(...)`, etc. — so it transparently reloads `image_names_and_setup.mat`, `data_objects.mat`/`data_objects_scaled.mat`, and `contrast_data_objects.mat` as needed. The generic template's `b_reload_saved_matfiles.m` is only relevant earlier/manually (e.g., if you want a fresh MATLAB session to have `DAT` in the workspace before editing `prep_1b` interactively).

## Shared data model

- **`DAT`** — one struct threaded through the whole pipeline: `conditions`, `contrasts`, `contrastnames`, `colors`/`contrastcolors`, `subfolders`/`functional_wildcard` (set in `prep_1`), `BEHAVIOR`/`BETWEENPERSON` (set in `prep_1b`), `SIG_conditions`/`SIG_contrasts`/`NPSsubregions`/`npsresponse`/`npscontrasts` (set in `prep_4`).
- **`DATA_OBJ`** / **`DATA_OBJsc`** — cell arrays (one cell per condition) of CanlabCore `fmri_data`/`fmri_data_st` objects, raw and z-scored respectively.
- **`DATA_OBJ_CON`** / **`DATA_OBJ_CONsc`** / **`DATA_OBJ_CONscc`** — same, per contrast: raw, z-scored-before-contrast, and l2norm-scaled-after-contrast.
- Persisted across scripts as `image_names_and_setup.mat` (DAT, DSGN, directory names, helper function handles `printhdr`/`printstr`), `data_objects.mat` + `data_objects_scaled.mat` (`DATA_OBJ`/`DATA_OBJsc`), and `contrast_data_objects.mat` (`DATA_OBJ_CON*`), all under `resultsdir`.

Because `secondlevel/model_x/results/` is DataLad/git-annex managed, some scripts explicitly `git annex unannex` a `.mat` file before overwriting it (e.g. `prep_1b_prep_behavioral_data.m`, `prep_4_apply_signatures_and_save.m`) to avoid write-permission errors on already-annexed files.

## `a2_set_default_options.m` and the capitalization convention

`a2_set_default_options.m` centralizes every option consumed by Group 2 scripts, one `%%` section per script. Per the file's own convention comment: *"If the title of the section below is capitalized, the scripts and their options have been revamped by @lukasvo76 already."* Capitalized sections mark scripts LaBGAS actively maintains — this is exactly how the Group 2 list above was derived. Lowercase sections (`prep_3d_run_SVMs_betweenperson_contrasts options`, `z_batch_publish_everything, z_batch_publish_analyses options`) are explicitly **not** revamped and out of scope here.

One exception worth knowing: `h_signature_responses_group_diff.m` has no dedicated section of its own — it reuses `keyword_sigs`, `myscaling_sigs`, and `similarity_metric_sigs` straight from the `PREP_4_APPLY_SIGNATURES_AND_SAVE` section.

## Script reference

### Group 1

| Script                                       | Role                                                                                                                                                                                                                                                                                                                                                                                               |
| -------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `a_set_up_paths_always_run_first.m`        | Must run first (directly or via another script calling it). Verifies/creates`rootdir`/`DSGN` via LaBGAScore functions, calls `a2_set_default_options`, ensures CanlabCore/CanlabPrivate/CANlab_help_examples/canlab_single_trials are cloned and on path, builds the standard `secondlevel/model_x` subdirectory tree, defines the `printhdr`/`printstr` helpers used everywhere else. |
| `a2_set_default_options.m`                 | All default options for every Group 2 script, in one file, organized by capitalized`%%` section per script (see above). Auto-called by `a_set_up_paths_always_run_first.m`; rarely run standalone.                                                                                                                                                                                             |
| `prep_1_set_conditions_contrasts_colors.m` | Defines`DAT.conditions` (with file-location wildcards for `prep_2` to use), `DAT.contrasts`/`contrastnames` (within-person contrast matrix), and `DAT.colors`/`contrastcolors`. Saves `image_names_and_setup.mat`. Study-specific by design — see [workflow step 3](#per-study-workflow) — the checked-in version is one worked example, not a generic template.                            |
| `prep_1b_prep_behavioral_data.m`           | Optional. Reads BIDS`phenotype/*.tsv` files, builds/z-scores behavioral variables, stores raw tables in `DAT.BEHAVIOR`, and builds `DAT.BETWEENPERSON.group` (single group vector, 1/-1 coded) and/or per-condition/per-contrast covariate tables in `DAT.BETWEENPERSON.conditions`/`.contrasts`. Study-specific by design — see [workflow step 3](#per-study-workflow).                   |

### Group 2

Always copied in alongside Group 1: `prep_2_load_image_data_and_save.m`, `prep_3_calc_univariate_contrast_maps_and_save.m`. The remaining 12 are added per model as needed.

| Script                                                 | Category         | Role                                                                                                                                                                                                                                                                                                                                                       | Key options (in`a2_set_default_options.m`)                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------------------------------------------------ | ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `prep_2_load_image_data_and_save.m`                  | `prep_`        | Loads first-level beta/con images per`DAT.conditions`/wildcards into `fmri_data_st` objects, runs QC and z-scoring, saves `data_objects*.mat`, publishes an HTML report.                                                                                                                                                                             | `dofullplot`, `omit_histograms`, `dozipimages`, `maskname_brain`, `subjs2exclude_data`                                                                                                                                                                                                                                                                                                                                                  |
| `prep_3_calc_univariate_contrast_maps_and_save.m`    | `prep_`        | Computes`DATA_OBJ_CON*` from `prep_2`'s condition objects per `DAT.contrasts`, l2norm-rescales, QC, saves `contrast_data_objects.mat`.                                                                                                                                                                                                             | (shares the`prep_2` section above)                                                                                                                                                                                                                                                                                                                                                                                                              |
| `prep_3a_run_second_level_regression_and_save.m`     | `prep_`        | Group-level regression per condition/contrast: voxelwise (`regress()`, optionally robust/TFCE — see the TFCE note below) or parcelwise (`robfit_parcelwise()`, atlas-defined parcels). Optional Bayes Factor conversion, ROI-average extraction, neurotransmitter-map similarity, and MVPA regression of covariates from between-subject brain data.                             | `maskname_glm`, `atlasname_glm`/`atlas_granularity`, `myscaling_glm`, `design_matrix_type`, `dorobust`, `dorobfit_parcelwise` (+ `csf_wm_covs`, `remove_outliers`), `doBayes`, `doTFCE` (+ perm/sidedness/tail), `doroi_analysis` (+ `roi_names`/`roi_modelname`/`roi_set_name`), `doneurotransmitter_maps`, `domvpa_reg_cov` (+ algorithm/holdout/folds)                                                       |
| `c2a_second_level_regression.m`                      | lettered (`c`) | Displays/thresholds`prep_3a` results (FDR-q, uncorrected-p, extent, Bayes Factor thresholds); optional bootstrapping of the MVPA-regression-on-covariates results.                                                                                                                                                                                       | `save_figures_glm`, `q_threshold_glm`, `p_threshold_glm`, `k_threshold_glm`, `BF_threshold_glm`, `dobootstrap_mvpa_reg_cov` (+ boot_n/parallel/cons2boot), `q_threshold_mvpa_reg_cov`, `k_threshold_mvpa_reg_cov`                                                                                                                                                                                                                 |
| `prep_3c_run_SVMs_on_contrasts_masked.m`             | `prep_`        | Cross-validated SVM per contrast (masked), via either`ooFmriDataObjML` or CANlab `predict()`; optional bootstrapping, stability selection, and searchlight SVM. Saves results.                                                                                                                                                                         | `ml_method_svm`, `holdout_set_method_svm`/`holdout_set_type_svm`/`nfolds_svm`, `maskname_svm`, `myscaling_svm`, `dosavesvmstats`, `dobootstrap_svm` (+ boot_n/cons2boot), `dostabilityselection_svm` (+ boot_n_ss/cons2ss/k_ss/threshold_ss), `dosearchlight_svm` (+ radius/cons2searchlight)                                                                                                                                 |
| `c2_SVM_contrasts_masked.m`                          | lettered (`c`) | Displays/thresholds SVM results from`prep_3c_`, with atlas-based region labeling; uses `LaBGAScore_atlas_binary_mask_from_atlas.m`-generated masks when a custom atlas is specified.                                                                                                                                                                   | `save_figures_svm`, `q_threshold_svm`, `p_threshold_svm`, `k_threshold_svm`, `atlasname_svm`                                                                                                                                                                                                                                                                                                                                            |
| `prep_3f_create_fmri_data_single_trial_object.m`     | `prep_`        | Builds a single-trial`fmri_data_st` object from single-trial con images produced by `LaBGAScore_firstlevel_s2_fit_model.m`, attaching per-trial ratings and VIF-based outlier flags.                                                                                                                                                                   | `cons2exclude_dat_st`, `behav_outcome_dat_st`, `subj_identifier_dat_st`, `cond_identifier_dat_st`, `group_identifier_dat_st` (optional), `vif_threshold_dat_st`                                                                                                                                                                                                                                                                       |
| `prep_3g_create_fmri_data_runwise_contrast_object.m` | `prep_`        | Builds an`fmri_data_st` object of runwise contrasts from condition betas, attaching runwise phenotype metadata read from a CSV in the BIDS subdataset.                                                                                                                                                                                                   | `phenofile_dat_rw`, `cons2include_dat_rw`, `behav_outcome_dat_rw`, `subj_identifier_dat_rw`, `run_included_dat_rw`, `group_identifier_dat_rw` (optional)                                                                                                                                                                                                                                                                              |
| `c2f_run_MVPA_regression_single_trial.m`             | lettered (`c`) | MVPA regression (default cross-validated PCR) predicting a continuous outcome from the`prep_3f_` single-trial object; optional bootstrapping, permutation testing, and source reconstruction ("structure coefficients").                                                                                                                                 | `ml_method_mvpa_reg_st`, `algorithm_mvpa_reg_st`, `holdout_set_method_mvpa_reg_st`, `nfolds_mvpa_reg_st`, `zscore_outcome_mvpa_reg_st`, `maskname_mvpa_reg_st`, `myscaling_mvpa_reg_st`, `dobootstrap_mvpa_reg_st` (+ boot_n/parallel), `doperm_mvpa_reg_st` (+ perm_n/sidedness), `dosourcerecon_mvpa_reg_st` (+ perm variant), `q_threshold_mvpa_reg_st`, `k_threshold_mvpa_reg_st`, `domultilevel_mvpa_reg_st` (WIP) |
| `c2g_run_multivariate_mediation_single_trial.m`      | lettered (`c`) | Multivariate (PDM) mediation analysis on a continuous outcome, on the`prep_3f_` single-trial object.                                                                                                                                                                                                                                                     | `save_figures_pdm`, `zscore_outcome_pdm`, `maskname_pdm`, `myscaling_pdm`, `nPDM`, `dobootstrap_pdm` (+ boot_n/k_threshold), `dosourcerecon_pdm`, `dosavepdmstats`                                                                                                                                                                                                                                                                |
| `prep_4_apply_signatures_and_save.m`                 | `prep_`        | Applies selected CANlab signature patterns (e.g. NPS, SIIPS1) to conditions and contrasts via`apply_all_signatures`, saving results into `DAT.SIG_conditions`/`DAT.SIG_contrasts`; computes NPS subregion responses when `nps`/`'all'` is among the selected signatures. Appends to `image_names_and_setup.mat` (`git annex unannex` first). | `myscaling_sigs`, `similarity_metric_sigs`, `keyword_sigs`                                                                                                                                                                                                                                                                                                                                                                                  |
| `d_signature_responses_generic.m`                    | lettered (`d`) | Plots and tests significance of signature responses from`prep_4_`, for individual signatures or groups, via `plugin_signature_condition_contrast_plot`.                                                                                                                                                                                                | `signatures_to_plot` (shared with `d10`)                                                                                                                                                                                                                                                                                                                                                                                                      |
| `d10_signature_riverplots.m`                         | lettered (`d`) | Cosine-similarity riverplots of signature responses vs. conditions/contrasts. Unlike`d_signature_responses_generic.m`, only works on signature *groups* as defined by `load_image_set`, not individual signatures.                                                                                                                                   | `signatures_to_plot` (shared with `d`)                                                                                                                                                                                                                                                                                                                                                                                                        |
| `h_signature_responses_group_diff.m`                 | lettered (`h`) | Two-sample t-test per contrast on signature responses (one figure per contrast, loops over signatures), plus NPS-subregion group differences. Group membership from`DAT.BETWEENPERSON.group` or the condition/contrast-specific fields set in `prep_1b_prep_behavioral_data.m`. Requires `prep_1b` to have been run with a real group variable.      | Reuses`keyword_sigs`/`myscaling_sigs`/`similarity_metric_sigs` from the `prep_4` section — no dedicated section (see [above](#a2_set_default_optionsm-and-the-capitalization-convention))                                                                                                                                                                                                                                                 |
| `e1_corr_patterns.m`                                 | lettered (`e`) | Pairwise searchlight correlation maps between all condition or contrast images (`searchlight_correlation()`), optionally masked/restricted to an atlas. Independent of `prep_4`/signatures.                                                                                                                                                            | `r_threshold_corr`, `corr_type`                                                                                                                                                                                                                                                                                                                                                                                                               |

## The `LaBGAScore` dependency

Several Group 1/2 scripts call into a separate sibling repo, **LaBGAScore** (local path `/data/master_github_repos/LaBGAScore` in this environment; upstream `github.com/labgas/LaBGAScore`), for study-specific setup that lives outside this repo:

| Function/script                                  | Called from                                                                                           | Purpose                                                                                                                                                          |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `LaBGAScore_prep_s0_define_directories`        | `a_set_up_paths_always_run_first.m`                                                                 | Defines`rootdir`, `githubrootdir`, `codedir`, `BIDSdir`, `spmrootdir`, etc. for the current study.                                                     |
| `LaBGAScore_firstlevel_s1_options_dsgn_struct` | `a_set_up_paths_always_run_first.m`                                                                 | Builds the first-level`DSGN` struct (model/condition definitions) this framework reads `DSGN.modeldir`/`DSGN.conditions`/`DSGN.contrastnames` from.      |
| `LaBGAScore_firstlevel_s2_fit_model.m`         | (upstream of`prep_3f_...`, not called by it)                                                        | Fits first-level models and produces the single-trial con images`prep_3f_create_fmri_data_single_trial_object.m` consumes.                                     |
| `LaBGAScore_atlas_binary_mask_from_atlas.m`    | referenced by`atlasname_glm`/`atlasname_svm` options in `a2_set_default_options.m`              | Generates custom`.mat` atlas/mask objects (`combined_atlas` variable) usable as `atlasname_glm`/`atlasname_svm`/`e1_corr_patterns.m`'s masking option. |
| `LaBGAScore_atlas_rois_from_atlas.m`           | referenced by`roi_names`/`roi_modelname`/`roi_set_name` options in `a2_set_default_options.m` | Generates per-ROI atlas objects for`prep_3a_...`'s `doroi_analysis` option.                                                                                  |

Three further LaBGAScore functions are called from these scripts and were missing from the
table above until the dependency tooling found them:

| Function                                        | Called from                                                                                          | Purpose                                                     |
| ----------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| `LaBGAScore_smart_parallel_pool_setup.m`      | `c2a_second_level_regression.m`, `prep_3a_...`, `prep_3c_run_SVMs_on_contrasts_masked.m`       | Sets up the parallel pool before bootstrapping/permutation. |
| `group_tfce_from_subject_maps.m`              | `prep_3a_run_second_level_regression_and_save.m`                                                   | Group TFCE from subject-level maps.                         |
| `thresholded_fmri_data_from_statistic_image.m` | `prep_3a_run_second_level_regression_and_save.m`, `c2_SVM_contrasts_masked.m`                    | Thresholded `fmri_data` object from a `statistic_image`.  |

This is intentionally a high-level summary — LaBGAScore's own internals are out of scope here.

## Dependency and provenance documentation

[`DEPENDENCIES.md`](DEPENDENCIES.md) documents what each script calls and which repository
each of those lives in. It covers exactly the **19 scripts listed above** (4 Group 1 + 15
Group 2) — the set LaBGAS actively uses and maintains — not the ~113 scripts in this
folder, the rest of which are generic CANlab machinery LaBGAS does not document.

That file, along with `dependencies.tsv` and `dependencies.yml`, is **generated** by
`LaBGAScore_dep_report` (in LaBGAScore's `clean/` folder) — regenerate rather than edit.

Because these templates are copied and renamed per study, the version of CanlabCore they
ran against is not recorded anywhere by default. LaBGAScore's `clean/LaBGAScore_prov_*`
tooling closes that gap:

- **Going forward** — publish with `LaBGAScore_prov_publish` instead of `publish`. The
  report gains a Provenance section naming the commit of every dependency the script
  reaches, plus the screen and figure dimensions it was produced at.
- **Looking back** — `LaBGAScore_prov_resolve_retrospective` reconstructs the same record
  from each artifact's embedded date and each clone's git reflog, covering the `.mat` files
  the `prep_` scripts write as well as the reports the others publish. It has been run over
  `proj_cfs` and `proj_discoverie`.

Two things follow for anyone editing these scripts:

- **`publish()` captures figures from the screen**, so a figure larger than the X2go
  session is captured at display size. `plugin_set_figure_size` fits the request to the
  display, preserving aspect ratio — see the note in
  [`CLAUDE.md`](CLAUDE.md) and the recommended X2go settings in
  [`LaBGAS_fMRI_analysis_workflow.md`](https://github.com/labgas/LaBGAScore/blob/main/LaBGAS_fMRI_analysis_workflow.md).
- **Renaming a script per study is fine** — the tooling maps a study's renamed copy back
  onto the template it came from, by token overlap on the step designator.

See [`clean/README_provenance.md`](https://github.com/labgas/LaBGAScore/blob/main/clean/README_provenance.md)
for the full guide.

## Out of scope

This README documents only the 4 Group 1 + 15 Group 2 scripts listed above. The rest of `core_scripts_to_run_without_modifying/` and `b_copy_to_local_scripts_dir_and_modify/` — including `prep_3d_run_SVMs_betweenperson_contrasts.m`, the `z_batch_publish_*` orchestration scripts, all `plugin_*` internal helpers, and the many other `b1`/`c3`–`c5`/`d1`–`d15` (except `d10`)/`f2`/`g2`/`h1`–`h3`/`j1`/`k1`–`k2` scripts — is generic CANlab machinery LaBGAS does not currently document here. See `list_of_scripts_and_workflow.m` for the full menu and `a0_begin_here_readme.m` for the generic template's own walkthrough.

---

## A note on TFCE

`doTFCE` in `prep_3a_run_second_level_regression_and_save.m` runs **classic
threshold-free cluster enhancement** (Smith & Nichols 2009) through LaBGAScore's
`group_tfce_from_subject_maps`, with a sign-flip permutation null for
`'onesample'` and a label-exchange null for `'twosample'`. Nuisance covariates,
when supplied, are handled by Freedman-Lane.

The searchlight path (`prep_3c_run_SVMs_on_contrasts_masked.m` →
`c2_SVM_contrasts_masked.m`) uses the same classic TFCE, via
`searchlight_disti_Lukas` in CanlabCore. Both paths therefore run one algorithm.

**Results produced before the 2026 LaBGAScore TFCE overhaul should be re-run.**
That work found the TFCE stack calling pTFCE — a different algorithm — with
mismatched arguments, and repaired both permutation schemes. The two-sample
null had been effectively degenerate (permuting the residuals and the group
labels by the same index changes nothing), and one-sample with covariates could
never reject, because the intercept sat in the nuisance design and
Freedman-Lane returned the effect of interest to every permuted dataset.
