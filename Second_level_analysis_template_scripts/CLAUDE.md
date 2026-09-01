# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A MATLAB template-script framework (part of `CANlab_help_examples`) for second-level (group) fMRI analysis on beta/contrast images. This working copy is a **LaBGAS lab fork** of the generic CANlab template, with a local customization layer on top.

This is not a conventional software package: there is no build system, linter, or automated test suite, and none should be added. It's a curated collection of runnable/copyable `.m` scripts, most of which produce figures/tables and can be run via MATLAB's `publish()` to generate timestamped HTML reports.

## Current objective in this repo

The active goal of work here is to write a single, extensive `README.md` (currently missing — the repo only has `a0_begin_here_readme.m`, a MATLAB-comment readme covering the generic, non-LaBGAS-specific CANlab template) documenting **how LaBGAS actually uses this framework**, including its dependency on the sibling `LaBGAScore` repo. A second goal is to improve the documentation in the Matlab help header (commented section in the beginning) of the scripts below.

A third goal is to optimize the *code* of these scripts — not just their header documentation — for more consistent `publish()` HTML output across the team. These scripts run on a shared lab server accessed via X2go from team members' own computers with differing screen resolutions. The prevailing pattern of creating a figure and then calling `set(gcf,'WindowState','maximized')` before `snapnow`/`saveas` ties the figure's pixel dimensions to whichever X2go client's screen happened to be active — and since these scripts use exclusively default, point-based font sizes (no explicit `FontSize`/`FontUnits` calls), that makes text in published figures look inconsistently too large or too small depending on who ran the script. The fix is to replace `WindowState maximized` with a size set in INCHES (not pixels — MATLAB font sizes are in points, a physical unit, so an inch-anchored canvas keeps the font-to-canvas ratio constant across sessions whose DPI differs) via a shared `plugin_set_figure_size.m` helper (`core_scripts_to_run_without_modifying/`), applied script by script. Each script's new approach is validated first on a `test_figs_`-prefixed copy (e.g. `test_figs_prep_3a_run_second_level_regression_and_save.m`) before being rolled into the canonical script. Started with `prep_3a_run_second_level_regression_and_save.m`.

**Status (2026-09-01):** the cross-machine validation has NOT been run yet, and `test_figs_prep_3a_...` is deliberately excluded from `DEPENDENCIES.md` until it has. Note that `plugin_set_figure_size.m` was revised on 2026-08-31, so the test harness now exercises different behaviour than when it was written:

- The requested size is a **maximum**, not a fixed value. `publish()` captures what is on screen, so a figure larger than the display was previously captured at display size *and at the wrong aspect ratio*, silently — `get(fh,'Position')` still reported the requested size. Measured on the LaBGAS server (1718x1360 at 133 DPI), the old fixed 16x10 in default needed 2128x1330 px, did not fit, and came out 1718x1254: aspect 1.37 instead of 1.60, i.e. exactly the `maximized` behaviour it exists to avoid. Both dimensions are now scaled by one factor, so the aspect ratio always holds, and the window is repositioned fully on screen.
- The **default changed from 16x10 to 12x7.5 in** (same 16:10 aspect). 16x10 is unreachable on any lab laptop — it would need 72 DPI on a 1366x768 client. 12x7.5 is reachable on every lab screen at 96 DPI. A default nobody can achieve guarantees the inconsistency this work is meant to remove.
- Run `LaBGAScore_check_display` (LaBGAScore `clean/`) to see what a given session can produce. `LaBGAScore_prov_publish` records the session's screen size, DPI and the resulting figure dimensions in every report, and flags figures whose size was set by the display — which is what makes the cross-machine comparison measurable when someone does run it. Recommended X2go settings per screen are in `LaBGAS_fMRI_analysis_workflow.md`.

**This documentation effort is scoped to a specific subset of scripts only** — not the full toolbox (~80 scripts). Do not analyze or document scripts outside this list unless explicitly asked to expand scope:

**Group 1 — user-edited entry points, in `b_copy_to_local_scripts_dir_and_modify/`:**
- `a_set_up_paths_always_run_first.m`
- `a2_set_default_options.m`
- `prep_1_set_conditions_contrasts_colors.m`
- `prep_1b_prep_behavioral_data.m`

**Group 2 — core analysis scripts, in `core_scripts_to_run_without_modifying/`.** These are exactly the scripts referenced by the ALL-CAPS `%%` section headers in `a2_set_default_options.m` (per that file's own convention: *"If the title of the section below is capitalized, the scripts and their options have been revamped by @lukasvo76 already"*). Two lowercase (non-revamped) headers in that file — `prep_3d_run_SVMs_betweenperson_contrasts` and the `z_batch_publish_*` options — are **excluded**, along with every other script in `core_scripts_to_run_without_modifying/` not listed below:

| Script | Role |
|---|---|
| `prep_2_load_image_data_and_save.m` | Loads first-level beta/con images into `fmri_data_st` objects, QC + z-scoring, saves to `.mat`, publishes HTML report |
| `prep_3_calc_univariate_contrast_maps_and_save.m` | Calculates contrast images from prep_2 condition images, l2norm-rescales, QC, saves |
| `prep_3a_run_second_level_regression_and_save.m` | Group-level regression per contrast/condition, voxel-wise (`regress()`) or parcel-wise (`robfit_parcelwise()`); optional Bayes Factor conversion and MVPA regression on covariates |
| `c2a_second_level_regression.m` | Displays/thresholds results from `prep_3a_...` |
| `prep_3c_run_SVMs_on_contrasts_masked.m` | Cross-validated SVM per contrast, masked; saves results |
| `c2_SVM_contrasts_masked.m` | Displays/thresholds SVM results from `prep_3c_...` |
| `prep_3f_create_fmri_data_single_trial_object.m` | Builds an `fmri_data_st` single-trial object from single-trial con images (produced by `LaBGAScore_firstlevel_s2_fit_model.m`), attaching ratings/VIFs metadata |
| `prep_3g_create_fmri_data_runwise_contrast_object.m` | Builds an `fmri_data_st` object of runwise contrasts from condition betas, attaching runwise phenotype metadata |
| `c2f_run_MVPA_regression_single_trial.m` | MVPA regression (default PCR) on a continuous outcome, on the single-trial object from `prep_3f_...` |
| `c2g_run_multivariate_mediation_single_trial.m` | Multivariate mediation (PDM) analysis on a continuous outcome, on the single-trial object from `prep_3f_...` |
| `prep_4_apply_signatures_and_save.m` | Applies selected CANlab signature patterns to conditions/contrasts, saves to `DAT.SIG_conditions`/`DAT.SIG_contrasts` |
| `d_signature_responses_generic.m` | Plots and tests signature responses from `prep_4_...` |
| `d10_signature_riverplots.m` | Riverplots of signature responses (cosine similarity) from `prep_4_...`; works only on signature *groups*, not individual signatures |
| `h_signature_responses_group_diff.m` | Group comparison of signature responses (cosine similarity) from `prep_4_...` |
| `e1_corr_patterns.m` | Pairwise searchlight correlation maps between condition/contrast images (`searchlight_correlation()`) |

Minor naming mismatches worth noting in the eventual README (not bugs to fix): `c2g_run_multivariate_mediation_single_trial.m`'s own section header in `a2_set_default_options.m` says "MULTILEVEL_MEDIATION"; `e1_corr_patterns.m`'s internal `%%` title says `e1_corr_patterns_conds.m`.

## Script-category convention (background)

- `prep_*` — one-time setup: load data, build contrasts, apply signatures/parcellations.
- Lettered scripts (`a_`/`b_`/`c_`/…/`k_`) — on-demand analysis/reporting, runnable in any order once prep has run.
- `z_batch_*` — orchestration; `*publish*` variants render timestamped HTML into `results/published_output`.
- `plugin_*` — internal helpers, not meant to be run or edited directly.

## Group 1 details and the shared data model

`a_set_up_paths_always_run_first.m` must be run first; it auto-calls `a2_set_default_options.m` and depends on two functions from the sibling **LaBGAScore** repo (`LaBGAScore_prep_s0_define_directories`, `LaBGAScore_firstlevel_s1_options_dsgn_struct`) — see below. `prep_1_set_conditions_contrasts_colors.m` defines `DAT` (conditions, contrasts, colors) and calls `a_set_up_paths_always_run_first` itself. `prep_1b_prep_behavioral_data.m` is optional and attaches behavioral/grouping data to `DAT.BEHAVIOR` / `DAT.BETWEENPERSON`.

Data model shared across all in-scope scripts:
- `DAT` struct — conditions, contrasts, contrastnames, colors, behavioral/group data, signature results.
- `DATA_OBJ` / `DATA_OBJ_CON` — cell arrays of CanlabCore `fmri_data`/`fmri_data_st` objects (per-condition / per-contrast).
- Persisted to `image_names_and_setup.mat` and `data_objects.mat`. Scripts from `prep_3a` onward, and the lettered on-demand scripts, **automatically reload these saved `.mat` files internally** — `b_reload_saved_matfiles.m` is not a required explicit step before every in-scope script, only earlier in the generic workflow.

`a2_set_default_options.m` centralizes the default options (masks, thresholds, scaling, cross-validation settings, etc.) consumed by each Group 2 script, organized into one section per script (capitalized section title = script is in scope / revamped for LaBGAS use).

## LaBGAScore dependency

Several in-scope scripts depend on a separate sibling repo, `LaBGAScore` (local path `/data/master_github_repos/LaBGAScore`; upstream `github.com/labgas/LaBGAScore`), for:
- Directory setup: `LaBGAScore_prep_s0_define_directories`
- First-level design options: `LaBGAScore_firstlevel_s1_options_dsgn_struct`
- Fitting first-level models (source of single-trial con images consumed by `prep_3f_...`): `LaBGAScore_firstlevel_s2_fit_model.m`
- Atlas/ROI mask generation referenced in `a2_set_default_options.m`: `LaBGAScore_atlas_binary_mask_from_atlas.m`, `LaBGAScore_atlas_rois_from_atlas.m`

- Parallel pool setup before bootstrap/permutation: `LaBGAScore_smart_parallel_pool_setup.m` (called from `c2a_second_level_regression.m`, `prep_3a_...`, `prep_3c_run_SVMs_on_contrasts_masked.m`)
- Group TFCE: `group_tfce_from_subject_maps.m` (called from `prep_3a_...`)
- Thresholded `fmri_data` from a `statistic_image`: `thresholded_fmri_data_from_statistic_image.m` (called from `prep_3a_...`, `c2_SVM_contrasts_masked.m`)

The last three were found by LaBGAScore's dependency tooling and were previously undocumented here.

Document this dependency at a high level (what's called and why) rather than diving into LaBGAScore's own internals.

`DEPENDENCIES.md` in THIS folder (not the repo root — it documents this folder, so it lives here) is the **generated**, authoritative version of the above, produced by `LaBGAScore_dep_report`. It covers exactly the 19 scripts in README.md's Script reference (4 Group 1 + 15 Group 2), not all ~113 in this folder. Regenerate with the file list from that table; do not hand-edit it, `dependencies.tsv` or `dependencies.yml`.

Provenance — which commit of CanlabCore et al. produced a given result — is recorded by LaBGAScore's `clean/LaBGAScore_prov_*` tooling, not by anything here. See `clean/README_provenance.md` in LaBGAScore.

## Where to look for more detail

`a0_begin_here_readme.m` and `list_of_scripts_and_workflow.m` cover the full generic CANlab toolbox (all ~80 scripts, out of scope for the current documentation effort) — useful background on conventions, but don't duplicate their content wholesale.
