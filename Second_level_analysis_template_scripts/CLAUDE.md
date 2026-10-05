# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A MATLAB template-script framework (part of `CANlab_help_examples`) for second-level (group) fMRI analysis on beta/contrast images. This working copy is a **LaBGAS lab fork** of the generic CANlab template, with a local customization layer on top.

This is not a conventional software package: there is no build system, linter, or automated test suite, and none should be added. It's a curated collection of runnable/copyable `.m` scripts, most of which produce figures/tables and can be run via MATLAB's `publish()` to generate timestamped HTML reports.

## State of the work in this repo

**Documentation.** `README.md` exists and documents how LaBGAS actually uses this
framework, including the dependency on the sibling `LaBGAScore` repo. Script header
comments (USAGE/OPTIONS/NOTES) are maintained alongside the code; when you add an option,
document it in the header of the script that consumes it, in the same `% * name` style.

**Figure consistency across sessions.** The prevailing pattern of `set(gcf,'WindowState',
'maximized')` tied a figure's pixel dimensions to whichever X2go client happened to be
connected, and since these scripts use exclusively default, point-based font sizes, that
made published text inconsistently large or small per person. It is replaced by
`plugin_set_figure_size.m` (which now lives in **LaBGAScore** `figures/`, since first level uses it too), which sizes in INCHES (not pixels: MATLAB font sizes are in
points, a physical unit, so an inch-anchored canvas keeps the font-to-canvas ratio constant
across sessions of differing DPI). Points worth not relitigating:

- The requested size is a **maximum**, not a fixed value. Interactively, `publish()`
  captures what is on screen, so a figure larger than the display was previously captured
  at display size *and at the wrong aspect ratio*, silently — `get(fh,'Position')` still
  reported the requested size. Both dimensions are now scaled by one factor, so aspect
  always holds, and the window is repositioned fully on screen.
- The **default is 12x7.5 in**, not 16x10 (same 16:10 aspect). 16x10 is unreachable on any
  lab laptop; a default nobody can achieve guarantees the inconsistency this exists to
  remove.
- `plugin_set_figure_size` also takes `'fig'` (size a figure a drawing call opened, since
  `canlab_results_fmridisplay 'multirow'`, `@region/montage` and `histogram(...,'byimage')`
  all open their own), `'titlescale'` (default 2/3), `'keepaspect'`, and `'minpanel'`
  (grow the canvas until each panel of a grid is legible, bounded at 20x30 in and 3:1).

**Headless is now the default way to run these scripts**, via LaBGAScore's
`clean/labgascore_run_headless.sh`. Consequences for anyone editing here:

- Headless, `publish()` **prints** figures rather than capturing them, so figure size is
  NOT limited by the 1024x768 / 72 dpi virtual screen. The screen-fitting logic is lifted
  headless only when `'minpanel'` is given, so no other figure changes size.
- `matlab -batch` **cannot** `publish()`. Use `-nodisplay` with `-r`, stdin from
  `/dev/null`.
- `LaBGAScore_check_display` and the X2go DPI table apply to the interactive route only.

**`publish()` swallows errors — this is the main hazard in this repo.** It catches a
script's error into the html and returns normally: no exception, exit status 0, report file
present. Chains of these scripts have repeatedly appeared to complete while one had died,
sometimes after an hour of computation and before anything was saved. Every failure found
in this repo recently was of that shape: unguarded optional variables (`idx_nuisance`,
`bayesian_results`, `tfce_results`, `summary()`), parcelwise branches assuming
Bayes/TFCE results exist, and `char()` on a multi-element `wh_interest` dying in `horzcat`.
So:

- Run scripts through `LaBGAScore_run_reports` / `labgascore_run_headless.sh`, which read
  the report back and fail on `<pre class="codeoutput error">`, and assert the expected
  `.mat` exists. Never trust a chain that merely "finished".
- When you add code that reads an option variable, assume it may not exist. `exist(...,
  'var')` guards are why several of these scripts now survive configurations they used to
  crash on.
- Detect errors from report **markup**, never by searching report text for "Error in" —
  that phrase occurs in ordinary comments, which `publish` renders as prose.

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
| `c2h_run_multivariate_mediation.m` | Single-level multivariate mediation (PDM) on second-level CONTRAST images: X = group from `DAT.BETWEENPERSON.group`, M = the subject's contrast image, Y = a between-person outcome. The single-level counterpart of `c2g_...`, needing no `prep_` step of its own because contrast images are already part of the standard pipeline |
| `prep_4_apply_signatures_and_save.m` | Applies selected CANlab signature patterns to conditions/contrasts, saves to `DAT.SIG_conditions`/`DAT.SIG_contrasts` |
| `d_signature_responses_generic.m` | Plots and tests signature responses from `prep_4_...` |
| `d10_signature_riverplots.m` | Riverplots of signature responses (cosine similarity) from `prep_4_...`; works only on signature *groups*, not individual signatures |
| `h_signature_responses_group_diff.m` | Group comparison of signature responses from `prep_4_...`, unadjusted and optionally covariate-adjusted; corrects across the signature family with the methods named in `corrections_i_want` and writes a summary table plus violin panels |
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

`DEPENDENCIES.md` in THIS folder (not the repo root — it documents this folder, so it lives here) is the **generated**, authoritative version of the above, produced by `LaBGAScore_dep_report`. It covers exactly the 20 scripts in README.md's Script reference (4 Group 1 + 16 Group 2), not all ~113 in this folder. Regenerate with the file list from that table; do not hand-edit it, `dependencies.tsv` or `dependencies.yml`.

Provenance — which commit of CanlabCore et al. produced a given result — is recorded by LaBGAScore's `clean/LaBGAScore_prov_*` tooling, not by anything here. See `clean/README_provenance.md` in LaBGAScore.

**Adapting these templates is where studies actually go wrong**, and the catalogue of how
lives in LaBGAScore too: *"Ten traps when adapting a template"* in
`LaBGAS_fMRI_analysis_workflow.md`. Every entry is a real wrong-but-clean run. Two of them
are now automated by checkers in LaBGAScore's `clean/`, which should be run on a study's
model script directory before any long chain:

- `use_before_def.py` — an option read ABOVE the line that defines it (the script dies
  late, after the expensive work, having saved nothing).
- `set_after_use.py` — an option set BELOW the line that already consumed it (nothing
  errors; the default silently wins). Advisory, not a gate.

`checkcode` reports zero messages on either checker's positive control, which is why they
exist alongside it. Three traps are specific to scripts in THIS folder and worth knowing
before you edit one:

- **`covs2use` also gates `roi_means_table`** in `prep_3a_...`, not just the design. A
  `prep_3a` run that exists only to generate features for the PLS-DA / Elastic Net pipeline
  must therefore name the covariates those pipelines residualise, or they fail with
  *"covariate_names not found in roi_stats table"*.
- **A constant column in a `custom` design** is read as a manual intercept: `prep_3a`
  prints *"Skipping this contrast"*, saves EMPTY results and exits 0. This bites whenever a
  subject filter makes a site dummy constant.
- **Region-table peaks are clipped at 7.0345 unless the vendored table is used.** CanlabCore's
  `@region/table` ends `get_signed_max` with `maxZ = norminv(1 - 1E-12)` and clips every value
  above it, not only the infinities its own comment describes. `c2a` passes a table function
  handle to `LaBGAScore_region_table_safe`; until v8.5 the FDR, uncorrected and Bayesian
  branches passed `@table` and only the TFCE branches passed `@LaBGAScore_region_table`. The
  giveaway is a column of identical `7.0345` values. Inference is never affected - it is the
  printed peak only - but the column is worthless above the ceiling, which t-maps and
  especially Bayes factor maps (stored as `2*ln(BF)`) do exceed. **For BF maps the two
  defects compound:** a 7.0345 ceiling on a `2*ln` scale caps reported evidence at
  BF10 = exp(3.517) ~ 34, so a real peak of BF10 = 39591 printed as 7.03.
- **The Bayes column was also MISLABELLED until 2026-09-25.** `estimateBayesFactor` sets
  `.type = 'BF'`, and the table builds its header as `['max' Z_descrip]`, so the column read
  `maxBF` while holding `2*ln(BF10)` - a printed 15.05 is a Bayes factor of **1854**, not 15.
  `LaBGAScore_region_table` now emits `max_2lnBF` plus `maxBF10 = exp(max_2lnBF/2)`.
  **Convert with `exp(v/2)`, never `exp(v)`:** the wrong conversion overstates evidence FOR
  THE NULL, turning model_2c_IOM's true "83% moderate, median BF10 0.16" into a spurious
  "82% strong, median 0.027". The thresholds themselves were always right - `c2a` uses
  `2*log(BF_threshold_glm)` and `prep_3a` hardcodes `2.1972`, which is `2*ln(3)` (labelled
  "|BF| > 3") and NOT ln(9), the misreading that caused this.
- **The JZS Bayes factor has a sample-size floor.** `t1smpbf(0, n)` bounds attainable
  evidence for H0: n=70 floors at BF10 = 0.131 (7.6:1), n=93 at 0.115, n=158 at 0.089. Below
  roughly n=100, **BF10 < 1/10 is unreachable no matter how null the data are**, so "0% strong
  evidence for the null" in a small sample is a design ceiling, not a weak result. Say so when
  reporting it.
- **`prep_2` and `prep_1b` must agree on whose sample the design describes.** `prep_2`
  subsets `DAT.BETWEENPERSON.group` itself but never touches
  `.conditions{}`/`.contrasts{}`. Since v2.6 it decides by length and errors with both
  counts named; before that, a model whose sample was itself a subset (patients only, one
  site only) either crashed inside `prep_2` or mis-aligned silently.

## The MVPA paths and `@predictive_model` (done for prep_3a/c2a, Sept 2026)

`fmri_data.predict` returns **no inferential quantity at all** for a continuous
outcome - `pred_outcome_r`, `mse`, `rmse`, `meanabserr`, `cverr`, and nothing
else. `prep_3a` supplied the missing pieces itself (`cv_seed_mvpa_reg_cov`,
`cv_strata_mvpa_reg_cov`, `nperm_mvpa_reg_cov`, `numcomponents_mvpa_reg_cov`).
That stopgap is now backed by a real migration for the `domvpa_reg_cov` path.
**`prep_3c`'s SVM still calls `predict` and has NOT been migrated.**

### `mvpa_engine` - three engines, and why the choice matters

| value | fit | inner CV that tunes the hyperparameter |
|---|---|---|
| `'legacy'` | `fmri_data.predict` | `estimateparams`: round-robin over ROW INDEX |
| `'predictive_model'` | `@predictive_model` via `mvpa_reg_cov_predictive_model` | `estimateparam`: the same round-robin |
| `'tuned_nested'` (**default** since 2026-10-05) | `mvpa_reg_cov_tuned_nested` | inner folds REBUILT from the training subset's strata, per outer fold - the pattern CanlabCore's own tutorials teach |

`'legacy'` was the default until 2026-10-05, so **a null MVPA result from an
earlier run may be the engine, not the data.** The default now comes as a
coupled SET (lassopcr + the `lasso_num` grid + `'strata'` + two distinct seeds),
documented in `prep_3a`'s header under *THE RECOMMENDED CONFIGURATION* and
shipped in `a2`; `cv_strata_mvpa_reg_cov` is the one member a study must supply,
and `tuned_nested` errors rather than fitting untuned without it.

**This is not a stylistic choice.** Measured on proj_discoverie model_2k
immune_PC1, identical data, identical outer folds, 149154 grey-matter voxels:

| engine | r |
|---|---|
| `tuned_nested` | **+0.1812** |
| ooFmriDataObjML (reference cross-check) | +0.1186 |
| legacy / `estimateparams` | +0.0158 |
| `predictive_model` / `estimateparam` | +0.0158 |

Engines that rebuild the inner splitter under the outer folds' structural
constraints find signal; those tuning by a round-robin over row index, blind to
that structure, find essentially none. The two `estimateparam` rows agree to four
decimals because they are one engine behind two front-ends. Do NOT rank +0.1186
against +0.1812 - inner-fold seed alone moves the estimate by ~0.06.

### Two traps that cost real time here

**MASK THE FEATURES.** `prep_3a` built the MVPA design from the UNMASKED
`cat_obj`. The univariate branch masks the STATISTIC IMAGE after fitting, which
is correct there because each voxel's test is independent; it is wrong for MVPA,
where the decomposition runs over every included voxel. 235807 voxels against
149154 in grey matter - 36% white matter, CSF and edge. Now behind
`domask_mvpa_reg_cov` (default true). The TFCE branch already did this correctly,
so the template encoded the principle; `mvpa_reg_cov` was the one path that missed it.

**ONE SEED CANNOT DRIVE TWO PARTITIONS.** `cv_seed_mvpa_reg_cov` seeded both the
outer CV partition and `tuned_nested`'s inner folds, so a run that used different
values for each could not be reproduced from the options. Now
`tuned_seed_mvpa_reg_cov`, defaulting to `cv_seed_mvpa_reg_cov`. The symptom was
quiet: the fit still ran and looked reasonable, landing on a different modal
hyperparameter and r = +0.1998 instead of +0.1812.

### Pattern inference lives in c2a, not prep_3a

`prep_3a` answers *is the model better than chance* (permutation test). Only if
that is significant is *which voxels does it lean on* worth paying for, so
bootstrap and stability selection are in `c2a`, opt-in, behind
`dobootstrap_mvpa_reg_cov`. `mvpa_engine` there selects which bootstrap runs -
the two are EXCLUSIVE, never both.

c2a prefers the object `prep_3a` saved as `mvpa_stats.pm` and bootstraps THAT,
rather than refitting and possibly landing on a different penalty. A legacy fit
has no `.pm` and is rebuilt from the same algorithm and folds, which the report
states. Folds are recovered from `teIdx`, a CELL of `nfolds` logical `[n x 1]`
masks - not a matrix, so it cannot be collapsed by multiplication.

**Bootstrap p can collapse.** On a strongly regularised model the weights are
near-identical across resamples, the empirical p floors at `2/(nboot+1)` for
every voxel, and the FDR mask becomes meaningless. c2a counts voxels at the floor
and says so. Stability selection is the recommended inference there.

**Stability selection: k and the threshold are coupled.** Meinshausen & Buhlmann
(2010) bound `E(V) <= k^2/((2*pi-1)*p)`, valid only for `pi > 0.5`, and `k` enters
SQUARED where `p` enters linearly. The class default `k = 2000, pi = 0.6` controls
nothing at brain scale - at `p = 149154` it bounds E(V) at 134, and no valid `pi`
rescues `k = 2000`. c2a inverts the formula instead,
`k = sqrt((2*pi-1)*p*E(V))`, giving **k = 345 at pi = 0.9, E(V) = 1**. Options:
`stab_threshold_mvpa_reg_cov` (0.9), `stab_EV_mvpa_reg_cov` (1),
`stab_k_mvpa_reg_cov` (empty = derive). The implied E(V) is always printed.
Caveat: the theorem assumes SUBSAMPLING at n/2; `stability_selection` resamples
WITH REPLACEMENT (Shah & Samworth 2013 give that case), so the number is the right
order, not an exact guarantee.

**Speed.** `stability_selection` refits per resample in a SERIAL loop; `bootstrap`
does the same resampling and RETAINS every weight vector in `pm.weights.boot_w`.
So stability is one sort per column of a matrix already in memory -
`stab_reuse_boot`, default true. At 1.02 s per fit on 93 x 149154 that is seconds
instead of ~85 min at `nstab = 5000`.

### The helper functions that ship alongside

Five files in `core_scripts_to_run_without_modifying/`, called BY `prep_3a` and
`c2a` rather than run directly, so they are not in the script table:

| file | role |
|---|---|
| `mvpa_reg_cov_predictive_model.m` | `mvpa_engine = 'predictive_model'`. Wraps crossval / permutation_test / bootstrap / weight_map_object. Returns `[pm, stab]` - stability comes back as a SECOND OUTPUT, not attached to `pm`, because `pm`'s diagnostics property is protected |
| `mvpa_reg_cov_tuned_nested.m` | `mvpa_engine = 'tuned_nested'`. Inner grid search rebuilt from the training subset's strata per outer fold. Also carries its own permutation test, which permutes the WHOLE nested procedure because the tuning is part of what is being tested |
| `mvpa_reg_cov_stability_from_boot.m` | Derives stability selection from `pm.weights.boot_w` instead of refitting. Returns a struct |
| `mvpa_reg_cov_oofmri.m` | ooFmriDataObjML engine, kept as a reference cross-check rather than a default - that package is UNMAINTAINED (last commit 2024-08-23) |
| `mvpa_reg_cov_benchmark_algorithms.m` | Compares algorithms under identical folds |

Two traps in `mvpa_reg_cov_oofmri` worth not rediscovering. `get_r` returns
**1 - r**, a LOSS for `gridSearchCV` to minimise, not a correlation - reading
`cvGS.scores` as r gives impossible values above 1. And `crossValScore` AVERAGES
per-fold scores while the other engines POOL all held-out predictions; those are
different statistics, because correlation is not an additive loss the way MSE is.
The wrapper now reports the pooled r, with per-fold values as a diagnostic and a
Fisher-z weighted mean as the defensible average. CanlabCore's `@pipeline` also
SHADOWS ooFmriDataObjML's, so that package must be added to the path LAST.

### predictive_model's fitted state is PROTECTED

`diagnostics`, `weights`, `fitted_values` and most other fitted state sit under
`properties (SetAccess = protected)`. Only class methods may write them; an
external function assigning `pm.diagnostics.stability_selection` fails at RUNTIME
with *"Unable to set the 'diagnostics' property ... because it is read-only"* -
`checkcode` passes, so nothing catches it until a job is hours in. The helpers
here therefore RETURN structs rather than mutating `pm`, and build frequency maps
by copying `pm.weights.weight_obj` and swapping its `.dat`.

### FDR: voxel and parcel level use BH, deliberately

`threshold(obj, q, 'fdr')` calls CanlabCore's `FDR.m`, which is plain
**Benjamini-Hochberg** (`pID`, the independence/PRDS form) with no pi0 estimation.
Every spatial threshold - voxelwise, parcelwise, TFCE, Bayes, MVPA weight maps -
flows through it, so they are already uniform. Storey appears ONLY in the
small-`m` tabular contexts (8 ROIs, 30 neurotransmitter maps) via
`LaBGAScore_Storey_FDR`, reported alongside BH rather than instead of it.

Do not "upgrade" the spatial path to Storey. The binding constraint at voxel
scale is DEPENDENCE, not m: BH controls FDR under positive regression dependency,
which is the standard justification for smooth neuroimaging data, whereas
Storey's pFDR assumes independence. A well-identified pi0 buys nothing if the
control it feeds is not valid under the dependence actually present.

### Still open

`prep_3c`'s SVM remains on `fmri_data.predict`. `@predictive_model` would also
remove two live hazards there: its `svr` runs on MATLAB's `fitrsvm` rather than
the unmaintained Spider copy vendored in `CanlabCore/External/spider`, and its
`pcr`/`lassopcr` take `{'numcomponents', k}` properly - see the `cv_pls` trap
below.

### The cv_pls trap (fixed, do not re-introduce)

`predict()`'s `cv_pls` calls `plsregress(X,Y)` with no `ncomp` when
`'numcomponents'` is absent, and MATLAB then uses the **maximum**,
`min(n-1, p)`. For PLS the component count *is* the regularisation, so the
default is no regularisation: the model interpolates the training fold and
generalises at chance. Measured on synthetic data with the signal in the leading
component, n = 60, 5-fold:

```
numcomponents   1     2     3     5    10    20   default(max)
r            0.70  0.72  0.73  0.73  0.73  0.73        -0.12
cv_pcr reference r = 0.73
```

Regularised PLS equals PCR from 3 components up; the default is worse than
chance. `prep_3a` now refuses `cv_pls` unless `numcomponents_mvpa_reg_cov` is
set. LaBGAScore's own ROI PLSR pipeline was checked and is clean - it always
passes an explicit `lv`, bounded by `capLV.m` and chosen by nested inner-fold CV.

A second finding from the same testing, worth keeping: at n = 60-90, when the
predictive direction is orthogonal to the dominant variance components, **no**
algorithm recovers it - `cv_pcr` r = 0.13, `cv_pls` -0.12, `cv_svr` 0.16, none
significant. A null from voxel-wise MVPA at these sample sizes therefore does not
license "no diffuse association exists".

## Where to look for more detail

`a0_begin_here_readme.m` and `list_of_scripts_and_workflow.m` cover the full generic CANlab toolbox (all ~80 scripts, out of scope for the current documentation effort) — useful background on conventions, but don't duplicate their content wholesale.

## Start from the repo, not from a sibling study

**Always build a study's scripts from the most recent version in the repo.** Use
an older study's copies only as an *example* of which study-specific adaptations
to make, never as the thing you copy and rename.

**Why:** study copies are frozen when they are made and then drift. Measured on
2026-10-01, `proj_bitter-reward`'s second-level copies were behind the current
templates by 66 lines (`a_set_up_paths`), 144 (`a2_set_default_options`), 460
(`prep_2`) and 628 (`prep_3`) - and the last two are
`core_scripts_to_run_without_modifying`, so more than half the current script was
simply missing. Two fixes already upstream were absent, including the
`scriptsdir`/`modelname_2nd` correction that sends a variant model's scripts dir
at the wrong model.

**How to apply:** copy the template, then port only the genuinely
study-specific values from the old copy. Distinguish those from stale defaults by
checking the template's own history - `git log -L '/^<option> /,+1:<file>'`. On
that same date, of five option differences in `a2_set_default_options`, two
(`atlasname_glm = 'canlab2023_fine_2mm'`, `similarity_metric_sigs = 'dotproduct'`)
were former template defaults rather than study choices, and only three
(`doroi_analysis`, `keyword_sigs`, `subjs2exclude_data`) were real.
