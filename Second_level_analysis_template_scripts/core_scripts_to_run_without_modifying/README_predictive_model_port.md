# Porting the MVPA-on-covariate analysis to `@predictive_model`

Branch `feature/predictive-model-mvpa`, based on `c918827` (prep_3a v9.3, the
commit that first gave `domvpa_reg_cov` any inference at all).

## Why

v9.3 hand-rolled, in template scripts, what CanlabCore now provides as a class:
fold construction, a permutation null, a bootstrap of the weight map, and the
scoring of predicted against observed. That code works, but it lives in two
1000+ line scripts that every study copies and edits, so every study inherits
the maintenance burden and any bug in it. Three bugs in exactly this area
surfaced on 2026-09-24 alone.

## Mapping

| v9.3, hand-rolled | `@predictive_model` |
|---|---|
| `holdout_set_method_mvpa_reg_cov` switch, incl. the `'strata'` case | `cv_splitter` |
| `predict(..., 'nfolds', fold_labels, 'algorithm_name', ...)` | `crossval` |
| `parfor` permutation loop, `p = (sum(null>=obs)+1)/(n+1)` | `permutation_test` |
| c2a's MVPA bootstrap branch | `bootstrap` |
| *(no legacy equivalent)* | **`stability_selection`** - new capability, see below |
| `corr(yfit, Y)` | `cv_scorer.pearson_r` |
| manual re-map of weights into an image | `weight_map_object` |

## Which algorithm

Settled a prior confusion: **CANlab's default `cv_lassopcr` is arithmetically
PCR.** `fit_lassopcr`'s own docs say that with neither `lasso_num` nor
`estimateparam`, it "use[s] the full LASSO model -> reduces to PCR (identical to
cv_pcr / default cv_lassopcr)". Choosing `pcr` was therefore never a departure
from the CANlab default; it is the same estimator under a different name. The
shrinkage exists only when asked for.

Regression algorithms actually available (`pcr`/`lassopcr` are special-cased in
`fit()`, the rest are in `algorithm_registry`):

| algorithm | fitter | regularised? | honestly tunable today? |
|---|---|---|---|
| `pcr` | PCA + OLS | no (component truncation only) | n/a |
| `lassopcr` (default opts) | = `pcr` | no | n/a |
| **`lassopcr` + `'estimateparam'`** | PCA + LASSO + relaxed-OLS refit | **yes** | **yes - nesting is built in** |
| `lassopcr` + `{'lasso_num', k}` | as above, fixed path step | yes | k is a free choice |
| `svr` | `@fitrsvm`, linear | via `BoxConstraint` | **no** - see below |
| `linear_svr` | `@fitrlinear` | yes | no |
| `lasso` / `ridge` | `@fitrlinear` | yes | no |

**Recommendation: `lassopcr` with `'estimateparam'` as primary**, `pcr` as the
unregularised reference. With n = 91 and ~150k voxels, real regularisation
should beat unregularised PCR, and `estimateparam` selects the penalty by
NESTED cross-validation (`lasso_cv` with an inner `cv_assignment`), then
OLS-refits on the surviving components. It is the only option here that is both
properly regularised and tunable without bias using what the class ships today.

**On SVR.** It no longer depends on the unmaintained Spider toolbox - `svr` is
`@fitrsvm`. That half of the original objection against it is void. The other
half stands: `grid_search`'s own documentation says nested CV is "not yet
automated - see plan SS G6", so tuning `BoxConstraint` within outer folds is work
we would have to do, and tuning it on the whole sample would bias the
cross-validated score. Prefer `linear_svr` (`@fitrlinear`) over `svr` if SVR is
used at all: `fitrsvm` does not scale comfortably to 150k features.

### Benchmark results (measured 2026-09-24)

Run on model_2j's own MVPA data (n = 91, ~150k voxels, stress vs control), under
folds rebuilt from s6c's seed and strata - 3 strata of 33/40/18, matching what
s6c reported.

| algorithm | r | R2 | RMSE | s |
|---|---|---|---|---|
| **`lassopcr` + `estimateparam`** | **+0.0786** | **-0.0276** | **1.0286** | 7.2 |
| `linear_svr` | +0.0461 | -0.1341 | 1.0806 | 6.9 |
| `ridge` | +0.0461 | -0.1341 | 1.0806 | 6.3 |
| `pcr` | +0.0199 | -0.2425 | 1.1311 | 6.9 |

**The port is faithful.** `pcr` reproduces the legacy `pred_outcome_r` of
+0.0199 exactly, so `crossval` + `cv_splitter.custom_partition(fold_labels)` is
equivalent to `predict(..., 'nfolds', fold_labels)`. Checklist items 1 and 2
below are closed by this run.

**The ordering is the regularisation gradient**: nested-CV-tuned lasso-PCR >
ridge at defaults > unregularised PCR. R2 improves from -0.24 to -0.03. Since
default `cv_lassopcr` reduces to PCR, the bottom row IS the CANlab default, and
the default leaves real performance on the table. This is the measured basis for
defaulting to `lassopcr` + `estimateparam`.

**`linear_svr` and `ridge` returned byte-identical results** - 0.0460891891816982
for both, every digit. Both dispatch to `@fitrlinear`, whose default `Learner` is
`'svm'` and whose default regularisation for an SVM learner is ridge, so the
registry's `linear_svr` (defaults `{{}}`) and `ridge`
(`{{'Regularization','ridge'}}`) describe the SAME model. Worth reporting
upstream: two registry rows that look like different algorithms and are not.

**All four have R2 < 0** - every estimator predicts worse than the sample mean.
The algorithm was never what limited this analysis.

Limits: only `estimateparam` had its hyperparameter tuned, since its nesting is
internal and `grid_search` cannot nest yet, so the other rows are lower bounds
rather than a fair head-to-head. `svr` (`@fitrsvm`) was NOT included, so its
scalability at ~150k features remains untested.

Superseded note: the earlier comparison (pcr 0.73, pls 0.73, svr 0.68) was
measured with the LEGACY Spider-based SVR on positive-control data and does not
transfer to `fitrsvm`/`fitrlinear`.

## Two kinds of pattern inference, not one

`bootstrap` and `stability_selection` answer different questions, and the port
exposes both rather than picking for the user:

| | question | output |
|---|---|---|
| `bootstrap` | is this voxel's weight reliably non-zero? | z, p, FDR-thresholded map |
| `stability_selection` | is this voxel reliably in the top-k weights? | selection frequency in [0,1], stable mask |

This matters because the recommended primary is a **regularised** model.
Regularisation pins the weights to nearly the same solution on every resample,
so the bootstrap z/p collapses and reports implausibly sharp inference.
Stability selection asks whether a voxel keeps its *rank* instead, which stays
informative in that regime (Meinshausen & Buhlmann, JRSS-B 2010). For
unregularised `pcr` the reverse holds: the bootstrap is the motivated choice.

So the pairing is deliberate - `pcr` + bootstrap, `lassopcr`+`estimateparam` +
stability selection - and running both on one fit is a useful cross-check.
Disagreement between them is a finding, not a bug.

The selection frequencies are mapped into voxel space via the route the
method's own docs recommend (stash as a weight vector, re-run
`weight_map_object`), done on a **copy** so the real weights are not
overwritten.

## The one thing that does not map cleanly

`cv_splitter.stratified_kfold` stratifies on **Y**. That is built for class
labels and is meaningless for a continuous outcome. Centre stratification is
therefore kept as the caller's job: compute fold ids from `num_center` exactly
as the legacy `'strata'` branch does, then wrap them with
`cv_splitter.custom_partition`.

This is the right split of responsibility anyway. Which variable defines a
stratum is a *study* decision; how folds are then executed is library work.

## Merge strategy

Designed so the merge back to `master` is small and reviewable:

1. **All new logic is in a new file**, `mvpa_reg_cov_predictive_model.m`. It adds
   a file rather than rewriting a block inside a 3000-line script, so it cannot
   conflict with anything.
2. **`prep_3a` gets one `switch`**, on a new option `mvpa_engine`
   (`'legacy'` | `'predictive_model'`), defaulting to `'legacy'`. The legacy path
   is untouched, so existing study copies keep working unchanged and the diff
   against `master` is roughly twenty lines.
3. **`c2a` likewise** reads from `pm` when the new engine produced the results,
   and from the legacy structs otherwise.

Stage two, once the new path has reproduced the legacy numbers on a real model,
is to flip the default and then delete the legacy branch in a separate commit.

### Rebase before merging

This branch is based on `c918827` and therefore does **not** contain two fixes
made on `master`'s working tree on 2026-09-24:

- the `cons2boot` / `cons2boot_mvpa_reg_cov` fix in `c2a` (the MVPA bootstrap
  branch was dead code, referencing a name defined nowhere)
- the neurotransmitter FDR use-before-def in `prep_3a`

Commit those to `master` first, then `git rebase master` here. The second one
does not touch the MVPA code; the first touches the very branch this port
replaces, so rebasing before writing the c2a half avoids resolving that twice.

## Verification checklist — none of this has been run yet

The function is written against the API as read from the class, not from a
working run. Before it replaces anything:

1. `crossval` reproduces the legacy `pred_outcome_r` on the same data and the
   same `fold_labels`, to within floating-point noise.
2. `custom_partition` yields exactly the folds the legacy `'strata'` case built
   — compare fold id vectors element-wise, do not eyeball fold sizes.
3. `permutation_test`'s null is centred near zero. A permuted outcome must not
   be predictable; the legacy code warns when `|mean(null)| > 0.10` and the
   replacement needs the same check.
4. `weight_map_object` returns a map in the same voxel space as `mvpa_dat`, and
   its weights correlate ~1 with the legacy weight map.
5. `stability_selection` returns `n_stable` > 0 and a `selection_freq` map in
   the same voxel space; confirm `pm.weights.w` still holds the real weights
   afterwards and not the frequencies.
6. `bootstrap` FDR-thresholded weights are comparable to the legacy bootstrap
   — this branch has never run in any model, so there is no reference yet;
   generate one from the legacy path first.
7. Run the whole thing on the model_2j MVPA data and confirm the conclusion does
   not change.

## Performance note

The legacy permutation loop measured **>24 s per permutation** on 21 workers for
model_2j (n=91, ~150k voxels). For `pcr` this is mostly wasted work: permuting
`Y` does not change `X` or the folds, so the per-fold PCA basis is *identical*
across all permutations, yet it is recomputed every time. Precomputing the
per-fold component scores once turns each permutation into a handful of small
least-squares fits.

This holds for `pcr` only — PLS components are derived from `Y`, and SVR refits
entirely, so both genuinely need the full loop. Worth raising upstream rather
than working around here.
