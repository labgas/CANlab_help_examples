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

**Status (2026-10-05): the default was flipped, but to a THIRD engine.** The plan
above anticipated `'predictive_model'` taking over from `'legacy'`. What happened
instead is that both turned out to tune by the same round-robin over row index
and to score identically (+0.0158 on proj_discoverie model_2k immune_PC1), while
`'tuned_nested'` — added later, and not part of this port — scored +0.1812 on the
same data, folds and mask. So `'tuned_nested'` is now the default in `prep_3a`,
`c2a` and `a2`, as part of a coupled recommended configuration documented in
`prep_3a`'s header.

**The legacy branch was NOT deleted, and should not be**: it is how existing
results stay reproducible, and `c2a` still needs it to bootstrap a legacy fit.
What did change in `c2a` is that `'tuned_nested'` now reaches the
`@predictive_model` bootstrap branch instead of being rejected. It needs nothing
new: `prep_3a`'s tuned_nested branch already attaches its full-data refit — at
the modal tuned hyperparameter — as `mvpa_stats.pm`, so that branch bootstraps
exactly the licensed model without refitting. This is the route proj_discoverie
model_2k already used, by pinning `mvpa_engine = 'predictive_model'` by hand in
`s7c` while reading `s6c1t`'s results; the dispatch change removes the manual
step. The legacy `predict()` path is **not** an alternative here — it refits at
`predict()`'s default shrinkage, which per the measurement at the top of this
file reduces to unregularised PCR, so a tuned result that somehow lacks `.pm`
now errors rather than being quietly rebuilt with `'estimateparam'`.

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

## Three engines, and why (added 2026-09-24)

The port now carries three implementations of the MVPA-on-covariate fit. They
differ in ONE thing that matters: how the inner CV for hyperparameter tuning is
constructed.

| engine | file | inner CV | structure-aware? |
|---|---|---|---|
| legacy `predict()` | (in prep_3a) | reuses the OUTER partition | yes, but inner k welded to outer k |
| `@predictive_model`, as shipped | `mvpa_reg_cov_predictive_model.m` | `estimateparam` -> round-robin over ROW INDEX | **no** |
| `@predictive_model`, tutorial pattern | `mvpa_reg_cov_tuned_nested.m` | rebuilt from the training subset's strata | yes |
| `ooFmriDataObjML` | `mvpa_reg_cov_oofmri.m` | partitioner HANDLE re-derived per level | yes, by construction |

### The defect in the shipped class path

`fit_lassopcr` tunes lambda with `cv_assignment = mod(0:n-1,5)'+1` - a
deterministic round-robin over row index. Under exchangeable rows that is a
valid partition. With GROUPED data it splits a dependent cluster across inner
folds and leaks; with STRATIFIED outer folds it selects lambda under a different
sampling model than the one being estimated. The function accepts a
`cv_assignment` argument that would fix this, but nothing supplies it - `fit.m`
calls it with three arguments. Note the docstring claims the implementation is
"faithful to the legacy fmri_data.predict cv_lassopcr"; on this point it is not.

CanlabCore's own tutorials do NOT use `estimateparam`. Part 3 builds the inner
splitter explicitly (`cv_splitter.stratified_group_kfold(4)` with `groups`
sliced to the training rows) and recommends tuning lasso-PCR via `lasso_num`
through `grid_search`. `mvpa_reg_cov_tuned_nested.m` implements that pattern.

### Two gaps the tutorial pattern does NOT close

1. **No composable tuned-estimator object.** `grid_search` returns a model with
   fixed `modeloptions`; there is no object representing "estimator + its tuning
   procedure". Because `permutation_test` and `bootstrap` call `crossval`
   internally, they cannot wrap a tuned model - a permutation test of a tuned
   model is not expressible in the API, which is why
   `mvpa_reg_cov_tuned_nested.m` has to permute the whole nested procedure
   itself.
2. **`select_features` is not refit per fold.** `crossval` clones and refits per
   fold, and `fit` standardises internally, so SCALING is not leaked. But
   `select_features` is applied to the full data and carried via
   `omitted_features`, so calling it before `crossval` leaks the outcome into
   feature selection. Neither engine here uses it.

`ooFmriDataObjML` closes both: `gridSearchCV(est, grid, innercv)` returns an
estimator, so `crossValScore(gs, outercv, scorer)` nests by composition, and its
`pipeline` refits every transformer per fold.

### Upstream proposal

Give `@predictive_model` a tuned-estimator wrapper holding (base model, inner
splitter, grid) whose `fit()` runs the inner search. That is `bayesOptCV`'s
design, it makes `crossval` / `permutation_test` / `bootstrap` compose over
tuned models for free, and it subsumes the narrower point that `cv_assignment`
should accept a SPLITTER rather than a vector - a vector is a partition of one
particular dataset and cannot re-derive itself for a subset.

### Status

All four have now been run on `proj_discoverie` model_2k, immune_PC1, 5 outer
folds stratified on centre. `ooFmriDataObjML` is UNMAINTAINED (last commit
2024-08-23) but is already a fork dependency via `prep_3c` and `c2f`; treat it
as a reference implementation and a cross-check, not as the default.

---

## Three findings from running them (2026-09-25)

### 1. `get_r` is a LOSS, and the ooFmri engine misread it

ooFmriDataObjML's scorers are losses to be MINIMISED — `get_r` returns
`1 - corr(yfit, Y)`, which is correct for `gridSearchCV`. `mvpa_reg_cov_oofmri.m`
reported `mean(cvGS.scores)` as r, giving an impossible "+0.9148" with per-fold
values of 1.1888 and 1.0661. Fixed: the engine now reports `1 - loss`.
**Nothing in ooFmriDataObjML was changed** — the bug was entirely in the wrapper.

### 2. Averaging per-fold r is the wrong summary

`crossValScore` AVERAGES PER-FOLD SCORES; the legacy and `@predictive_model`
engines POOL all held-out predictions and correlate once. These are different
statistics, and the difference is not cosmetic.

Correlation is **not an additive loss**. For MSE — or R² against a global mean —
"mean of per-fold values weighted by n_k" and "computed once over the pool" are
algebraically the same number, so no choice arises. For r they differ, because r
re-centres and re-scales within whatever set it is computed on.

Prefer the **pooled** estimate:

* r is a biased estimator of rho (bias ~ `-rho(1-rho^2)/2n`) and skewed, so
  averaging raw r averages the bias in rather than cancelling it.
* Variance. `SE(r) ~ (1-r^2)/sqrt(n-3)`. At n = 93 with k = 5, per-fold n ~ 18.6
  gives SE ~ 0.25 against ~0.105 pooled. The observed per-fold spread
  (-0.189 to +0.389) is exactly what that predicts, and carries almost no
  information.
* Per-fold r is undefined when a fold's Y has little variance — which stratified
  partitions make MORE likely, not less.

Where an average is genuinely wanted, average in **Fisher z**:
`tanh( sum((n_k-3) * atanh(r_k)) / sum(n_k-3) )`. z is near-normal with variance
~`1/(n-3)` INDEPENDENT OF RHO, which is what makes `(n_k - 3)` the right weights.
The engine now reports all three.

**The caveat on pooled r, which is real.** It mixes predictions from k different
models, so a dataset where every fold has r ~ 0 internally, but fold means of
yfit track fold means of Y, yields a large and entirely spurious pooled r —
between-fold variance read as prediction. Per-fold r is immune to this; pooled r
is not. The precondition is that fold membership must not predict Y, and it
should be CHECKED, not assumed. On model_2k: `F(4,88) = 0.579, p = .68,
eta^2 = .026`, so the pooled figures are clean. Note the outer partition is
stratified on CENTRE, not on Y, so this was not guaranteed by construction.

Pooled points are also NOT independent (overlapping training sets), so the
ordinary r-based p-value and CI are wrong. Inference comes from the permutation
test.

### 3. `prep_3a` did not mask the MVPA features (fixed)

`mvpa_data_objects{covar} = cat_obj` took the UNMASKED contrast object.
`prep_3a` applies `glmmask` only to the statistic image after the univariate fit
— correct and sufficient there, because each voxel's test is independent, so
restricting the map afterwards leaves every surviving statistic identical.

**That reasoning does not carry to MVPA.** lasso-PCR decomposes over every
included voxel, so out-of-mask voxels shape the components, the weights and the
cross-validated prediction, and no post-hoc mask undoes it. On model_2k:
**235807 voxels unmasked against 149154 in the canlab2023 grey-matter mask — 36.7%
of the features were white matter, CSF and edge.**

Note the TFCE branch already did this correctly
(`cat_obj_tfce = apply_mask(cat_obj_tfce, glmmask_tfce)`), so the template
already encoded the principle; the `mvpa_reg_cov` block was the only
feature-building path that missed it. All ten `apply_mask` call sites were
audited to confirm that.

Fixed behind `domask_mvpa_reg_cov` (default `true`), documented in the a2
template. **Every engine number measured before this fix was computed on the
unmasked object and is superseded.**


---

## The engines measured on masked data (2026-09-25)

`proj_discoverie` model_2k, immune_PC1, n = 93, 5 outer folds stratified on
centre, 149154 grey-matter voxels. Same data, same outer folds, same mask
throughout. The fold reconstruction was verified against the value s6c1 wrote to
disk (+0.1668, reproduced to 1e-6), so "same folds" is measured, not assumed.

| engine | inner CV design | r |
|---|---|---|
| `tuned_nested` | inner folds rebuilt from the TRAINING subset's strata | **+0.1812** |
| `ooFmriDataObjML` | partitioner HANDLE re-derived per level | **+0.1186** |
| legacy `predict()` / `estimateparams` | round-robin over ROW INDEX | +0.0158 |
| `@predictive_model` / `estimateparam` | round-robin over ROW INDEX | +0.0158 |
| *legacy, UNMASKED control* | *round-robin over ROW INDEX* | *+0.1668* |

**The split is by inner-CV design, not by the mask.** Both engines that rebuild
the inner partition under the outer folds' constraints find grey-matter signal;
both engines that tune through the structure-blind round-robin return ~0.016 on
exactly the same data. `tuned_nested` finds MORE signal inside grey matter than
the legacy engine found on the whole unmasked brain.

This is the defect described in "The defect in the shipped class path" above,
now with a measured cost attached. Note the two `estimateparams` engines agree
to four decimals: they are one engine behind two front-ends.

### What was ruled out

The unmasked legacy predictions carry a long negative tail (min -1.774 against
Y sd 0.7745), so the obvious suspicion was that +0.1668 was a few extreme
predictions inflating a Pearson correlation. **It is not:**

| | Pearson | Spearman | drop 1 extreme | drop 3 extreme | \|z\|>3 |
|---|---|---|---|---|---|
| unmasked | +0.1668 | +0.1132 | +0.1628 | +0.1688 | 1 |
| GM-masked | +0.0158 | -0.0187 | +0.0627 | +0.0598 | 3 |

The unmasked association survives in ranks and is unchanged by dropping the
extremes; the masked predictions are in fact MORE outlier-ridden. So the tail is
real but carries none of the association, and the collapse is not an outlier
artefact.

### How much of this is noise

The inner-fold seed alone moves `tuned_nested` by ~0.06 (+0.2198 at seed
20260923 against +0.1609 at 20260925, same data and same outer folds). So the
gap between +0.1186 and +0.1812 is WITHIN seed noise and the two should not be
ranked against each other on these numbers. What is robust is the gap between
~+0.15 and +0.016.

For the same reason a single point estimate should not be reported alone. The
5000-permutation run is structured as 5 blocks of 1000 with different seeds, and
reports the spread of the observed r across blocks alongside the pooled null.

### The mask effect itself, measured on identical folds

| model | legacy (disk) | tuned_nested unmasked | tuned_nested GM-masked | p(40 perm) | null mean |
|---|---|---|---|---|---|
| model_2k immune_PC1 | +0.1668 | +0.1609 | +0.1812 | 0.073 | -0.0636 |
| model_2j comorbidPCA | +0.0199 | +0.0013 | +0.0414 | 0.146 | -0.0326 |

Masking slightly HELPS the properly-nested engine in both models. model_2j is
null either way.

Both null means sit just below zero, which is the expected, healthy result: the
null distribution of a CROSS-VALIDATED r is not centred at zero, because under
the null the model still fits training noise that does not generalise, so
held-out predictions are mildly anti-correlated with Y. The null-centring check
in `mvpa_reg_cov_tuned_nested` was originally `abs(mean) > 0.10`, which would
have fired on every correct run; it now flags a POSITIVE null mean (fold
membership carrying signal) or a strongly negative one.

### Wrapper fixes made along the way

- `mvpa_reg_cov_oofmri` ran the whole nested CV TWICE - `crossValScore` for the
  per-fold losses and `crossValPredict` for the pooled predictions. It now makes
  a single `crossValPredict` pass and derives everything from `.yfit` and
  `.cvpart`, including the true fold sizes for the Fisher weights (previously
  assumed equal). It also no longer depends on the `get_r` sign convention at
  all, since correlations come straight from the held-out predictions.
- Running several MATLAB sessions that each call `saveProfile` on the shared
  `local` cluster profile deadlocks pool startup - one session hung 15 minutes
  at "Starting parallel pool" with the machine idle. Set `NumWorkers` on the
  cluster OBJECT and pass it to `parpool`; never persist it.


---

## Pattern inference wired into c2a (2026-09-26)

`prep_3a` answers *is the model better than chance* (permutation test). Only if
that is significant is it worth asking *which voxels does it lean on*, which is
what bootstrap and stability selection answer. So the pattern inference lives in
**c2a**, not prep_3a: the expensive step is opt-in once the cheap model-level
result is known. c2a's previous comment recommended the opposite and has been
corrected.

**The two engines are exclusive.** `mvpa_engine` selects which bootstrap c2a
runs — `'legacy'` (`fmri_data/predict` with `bootsamples`) or
`'predictive_model'`. Never both, so there is one set of weight maps to report.

### The `'predictive_model'` path

Routes to `mvpa_reg_cov_predictive_model`, which already wraps `bootstrap`,
`stability_selection` and `weight_map_object`. c2a adapts the result back into the
legacy `mvpa_bs_stats{j}.weight_obj` shape so the existing montage and threshold
code runs unchanged, and keeps the full object on `.pm`.

**Folds are recovered, not reinvented.** `prep_3a` saves `teIdx` as a **cell of
`nfolds` logical `[n x 1]` test masks** — not a matrix, so it cannot be collapsed
by multiplication (an earlier draft tried, and would have failed). c2a rebuilds
fold labels from it and asserts every subject lands in exactly one test fold;
otherwise the bootstrap would not match the cross-validated fit it follows up.
Verified on model_2k: `{1x5}` of logical `[93x1]` giving `[18 19 19 19 18]`.

**The bootstrap-p collapse is detected and reported.** On a strongly regularised
model the weights are near-identical across resamples, so the empirical p floors
at `2/(nboot+1)` for every voxel and the FDR mask becomes meaningless — the class
documentation is explicit about this. c2a counts voxels at the floor and, when
that is most of them, says in the report that the z/p and FDR mask are not
interpretable and that stability selection is the inference to read.

### Choosing the stability-selection threshold

Meinshausen & Buhlmann (2010) Thm 1: `E(V) <= k^2 / ((2*pi - 1) * p)`, valid only
for `pi > 0.5`. **k and pi cannot be chosen separately**, and `k` enters squared
while `p` enters linearly — which is what breaks the class defaults at brain
scale. Measured at `p = 149154`:

| k | pi = 0.6 | pi = 0.9 |
|---|---|---|
| **2000** (class default) | **134.1** | **33.5** |
| 500 | 8.4 | 2.1 |
| 345 | 4.0 | **1.0** |

At `k = 2000` no valid `pi` reaches `E(V) <= 1` — the algebra demands
`pi >= 13.9`. So c2a inverts the formula instead:
`k = sqrt((2*pi - 1) * p * E(V))`, giving **k = 345** at `pi = 0.9`, `E(V) = 1`.

| option | default | behaviour |
|---|---|---|
| `stab_threshold_mvpa_reg_cov` | **0.9** | must exceed 0.5, else hard error |
| `stab_EV_mvpa_reg_cov` | **1** | false-selection budget, used to derive k |
| `stab_k_mvpa_reg_cov` | **empty** | empty derives k; an explicit value is used and its implied E(V) reported back |

Either way the implied `E(V)` is printed, so a threshold never appears without its
error characteristics; `E(V) > 10` also warns and names the k that would reach 1.
The montage is **thresholded at pi** (c2a reports thresholded results), with
frequency kept as blob intensity so it still shows how stable each survivor is.

**Caveat for a Methods section:** the theorem assumes SUBSAMPLING at n/2
(complementary pairs); `stability_selection` resamples WITH REPLACEMENT, as does
the `boot_w` reuse below. Shah & Samworth (2013) derive the bound for that case
and it differs. The printed `E(V)` is the right order, not an exact guarantee.

### Speed: stability selection reuses the bootstrap

`stability_selection` refits a model per resample in a **serial** `for` loop.
`bootstrap` does the *same* resampling — the two blocks are line-for-line
identical, `randi(n,[n,1])` or whole-group sampling, then `clone`+`fit` — and
retains every weight vector in `pm.weights.boot_w` as `[p x nboot]`. Stability
selection is then one sort per column of a matrix already in memory.

`mvpa_reg_cov_predictive_model` therefore gains `stab_reuse_boot` (default
`true`): when a bootstrap of at least `nstab` samples has run, stability is
derived from `boot_w` rather than refitting. Same estimator, same scheme, same
fits — seconds instead of `nstab` model fits. Set it `false` to validate against
the method's own path.

**Note on the parallelism.** On this machine `bootstrap.m` runs its loop under
`parfor`, but that is an **uncommitted working-tree edit** in CanlabCore
(`for` -> `parfor`, 2026-07-01; `permutation_test.m` carries the same edit from
2026-07-13; `stability_selection.m` never got it). Upstream both are serial.
Reuse helps either way, and more so upstream. Do not assume the `parfor` is
present — and note that any result produced through these methods currently
depends on a working-tree change the provenance tooling cannot record.
