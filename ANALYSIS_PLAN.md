# Analysis plan: calibrating plate growth against relatedness to cultured strains

Owen Swanson. Registered 2026-10-22, before any Phase 3 model was fitted. This file is committed
to the repository so that its date and every later change are public. Any change after
registration is listed in section 10 with its reason.

## 1. Question and hypotheses

Estimates of how much of the microbial world can be cultured often rest on a shortcut: if an
environmental 16S sequence closely matches a strain that has been grown in a lab, the organism
is treated as culturable (Martiny 2019, ISME J 13:2125; disputed by Steen et al. 2019, ISME J
13:3126). This analysis tests the shortcut against what actually grew on plates, across
independent studies that sequenced both an environmental sample and plate scrapings from it.

**Primary question.** How well does 16S identity to the nearest cultured relative predict whether
an environmental taxon is enriched on plates, and is that prediction calibrated (correct in
absolute terms) in a study it was not fitted on?

**Hypotheses, stated before fitting.**
- H1 (ranking). Higher identity to a cultured relative means a higher probability of enrichment,
  in nearly every study, when the model is fitted on the other studies.
- H2 (calibration). Absolute probabilities transfer between studies only approximately: I expect
  good ranking but study-to-study differences in level, because media, depth and handling differ.
  I will report whichever way this comes out.
- H3 (database shortcut). The share of environmental reads with a relative at or above 94.5%
  identity overstates the share of reads from taxa that were actually enriched, in most studies.
- H4 (added value). Identity improves held-out prediction beyond source abundance, conditions and
  phylum.

## 2. Data and inclusion

- **Studies.** The 11 study datasets (10 independent study groups) processed into the project
  database in Phase 2: soil (4), activated sludge (3), marine sponge, oak, broiler carcass and
  rumen. The two Zheng et al. 2024 datasets (soil and sludge, one BioProject) form one group. The
  seawater study is held out because its reads were trimmed by the provider so that identity
  cannot be measured over the same window as in the other studies.
- **Reads to ASVs.** One pipeline for every study (primer removal to the study's amplified window,
  merging, expected-error filter, UNOISE3 denoising, mapping at 97%), as committed in `scripts/16`.
  For the Zheng datasets, whose reads were partly deposited without mates, counts come from every
  read (the paired-read counts are a sensitivity analysis).
- **Usable plating conditions.** A plating condition (study, source sample, medium, oxygen,
  temperature, incubation time and any study-specific detail) is used if its pooled plate reads
  and its source sample's pooled environmental reads both reach 1,000 prokaryotic reads.
  Chloroplast, mitochondrial and non-prokaryotic ASVs are removed before anything is counted.
- **No other exclusions.** No study, condition or taxon is removed for any reason not written here.

## 3. Definitions

- **Unit of analysis.** One row per (ASV, usable plating condition), restricted to ASVs present in
  the condition's source sample (at least 2 environmental reads pooled over that sample's
  environmental libraries). Taxa found on plates but absent from the source sample are counted and
  described separately and are not modelled.
- **Outcome (primary): enriched on plates.** At least 5 reads in the condition's pooled plates AND
  plate relative abundance at least 2 times the source sample's relative abundance. Rationale:
  DNA carried over from the inoculum should track its source abundance; growth should raise it.
- **Outcome (secondary): detected on plates.** At least 2 reads in the condition's pooled plates.
- **Predictor: nearest-relative identity.** Global alignment identity, gaps counted (VSEARCH
  `--iddef 1`), between the ASV and its best hit among cultured-strain 16S sequences cut in silico
  to the same primer pair (BacDive strains, RefSeq type material and 16S from GTDB isolate genomes,
  after removing genomes with foreign 16S). No hit at 75% or above is "no hit".
- **Database shortcut.** An ASV "has a cultured relative" if its identity is at least 94.5%.
- **Recovery shares.** Within a source sample, the share of environmental reads (relative
  abundance) belonging to ASVs that (a) have a cultured relative (shortcut), (b) were enriched in
  at least one of that sample's conditions (observed), or, per condition, (c) the model's held-out
  predicted probability of enrichment (predicted).

## 4. Primary model

Mixed-effects logistic regression (R, `glmmTMB`, binomial, logit link):

```
enriched ~ ns(identity; boundary 75, 100; interior knots 90, 94.5, 97, 99) + no_hit
         + log10(source relative abundance) + log10(plate pool reads)
         + medium class + temperature band + oxygen + log2(incubation days) + phylum
         + (1 | study) + (1 | source sample) + (1 | plating condition)
```

Why each term is there:
- **Identity spline and no-hit indicator.** The quantity being calibrated. A natural cubic spline
  lets the curve bend without forcing a shape; the knots sit at the identity bin edges used since
  the pilot, fixed in advance, so the basis does not depend on the data. No-hit ASVs are set to 75%
  and given their own indicator, so "no relative at all" is not forced onto the curve.
- **Source abundance.** Abundant taxa are inoculated in larger numbers, and the fold rule behaves
  differently for rare and common taxa. Without it, identity could stand in for abundance.
- **Plate pool reads.** Deeper plate pools make 5 reads easier to reach.
- **Medium class, temperature band, oxygen, incubation time.** The plating conditions the question
  asks about. Medium class is computed from recipe facts by one ordered rule (selective, habitat
  extract, defined, dilute complex at most 4 g/L organic matter, rich complex). Temperature bands:
  cold (20 C or below), moderate, warm (33 C or above).
- **Phylum.** The eight phyla with at least 1,000 ASVs in the database are separate levels; all
  others and unassigned are "other". This tests whether relatedness works beyond broad lineage.
- **Random intercepts.** Rows from the same study, source sample and plate pool share conditions
  the fixed terms do not capture (handling, crowding, how many colonies grew). Ignoring this
  clustering would make the model overconfident.

**Missing or unseen condition levels.** Factors use sum-to-zero (deviation) coding over their known
levels. A condition whose level is unknown (oak temperature and oxygen, the Hinsu plates pooled
over 27 C and 37 C, media whose recipe is unresolved or not yet found) is coded as zero on every
contrast, which means "the average of the known levels". The same rule applies to a level that
does not occur in the training data of a validation fold (for example microaerophilic plates,
which occur only in the Zheng group). This keeps every condition in the analysis without
inventing a value.

**Predictions for a study the model has not seen** average over the study, source-sample and
plate-pool random effects (Gauss-Hermite quadrature over a normal with the summed variance). A new
study's own effect is unknown, so predicting at a random effect of zero would describe a "median
study" rather than the expected rate.

**Companion models, all evaluated the same way:**
- *Binned rates:* the enrichment rate per identity bin (none, below 90, 90-94.5, 94.5-97, 97-99, 99
  and above) pooled over the training groups, smoothed as (k + 0.5) / (n + 1). The simplest
  transparent calibration.
- *Identity only:* the primary model with only the identity terms and the random intercepts.
- *No identity:* the primary model without the identity terms (for H4).
- *Gradient boosting (LightGBM)* with the same features and fixed settings (`config/phase3.yaml`),
  without tuning, to check whether the mixed model misses structure. The mixed model stays primary.

## 5. Validation and metrics

- **Leave-one-study-group-out.** Each of the 10 groups is held out in turn; every model is refitted
  on the other nine, and the held-out group is predicted. No performance number comes from data a
  model was fitted on.
- **Primary metric: held-out calibration per study group.** Observed-to-expected ratio (O/E, total
  observed enrichments over the sum of predicted probabilities) and calibration slope (the slope of
  a logistic regression of the outcome on the logit of the prediction; 1 is ideal, below 1 means
  predictions are too extreme). Also calibration curves: mean predicted vs observed rate per
  identity bin and per decile of predicted probability.
- **Discrimination: AUC** per held-out group, summarized as the unweighted mean over groups (each
  study counts once, however many rows it has); the pooled AUC is also reported.
- **Brier score** per held-out group.
- **Uncertainty.** Summary metrics: percentile intervals from 1,000 bootstrap resamples of whole
  study groups (held-out predictions resampled, models not refitted). With 10 groups these
  intervals are approximate and I will say so. Curve: 95% intervals from 1,000 draws of the fixed
  effects from their estimated covariance, plus the spread expected for a new study (10th to 90th
  percentile of the study random effect). Recovery shares: bootstrap over source samples within a
  study; single-source studies get no interval, and this is reported.

## 6. Outputs

- The identity-growth curve (marginal, at average conditions) with its uncertainty, overall and
  for environmental vs host-associated studies.
- Per-group held-out O/E, calibration slope, AUC and Brier for every model.
- Fixed-effect estimates (odds ratios with 95% intervals) for source abundance and conditions, and
  how they shift the curve. Most condition contrasts vary mainly between studies; the study random
  effect absorbs between-study differences, so these effects are estimated mostly from studies
  that varied conditions internally. I will report which studies inform each contrast.
- Per study: shortcut share vs observed share vs predicted share of environmental reads, and the
  overestimate factor (shortcut / observed). Observed share uses "enriched in at least one
  condition", which is the most generous reading for the shortcut; the per-condition median is
  also reported.
- Results are reported separately for environmental (soil, sludge) and host-associated (sponge,
  oak, broiler, rumen) studies, so that neither set dominates.

## 7. Sensitivity analyses

Each reruns the full leave-one-group-out evaluation of the primary model (and the binned rates)
with one change, and reports the same metrics and the curve:

1. Outcome: detected (2 reads) instead of enriched.
2. Enrichment read threshold 2 and 10 (fold 2).
3. Fold threshold 1.5 and 5 (5 reads).
4. Identity ignoring terminal gaps (`--iddef 2`), and gap-blind identity.
5. Nearest relative from BacDive culture-collection strains only, and from RefSeq type material
   only (does a genome-only relative predict growth as well as a strain kept in a collection?).
6. Minimum pool depth 5,000 and 7,500 reads (7,500 is the strictest round level at which every
   study keeps all its environmental samples).
7. Zheng datasets with paired-read counts instead of every-read counts.
8. The long V1-V3 sludge study (Verhoeven et al. 2023) reprocessed with a looser expected-error
   filter (maxEE 2), because the strict filter removed more reads from some taxa than others.
9. Relative abundances corrected for 16S copy number (rrnDB).
10. Rows weighted so that each study group carries equal total weight.
11. Without the three studies graded weakest in the census (one source sample each: the two
    single-source sludge studies and the rumen study).

Fits run in this order: primary and companion models first, then the sensitivity analyses in the
order listed. A sensitivity analysis that cannot be completed is reported as not run, with the
reason.

## 8. Inference criteria

Set now, applied as written:
- **H1 supported** if the identity-only model's held-out AUC exceeds 0.5 in at least 8 of 10
  groups and the lower bound of the bootstrap interval for the mean AUC is above 0.5.
- **H2: "calibrated across studies"** if, for the primary model, both the calibration slope (0.8
  to 1.25) and O/E (0.67 to 1.5) fall in range in at least 7 of 10 groups. Otherwise I will
  describe the relationship as ranking taxa well but needing local recalibration of its level,
  and quantify the between-study spread.
- **H3 supported** if shortcut / observed exceeds 1 in at least 8 of 10 groups and the bootstrap
  interval of the median ratio excludes 1.
- **H4 supported** if the mean held-out AUC of the primary model exceeds that of the no-identity
  model and the bootstrap interval of the difference excludes 0.
- A sensitivity analysis "changes the conclusion" if it moves any hypothesis across its criterion.
  Such cases are reported in the main text, not only in a supplement.
- Null or contrary results are reported as results.

## 9. Known limitations, stated in advance

- Ten study groups is few for estimating between-study variation; intervals on study-level
  quantities are wide.
- Plating conditions in published studies are not a random sample of possible conditions.
  Predictions apply to conditions like those represented.
- Amplicon data are biased toward cultured taxa, and the studies amplify different regions (V4,
  V3-V4, V1-V3, full length cut to V4). Region is confounded with study.
- Several studies contribute few enrichment events (rumen, sponges, oak), so their held-out
  metrics are noisy.
- Plate scrapings measure enrichment of DNA on a plate, not colony counts or isolation success.
- Row counts differ greatly between studies (one desert-soil study has 360 conditions). The
  primary fit is unweighted; the equal-weight sensitivity analysis checks whether large studies
  drive the curve.

## 10. Deviations from this plan

Each deviation is listed with its date, what changed and why. The text of sections 1-9 is left
as registered.

### Deviation 1 (2026-10-22): phylum as a random intercept; curve averaged over the data

**When and what was known.** Decided after the registered primary model had been fitted in all
ten validation folds and on the full data, from the fits' convergence diagnostics and
coefficients. No held-out performance metric (AUC, O/E, calibration slope, Brier) had been
computed or viewed.

**What happened.** Two of the eight phylum levels have almost no enrichment events:
Planctomycetota 0 of 32,023 rows in every study, and Chloroflexota 2 of 27,412, both in the
Zheng group. A fixed effect for a level with no events has no finite maximum-likelihood estimate
(complete separation): the estimate drifts towards minus infinity. In the fold that holds out
the Zheng group, both levels reached about -23 on the log-odds scale with standard errors near
300,000. Four of the ten folds reported convergence problems (non-positive-definite Hessian or
singular convergence).

**Why it matters.**
- With sum-to-zero coding, the "average phylum" includes these runaway levels. Every other phylum
  coefficient and the intercept shift to compensate (intercept -17.6 in that fold), so the
  coefficient table cannot be interpreted.
- The registered curve, drawn "at average conditions" with all contrasts at zero, inherits the
  same artefact.
- Held-out predictions are much less affected (rows of those phyla are predicted near zero, which
  is also what was observed), but intervals that depend on the covariance matrix are unreliable.

**Change 1: phylum enters as a random intercept, not as fixed levels.**
- The same nine levels (the eight phyla with at least 1,000 ASVs, and "other") are kept.
- Their effects are drawn from a common normal distribution whose spread is estimated. A phylum
  with no events is pulled towards the others by an amount set by the data, so its estimate stays
  finite and very low, rather than being removed or left unbounded.
- Phylum is always known for a held-out row, so held-out predictions include the estimated
  phylum effect. The average over unknown random effects for a new study still covers study,
  source sample and plate pool only.
- The identity-only model keeps no phylum term, as registered.

**Change 2: the curve is averaged over the observed rows.**
- At each identity value, every analysed row is given that identity (with a hit) and the
  new-study predictions are averaged over a fixed random subsample of rows (seeded).
- This gives the expected enrichment rate at that identity for the mix of source abundances,
  depths, conditions and phyla actually studied. It is the standard way to draw a covariate-
  adjusted curve, and it does not depend on how factor levels are coded.
- The habitat curves are averaged over the rows of each habitat group (environmental vs
  host-associated).

**Alternatives considered.**
- Merging the phyla with almost no events into "other": rejected because it chooses levels from
  outcomes.
- A weakly informative prior or bias-reduced (Firth-type) estimation on all fixed effects: valid,
  but it changes every coefficient a little, where the problem sits in one term.
- Dropping phylum: loses a registered question, whether relatedness works beyond broad lineage.

**Reporting.** Results from the model exactly as registered (held-out predictions and metrics)
are kept and reported alongside, so readers can see the change does not drive the conclusions.
All sensitivity analyses use the revised model. The pass/fail criteria in section 8 are unchanged.

### Deviation 2 (2026-10-23): group-count criteria scaled to the number of groups analysed

**When and what was known.** Decided after all sensitivity analyses had been evaluated, so the
results were known. It changes no verdict in the primary analysis or in any variant with all ten
groups.

**What happened.** Section 8 states three criteria as counts out of ten validation groups (H1 and
H3: at least 8; H2: at least 7). The analysis that drops the three single-source studies has
seven groups, and the plan did not say how a count out of ten applies to fewer groups. Read
literally, H3 fails there even though all seven groups pass, because eight cannot be reached.

**Change.** Each count is applied as the same share of the groups analysed, rounded up: 8 of 10
becomes 6 of 7, and 7 of 10 becomes 5 of 7. With ten groups the criteria are exactly as
registered. The interval parts of H1 and H3 are unchanged.

**Effect.** In the seven-group analysis H3 is supported (7 of 7 groups). H1 is still not
supported: 5 of 7 groups have AUC above 0.5, one short of 6, although the lower bound of the
mean AUC stays above 0.5. Both the literal and the scaled readings are reported.

**Alternatives considered.**
- Applying the counts literally: makes any analysis with fewer than eight groups fail H1 and H3
  by construction, whatever the data show.
- Calling the count criteria not applicable for such analyses and describing them in words:
  leaves the analysis without a verdict, where the intent of the rule (a large majority of
  studies) can still be applied.

### Exploratory analysis 1 (2026-10-23, not registered): local recalibration

**Status.** Added after H2 was found unsupported. It is not a deviation: no registered model,
metric or verdict changes, and it carries no pass/fail verdict of its own. It is reported as
exploratory.

**Question.** H2 failed because a new study's overall growth rate is not predictable from the
training studies. If a small local pilot measured that rate, would the predictions for the rest
of the study fall inside the H2 ranges?

**Method.**
- Input: the held-out predictions of the primary model (each study predicted by a model that
  never saw it).
- Within each study group, a random share of plating conditions (10%, 25% or 50%; at least one
  condition, and at least one left over) is the pilot. The rest is scored.
- Intercept recalibration: one shift on the log-odds scale chosen so the pilot's predicted
  events equal its observed events. A pilot with no events is treated as half an event.
- Logistic recalibration: regress the pilot's outcomes on the logit of the predictions (an
  intercept and a slope), when the pilot has at least 10 events.
- Score O/E and calibration slope on the remaining conditions only, with the H2 ranges, and
  repeat over 100 random splits. A group counts as calibrated when the median O/E and median
  slope over its splits are both in range.

**Limits stated in advance.** An intercept shift cannot change the calibration slope, so it can
help only the groups whose slope is already in range. Pilot and scored conditions share source
samples and taxa, so this measures what a pilot plating of the same sample would give, not a
different sample from the same site.

## 11. Phase 4: medium classes, the curation hypothesis and under-recovered taxa

Registered 2026-10-23, after the Phase 3 results were final and before any Phase 4 model was
fitted. Known at registration: every Phase 3 result (including the overall medium-class and
phylum effects in the primary model); for the Phase 4 data, counts of rows, taxa, plating
conditions and linked culture-collection records. No enrichment rate by taxon and medium class,
and no association between culture-collection records and enrichment, had been computed.
Sections 1-10 are unchanged. Settings are in `config/phase4.yaml`.

**Questions.**
- Q5 (medium class by taxon). Which medium classes favour which families and phyla, beyond
  each taxon's overall growth rate and beyond what differs between studies?
- Q6 (curation hypothesis). Culture collections record the media each strain is kept on. Does
  knowing that a taxon's nearest cultured relative is kept on medium class X predict that the
  taxon is enriched on plates of class X specifically, once relatedness is accounted for?
- Q7 (under-recovered taxa). Which taxa are enriched less often (or more often) than
  relatedness, source abundance and plating conditions predict, consistently across studies?

**Hypotheses, stated before fitting.**
- H5 (taxon-specific medium effects). (a) Families differ in which medium classes favour them:
  the spread of family-by-medium-class effects is above zero. (b) These preferences transfer:
  estimated in other studies, they improve prediction in a study held out.
- H6 (curation). The claim tested is that a relative kept on medium class X predicts enrichment
  on X specifically. From the soil pilot I expect it not to hold: collection records say where a
  strain is maintained, not where it can grow from a mixed community. I will report whichever
  way it comes out.
- Q7 is descriptive: a ranked list with a rule, fixed here, for calling a taxon systematically
  under- or over-recovered. It has no overall verdict.

## 12. Phase 4 data

- **Rows.** The Phase 3 primary rows (section 3 definitions and outcome unchanged).
- **Medium classes.** Recomputed from `config/media.yaml` by the section 4 rule. Four media of
  Hinsu et al. 2021, unknown in Phase 3, now have recipes from the supplier's technical data
  sheets: actinomycete isolation agar (rich complex), Azotobacter agar with mannitol (habitat
  extract: it contains soil extract), yeast extract mannitol agar and yeast mannitol agar with
  Congo red (both rich complex). Phase 3 keeps them unknown, as registered. Two media of the
  broiler study (Yu et al. 2020) still lack a recipe; that study used only rich media, so it
  informs no medium contrast either way.
- **Taxa.** SILVA family (SINTAX, cutoff 0.8). An ASV without a family gets "unassigned family
  (<phylum>)"; these groups are modelled but not interpreted. Phylum: the full SILVA phylum.
- **Q5 and Q6 studies.** The study groups with at least two known medium classes among their
  usable conditions: Hinsu et al. 2021 (rich, dilute, habitat extract), Indraningrat et al.
  2019 (dilute, defined, selective), Li et al. 2023 (rich, dilute, defined, habitat extract),
  Villanelo et al. 2026 (selective, habitat extract) and Zheng et al. 2024 (dilute, habitat
  extract). Only their rows with a known class (rich complex, dilute,
  defined, habitat extract, selective) are used: 424,879 rows, 442 families (including the
  unassigned groups), 68 phyla. Reason: a medium contrast is
  informative only inside a study, where the sample, handling and sequencing are shared.
- **Q7 studies.** All ten study groups and all Phase 3 rows.

## 13. Q5: medium class by taxon

**Model** (glmmTMB, binomial, logit): the Phase 3 primary fixed effects (identity spline, no-hit
indicator, log source abundance, log plate depth, log incubation days, medium class, temperature
band, oxygen; same coding) and random intercepts for

```
study + source sample + plating condition + family + family:study + family:medium class
```

- *family:medium class* is the quantity of interest: how much better or worse a family does on
  one class than its own average and than the class average. As a random effect it is shrunk
  towards zero when a family has little data on a class, which avoids the separation problem of
  Deviation 1 (1,457 family-by-class cells, many with few events).
- *family:study* stops a family that happens to grow well in one study, which used only some
  classes, from looking like a medium preference.
- Relatedness and conditions are adjusted as in Phase 3, so a family's class preference is not
  just its relatives' culture-collection coverage.
- The same model with phylum in place of family is fitted as a coarser companion.

**Inference.**
- H5a: likelihood-ratio test of the model against the same model without family:medium class.
  The variance is zero under the null, on the edge of its range, so the p-value is half the
  chi-square(1) tail. Supported if p < 0.05. The estimated SD is reported.
- H5b: leave-one-study-group-out over the Q5 groups. Both models (with and without
  family:medium class) are refitted without the held-out group and predict it. The family and
  family:medium class effects are known for a held-out row when the family (and the pair) was
  seen in training; study, source, condition and family:study effects are averaged over, as in
  Phase 3. Metric: mean held-out log loss per group (and AUC). Supported if log loss is lower
  with family:medium class in at least 4 of the 5 groups and the bootstrap interval (groups
  resampled, 1,000 times) of the mean difference excludes 0. With five groups this interval is
  crude.
- Reported: family-by-class effects (conditional modes with 95% intervals from their
  conditional SDs) for families with at least 2,000 Q5 rows in at least two classes (row counts,
  not outcomes); the same for phyla. With hundreds of cells, some intervals will exclude zero by
  chance; this table is descriptive.

## 14. Q6: the curation hypothesis

**Collection records.**
- Culture-medium entries of every BacDive strain record (cached full crawl). Only entries that
  link a MediaDive medium are used, because only those have a structured recipe; growth must
  not be recorded as "no".
- Each MediaDive medium is classed by the section 4 rule, with the recipe facts computed from
  its MediaDive ingredients by fixed lists (`config/phase4.yaml`): antibacterial agents
  (antibiotics and bile; antifungals excluded); habitat material (soil, soil extract, natural
  seawater, rumen fluid, manure and dung, sludge); complex nutrient source (MediaDive's
  complex-ingredient flag, minus gelling agents, single polysaccharides, lipids, salts,
  charcoal, enzymes, detergents and defined media sold as one ingredient); organic g/L (every
  ingredient MediaDive flags complex, so soluble starch and glycerol count as for the plate
  media, and every carbon-containing ingredient, minus gelling agents, insoluble polymers,
  salts, inorganic carbon, buffers, chelators, dyes and detergents). Of 3,342 MediaDive media,
  1,434 are rich, 487 dilute, 718 defined, 142 habitat extract, 15 selective, 242 unresolved and
  304 have no recipe. A complete dehydrated agar medium listed as one ingredient counts its
  weight minus 15 g/L of agar (this moves 8 media, none with a strain record). A complex
  medium with an organic
  ingredient of unstated amount is rich if its stated organic part already exceeds 4 g/L and
  unresolved otherwise.
- **Link from ASV to strains.** Every tied best hit of the ASV (primary identity, all-sources
  reference) is followed to BacDive strains: a strain's own 16S links to that strain; a GTDB
  isolate genome links to every BacDive strain whose record lists that assembly. RefSeq
  type-material 16S records carry no strain identifier here and are not linked.
- **Recorded shares.** For each ASV, rec_X is the share of its linked strains (those with at
  least one classed growth record) that have a growth record on class X. ASVs with no such
  strain are not in Q6. Q6 has 223,800 rows from 10,918 ASVs (50.1% of the Q5 ASVs). Their
  linked strains are recorded on rich media for 45.8% of ASVs, dilute 55.4%, defined 18.0%,
  habitat extract 1.5% and selective 0.1% (so no selective term is fitted).

**Model:** the Phase 3 primary model (fixed effects and random intercepts for study, source
sample, plating condition and phylum) plus

```
rec_X for each class X recorded for at least 1% of Q6 ASVs + match
```

where *match* is rec for the class of the row's own plating medium. The rec_X terms absorb
"taxa whose relatives are kept on X grow more (or less) everywhere"; *match* measures what is
specific to the medium actually used.

**Inference.** H6 supported if both:
- the odds ratio for *match* (all linked strains recorded on the plate's class vs none) has a
  95% interval above 1 in the full fit; and
- leave-one-study-group-out over the Q6 groups: mean held-out log loss is lower with *match*
  than without in at least 4 of the 5 groups.
Otherwise H6 is not supported, and the odds ratio's upper bound is reported as the largest
effect compatible with the data.

**Companion: matched plates.** Pairs of usable plating conditions from the same study and source
sample, with the same oxygen, temperature band and incubation time but different medium classes.
For an ASV whose linked strains are recorded on exactly one of the two classes (any strain) and
which was enriched on exactly one of the two plates, the pair falls in one cell of a 2 x 2 table
(recorded on class x or y; grew on x or y). Tables are pooled over strata (study, class pair) by
the Mantel-Haenszel odds ratio, with a 95% interval from 1,000 bootstrap resamples of ASVs within
study groups. An odds ratio above 1 means taxa grew on the class their relatives are kept on more
often than on the other. Reported with the model result; it has no verdict of its own.

**Sensitivity analyses** (full fit with and without *match*; odds ratio and likelihood-ratio
test reported): (i) rows with identity at least 97%, where a relative's records are most
relevant; (ii) indicators (any linked strain recorded on X) instead of shares; (iii) links
through BacDive strains' own 16S only.

## 15. Q7: under-recovered taxa

- **Model:** the Phase 3 primary model without the phylum term (identity, source abundance,
  depth and conditions; random study, source sample and plating condition), fitted leave-one-
  study-group-out over all ten groups with the Phase 3 rows and coding.
- **Level.** Phase 3 showed that a new study's overall rate is not predictable. Within each
  held-out group, the new-study predictions are shifted by one constant on the log-odds scale so
  that predicted events equal observed events. A taxon's ratio then compares it with the rest of
  its own study.
- **Per taxon** (family; also phylum): O = observed enriched rows, E = sum of shifted
  predictions, pooled over groups as O/E, with a 95% interval from 1,000 bootstrap resamples of
  study groups.
- **Evaluated:** taxa with E of at least 5 in at least 3 groups (E is a model quantity, not an
  outcome); only those groups count for the taxon.
- **Call.** *Systematically under-recovered:* pooled O/E below 1, interval upper bound below 1,
  and O/E below 1 in at least two thirds of the taxon's groups. *Over-recovered:* the mirror
  image. Otherwise neither.
- **Sensitivity:** the identity-only model's Phase 3 held-out predictions in place of the Q7
  model; no level shift.
- Many taxa are tested and the intervals are not corrected for multiplicity; the consistency
  rule guards against a call driven by one study. The list is descriptive.

## 16. Phase 4 outputs

- Q5: likelihood-ratio test and SD; held-out log loss and AUC per group, both models; a
  family-by-class table and heat map; the phylum version.
- Q6: odds ratio for *match* with interval, likelihood-ratio test, held-out log loss per group,
  the matched-plate odds ratio, the three sensitivity analyses; coverage of collection records
  (share of ASVs linked, classes recorded).
- Q7: per-taxon O/E with interval, groups and call; the sensitivity analyses.
- Environmental and host-associated studies are listed separately in every per-group table.

## 17. Phase 4 limitations, stated in advance

- Five study groups inform Q5 and Q6. Transfer tests and bootstrap intervals over five groups
  are crude.
- Medium classes are broad; media within a class differ.
- Culture-collection records are dominated by DSMZ media and by rich media; a class with few
  records gives the test little power on that class.
- At the amplified region, tied best hits often span several species, so an ASV's "relative" can
  be a set of strains. The shares average over them.
- The Q7 level shift uses each group's own outcomes (one number per group), so Q7 compares taxa
  within a study and says nothing about a study's overall level.
