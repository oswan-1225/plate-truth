---
title: "Having a cultured relative is weak evidence of growth on a plate: a cross-study calibration against plate-grown communities"
author: "Owen Swanson"
date: "{{date}}"
---

<!-- Source of manuscript/preprint.md. Every number in double braces and every citation key in
square brackets is filled in by scripts/47_manuscript.py from results/ and from Crossref
records. Edit this file, not preprint.md. -->

# Abstract

Estimates of how much of the microbial world can be cultured often rest on a shortcut: an environmental 16S rRNA gene sequence that closely matches a cultured strain is counted as culturable. Whether that inference holds has been argued in print but not tested against what actually grows. I assembled {{n_datasets}} public datasets from {{n_groups}} independent studies in which both an environmental sample and the communities grown from it on agar plates were sequenced, spanning soil, activated sludge, a marine sponge, oak, broiler carcasses and rumen. All reads were reprocessed with one pipeline, and for {{n_rows}} combinations of sequence variant and plating condition I asked whether identity to the nearest cultured relative predicts enrichment on plates, with every performance number computed on a study the model was not fitted on. Identity ranked taxa well in new studies, with a mean held-out AUC of {{h1_auc}}, but the overall level of growth differed between studies more than tenfold and could not be predicted (observed-to-expected ratios {{oe_min}} to {{oe_max}}). The shortcut "has a cultured relative at 94.5% identity or above" overstated the share of the community that was enriched on plates in every study, by a median factor of {{h3_ratio}}. Most bacterial families grew less often than relatedness predicts, and culture-collection records of the media strains are kept on carried no medium-specific information once studies were allowed their own baselines. Applied to {{n_survey_studies}} public surveys from seven habitats, the shortcut counted {{shortcut_min}} to {{shortcut_max}} of reads as culturable, above what one plating is predicted to recover even in a high-recovery study in every habitat. A cultured relative is useful for ranking which taxa are likely to grow, but it is not a measure of how much of a community can be grown.

# Introduction

Many bacterial lineages detected by sequencing have no cultivated representative, and in soil, seawater and freshwater many abundant organisms are only remotely related to any strain that has been characterized [@rappe_review]. How large this uncultured majority is matters for how much of microbial physiology can be studied in the laboratory. Two kinds of estimate dominate. One counts taxa or cells whose 16S rRNA gene sequences fall in lineages without any cultured member [@lloyd2018]. The other counts the share of sequences that have a close cultured relative and treats those as culturable [@martiny2019]. The second reading led Martiny to conclude that a median of 97.3% similarity to cultured relatives, and about half of sequences with a relative above 97%, meant the "1% culturable" paradigm was no longer correct [@martiny2019]. Steen and colleagues replied that the conclusion rested on technical errors, notably identity calculations that ignore gaps and amplicon primers that favour cultured organisms, and on a conceptual one: similarity to a cultured strain cannot be used to infer culturability, because physiology cannot be read from a 16S sequence [@steen2019].

Both positions rest on a claim that can be tested directly but, as far as I can find, has not been: that a sequence with a cultured relative will grow when the community is plated. Many studies have sequenced plate-grown communities alongside the source community [@lau2016; @zheng2024], and some report that cultivated taxa are closer to type strains than uncultivated ones while noting that closely related taxa often fail to grow [@zheng2024]. Others predict cultured status from genome content, with cultured status taken from databases rather than plates [@oduwole2025]. What is missing is a calibration: across independent studies, how does the probability that a taxon grows on plates change with its identity to the nearest cultured relative, and is that probability transferable to a study the model has not seen?

I tested this with public datasets in which an environmental sample and the communities grown from it on plates were both sequenced. I reprocessed every dataset with one pipeline, measured identity to cultured strains over each study's own amplified region with gaps counted, and defined growth as enrichment on plates relative to the source community, so that DNA carried over from the inoculum is not counted as growth. Hypotheses, models, metrics and criteria were registered before any model was fitted, and every performance number comes from a study held out of the fit. I then asked which plating conditions and taxa shift the relationship, whether culture-collection records of the media strains are kept on predict the medium a taxon grows on, and what the calibrated model implies for public surveys of habitats in which the shortcut has been used to estimate culturability.

# Results

## Ten studies with paired environmental and plate-grown communities

A search of {{n_articles}} articles and {{n_projects}} sequencing projects identified {{n_census_groups}} independent studies that sequenced an environmental sample and plate-grown communities from it, with raw reads public and the two linkable. One seawater study was held out because its deposited reads had been trimmed so that identity could not be measured over the same window as in the others. The remaining {{n_datasets}} datasets from {{n_groups}} study groups ({{n_runs}} sequencing runs: {{n_env_runs}} environmental libraries and {{n_plate_runs}} plate communities) gave {{n_asvs}} amplicon sequence variants (ASVs) in {{n_conditions}} plating conditions (Fig. 1). After the depth rule, {{n_usable_conditions}} conditions and {{n_rows}} combinations of ASV and condition remained, of which {{n_enriched}} were enriched on plates (at least 5 reads and at least twice the source relative abundance). The studies differ in nearly every respect a plating experiment can differ: from {{min_conditions}} to {{max_conditions}} conditions per study, {{min_sources}} to {{max_sources}} source samples, {{n_medium_classes}} medium classes, incubation from cold to warm and from aerobic to anaerobic (Fig. 1b).

## Identity ranks taxa in a new study, but does not set the level

The probability that an ASV was enriched rose with identity to its nearest cultured relative (Fig. 2a). For a study the model had not seen, the predicted rate was {{curve_90}} at 90% identity, {{curve_945}} at 94.5%, {{curve_97}} at 97%, {{curve_99}} at 99% and {{curve_100}} for an exact match. No ASV without a cultured relative at 75% identity or above was enriched ({{nohit_k}} of {{nohit_n}}). Even at 99% identity, then, most taxa with a cultured relative did not grow on the plates used, and the spread between studies was wide: at 99% identity the 10th to 90th percentile across studies ran from {{curve_99_p10}} to {{curve_99_p90}}.

Held out one study group at a time, identity alone ranked taxa better than chance in {{h1_groups}} of {{n_groups}} groups, with a mean AUC of {{h1_auc}}; H1 was supported (Table 1; Fig. 2b). Identity added to source abundance, plating conditions and phylum: the full model's mean held-out AUC exceeded that of the same model without identity by {{h4_gain}}, supporting H4. A gradient-boosting model with the same features reached a mean AUC of {{lgbm_auc}}, against {{primary_auc}} for the mixed model, so the simpler model does not miss much structure. Ranking was poor in two studies: the broiler-carcass study (AUC {{auc_yu2020}}), where most reads already have a close cultured relative, and the desert-soil study (AUC {{auc_li2023}}).

The level of growth did not transfer (H2 not supported). The ratio of observed to expected enrichments in held-out studies ranged from {{oe_min}} to {{oe_max}}: {{oe_low_n}} groups fell below the registered range (0.67 to 1.5), where the model over-predicted, and {{oe_high_n}} above it, where it under-predicted (Fig. 2c). No group had both a calibration slope and an observed-to-expected ratio inside the registered ranges. The overall growth rate of a new study is the main unknown, not the shape of the curve.

| Hypothesis | Criterion (registered) | Result | Verdict |
|---|---|---|---|
{{table_hypotheses}}

Table 1. Registered hypotheses and held-out results. Intervals are 95% bootstrap intervals over study groups. The medium and curation hypotheses (H5, H6) are reported in their own sections below.

## The shortcut overstates recovery in every study

In each study I compared three read-weighted shares of the source community: reads with a cultured relative at 94.5% or above (the shortcut), reads in ASVs enriched in at least one plating condition (observed), and the held-out prediction (Fig. 3a). The shortcut exceeded the observed share in {{h3_groups}} of {{n_groups}} groups, by a median factor of {{h3_ratio}}, supporting H3 (Fig. 3b). The factor ranged from {{ratio_min}} in the broiler-carcass study to {{ratio_max}} in the sponge study. Using "enriched in any condition" as the observed share is the most generous reading for the shortcut; per condition the gap is larger.

Of {{n_sens}} sensitivity analyses, {{n_sens_same}} left every verdict unchanged, including different read and fold thresholds, identity definitions, cultured references, depth rules, count tables and a copy-number correction. Two changed a verdict. When "grew" is replaced by mere detection on plates (at least 2 reads, no fold rule), the shortcut no longer clearly overstates recovery (median factor {{detected_ratio}}, lower bound {{detected_ratio_lo}}): detection counts DNA carried over from the inoculum, which tracks the source community. Dropping the three single-source studies left {{notierc_groups}} groups, and identity ranked better than chance in {{notierc_auc_groups}} of them, one short of the scaled criterion; the shortcut verdict held there.

## What shifts the curve

In a model that let the curve differ between habitat types, host-associated habitats had a predicted rate of {{host_90}} at 90% identity against {{env_90}} in environmental habitats, and {{host_99}} against {{env_99}} at 99%, with wide, overlapping intervals (Fig. 4a). How sharply identity separated growers from non-growers differed between studies: held-out calibration slopes ranged from {{slope_min}} to {{slope_max}} (Fig. 4b). Above 1, identity separated growers more sharply in that study than the model expected (the rumen and sponge studies); near 0, identity barely ranked (the broiler and desert-soil studies). Studies therefore differ in the strength of the identity effect as well as in level.

Plating conditions shifted the curve much less than the study did (Fig. 4c). Relative to the average medium class, dilute media had an odds ratio of {{or_dilute}} and habitat-extract media {{or_extract}}, against {{or_selective}} for selective media, which contain antibiotics. Warm incubation (33 °C or above) had an odds ratio of {{or_warm}} relative to the average temperature band, and each doubling of incubation time {{or_days}}. Abundant taxa passed the fold rule less often: each tenfold increase in source abundance had an odds ratio of {{or_abundance}}, as expected for a rule that asks a taxon to double its share. Because the study random effect absorbs differences between studies, these condition effects are estimated mainly from studies that varied conditions internally, and the warm band rests on {{warm_studies}} studies.

If the level of a new study is the main unknown, a small local pilot could fix it. In an exploratory analysis (not registered), a share of each held-out study's plating conditions recalibrated the predictions and only the remaining conditions were scored. Shifting the level alone brought the observed-to-expected ratio into range in {{recal_level_oe}} of {{n_groups}} groups even with a 10% pilot, but the slope stayed out of range in {{recal_level_slope_out}}; recalibrating level and slope met the registered calibration bar in {{recal_logistic_25}} of {{n_groups}} groups with a 25% pilot and {{recal_logistic_50}} with a 50% pilot.

## Medium preferences exist within studies but do not transfer

In the {{q5_groups}} studies that plated on two or more medium classes, families differed in which classes favoured them: the standard deviation of family-by-medium effects was {{q5_sd}} on the log-odds scale (likelihood-ratio test, p {{q5_p}}; H5a supported; Fig. 5). The clearest signal was Bacillota on rich media (Bacillaceae {{bacillaceae_rich}}, Planococcaceae {{planococcaceae_rich}} on the log-odds scale); only {{q5_cells_sig}} of {{q5_cells}} reported family-by-class cells had an interval excluding zero, all of them in Bacillota. These preferences did not help predict a held-out study: held-out log loss improved in only {{q5_improved}} of {{q5_groups}} groups (mean change {{q5_dll}}; H5b not supported).

Culture collections record the media each strain is kept on. In the registered test, a taxon whose linked collection strains were all recorded on the plate's own medium class was more likely to be enriched there than one whose strains were recorded on none with an odds ratio of {{q6_or}}, and held-out log loss improved in {{q6_improved}} of {{q5_groups}} groups (H6 supported as registered; Supplementary Fig. S1). An exploratory check, added after this result was seen, allowed each study its own effect of the recorded classes. The match effect then disappeared, with an odds ratio of {{q6_bystudy_or}}, and within single studies it pointed in opposite directions: {{q6_li}} in the desert-soil study and {{q6_zheng}} in the Zheng soils and sludge. The registered result mostly reflects differences between studies: taxa whose relatives are kept on dilute media did well on every plate in the study with the most dilute plates. Whether collection records carry medium-specific information is not established by these data.

## Most families grow less often than relatedness predicts

Compared with a model of identity, source abundance and conditions without taxonomy, shifted to each study's own overall rate, {{q7_under}} of {{q7_families}} evaluated families were systematically under-recovered and {{q7_over}} over-recovered: Pseudomonadaceae, Aeromonadaceae and Rhizobiaceae, with observed / expected ratios of {{oe_pseudomonadaceae}}, {{oe_aeromonadaceae}} and {{oe_rhizobiaceae}} (Supplementary Fig. S2). {{q7_zero}} under-recovered families had no enriched rows at all. At phylum level, {{q7_phyla_under}} of {{q7_phyla}} evaluated phyla were under-recovered, including Acidobacteriota at {{oe_acidobacteriota}}, Verrucomicrobiota at {{oe_verrucomicrobiota}} and Planctomycetota at {{oe_planctomycetota}}. Having a cultured relative is far weaker evidence of growth for an acidobacterium or a planctomycete than for a pseudomonad.

## Public surveys: the shortcut against what one plating could recover

Because a new study's level cannot be predicted, no single recoverable fraction can be given for any habitat. What can be tested without knowing the level is whether the shortcut exceeds even the upper end of the range a plating experiment is predicted to recover. I applied the full model to {{n_survey_studies}} MGnify survey studies with usable samples ({{n_amp_studies}} amplicon studies with {{n_amp_samples}} samples, {{n_shot_studies}} shotgun studies with {{n_shot_samples}} samples), selected in a fixed random order by rules registered before any sample was processed, averaging over the plating conditions of the training studies.

The procedure was first applied to each training study with a model that had not seen it. The predicted upper end of the one-plate range was at or above the observed recovery in {{tc_upper}} of {{n_groups}} studies, so it works as an upper reference. The range itself was too narrow at the low end, containing the observed recovery in only {{tc_range}} of {{n_groups}} studies (registered bar {{tc_range_needed}}), so it is not read as the spread of real studies (Supplementary Fig. S3).

In every habitat and in both survey types, the shortcut exceeded the upper end of what one plating is predicted to recover in a high-recovery study (Fig. 6; Table 2). The difference ranged from {{d_min}} in {{d_min_habitat}} to {{d_max}} in {{d_max_habitat}} in the amplicon surveys, and it was positive in every study. A more generous reference, recovery on any of eight plates assuming the plates are independent, was still below the shortcut on average in every habitat, but the difference was clear only in the human gut, where it was {{d8_gut}}. The registered verdict held under every sensitivity analysis (identity definitions, condition weighting, sample depth, V4 studies only, and without a laboratory-reactor study).

| Habitat | Studies / samples | Shortcut | One plate, 10th-90th percentile study | Shortcut minus one-plate upper end |
|---|---|---|---|---|
{{table_habitats_amplicon}}

Table 2. Amplicon surveys. Shares of reads; intervals are 95% intervals over studies and fixed-effect draws. Habitats marked * have no training study (extrapolation). Shotgun surveys: {{shotgun_sentence}}

The composition of the surveys adds a level-free reason for caution. In soil surveys, {{comp_soil}} of reads belonged to families that were under-recovered in the training studies, and {{comp_shortcut_soil}} of reads counted as culturable by the shortcut sat in such families; in the human gut the shares were {{comp_gut}} and {{comp_shortcut_gut}}. In the marine habitats most reads were in families the training data could not evaluate ({{comp_ne_marine}}), so composition says little there (Supplementary Fig. S4).

The shortcut shares here are much higher than the metagenome-based shares of cells in cultured species to genera reported by Lloyd et al. (for soil {{lloyd_soil}}, against {{shortcut_soil}} here) and similar to Martiny's {{martiny_97}} of sequences with a cultured relative above 97% [@lloyd2018; @martiny2019]. Primer bias towards cultured taxa [@steen2019], a reference set that includes isolate genomes and strain collections, and different definitions of "cultured relative" are likely reasons; these data cannot separate them. The test above does not depend on which estimate is right: whichever is used, it is several times what one plating is predicted to recover.

# Discussion

The shortcut at the centre of the culturability debate makes a testable claim, and against plate data it fails in a specific way. Identity to the nearest cultured relative is real information: it ranks which taxa in a new study are likely to grow, adds to what abundance and conditions predict, and nothing without a cultured relative at 75% identity grew. But a cultured relative is not a measure of how much of a community grows. Most taxa with a close relative were not enriched, the shortcut overstated recovery in every study, and the overall level of growth varied between studies more than tenfold for reasons the recorded conditions did not explain.

This bears on both sides of the debate. Steen and colleagues argued that similarity cannot be used to infer culturability [@steen2019]; the results support that as a statement about levels, while showing that similarity does carry ranking information. Martiny's estimate [@martiny2019] is not contradicted as a count of sequences with cultured relatives, which the survey shares here resemble, but it cannot be read as the share of a community that grows when plated. Gap handling, one of the technical criticisms, made little difference once references were cut to the amplified region.

The finding that most families grow less often than their relatedness predicts, and that recovery concentrates in a few Pseudomonadota families, fits the long-standing observation that standard plates favour fast-growing copiotrophs. It also means that culturability inferred from relatives is least reliable for exactly the lineages, such as Acidobacteriota, Planctomycetota and Verrucomicrobiota, that dominate many soils.

Several limitations apply. Ten study groups are few for estimating between-study variation, and plating conditions in published studies are not a random sample of possible conditions; predictions apply to conditions like those studied. Plate scrapings measure enrichment of DNA on a plate, not colony counts or isolation success. Amplicon regions differ between studies and are confounded with study. Culture-collection coverage is biased towards clinical and industrial taxa. Five of the survey habitats have no training study, the survey studies are those deposited in MGnify rather than a random sample of each habitat, and the predicted range is too narrow at its low end. "Any of K plates" assumes independence between plates and overstates recovery. Finally, the enrichment definition guards against carryover DNA but cannot rule out every artefact; with detection alone the shortcut verdict does not hold, which is itself the reason growth was defined by enrichment.

A calibrated estimate of how much of a habitat can be cultured would need what these data lack: many more studies with shared protocols, or a local pilot plating to fix the level, which in the exploratory analysis brought most studies into calibration. Until then, "has a cultured relative" is best treated as a ranking of candidates for cultivation, not as an estimate of how much of a community can be grown.

# Methods

## Pre-registration

The analysis plan (hypotheses, models, metrics, inference criteria and sensitivity analyses) was committed to the public repository before any model was fitted, and extended before each later stage: the medium and taxon analyses before any of their models was fitted, and the survey application before any survey sample was processed. Three deviations and three exploratory analyses, each dated and justified, are listed in the plan; none changed a registered verdict.

## Study census and inclusion

Studies were found by Europe PMC searches combining plate-scraping and culture-enriched profiling terms with 16S amplicon sequencing ({{n_articles}} unique articles), followed by a check of every linked sequencing project against its ENA run records ({{n_projects}} projects). Inclusion required raw reads in a public archive, environmental and plate samples linkable to the same source sample, medium and temperature reported, and 16S amplicon data. Every study's DOI, first author, year and accession were checked against its Crossref record and the paper's full text by script.

## Read processing

All studies were processed with one pipeline: primer removal and cutting to the study's amplified window with cutadapt [@cutadapt], read merging, an expected-error filter (maxEE 1) and region-specific length windows with VSEARCH [@vsearch], pooled denoising per study with UNOISE3 [@unoise] and de novo chimera removal, and mapping of every run's reads to the ASVs at 97% identity. Taxonomy was assigned with SINTAX [@sintax] against SILVA 138.2 [@silva] with a cutoff of 0.8. Chloroplast, mitochondrial and non-prokaryotic ASVs were removed before anything was counted. For the Zheng et al. datasets, whose reads were partly deposited without mates, counts came from every read; the paired-read counts are a sensitivity analysis.

## Cultured reference and nearest-relative identity

The cultured reference combined 16S sequences of culture-collection strains from BacDive [@bacdive], NCBI RefSeq targeted-loci 16S sequences of type material, and 16S genes from isolate genomes in GTDB release 232 [@gtdb]; sequences from metagenome-assembled or single-cell genomes were excluded. {{gtdb_dropped}} of {{gtdb_genomes}} isolate genomes ({{gtdb_dropped_pct}}) carried 16S genes from unrelated phyla and were removed. For each study, the reference was cut in silico to the study's primer pair, and each ASV's identity to its best hit was computed by global alignment with all gaps counted (VSEARCH `--iddef 1`). ASVs without a hit at 75% or above were "no hit". Identity ignoring terminal gaps and gap-blind identity were sensitivity analyses.

## Outcome and unit of analysis

The unit of analysis was one ASV in one usable plating condition (study, source sample, medium, oxygen, temperature, incubation time), restricted to ASVs present in the condition's source sample (at least 2 environmental reads). A condition was usable if its pooled plate reads and its source's pooled environmental reads both reached 1,000 prokaryotic reads. The primary outcome, enriched on plates, required at least 5 reads in the condition's plates and a plate relative abundance at least twice the source relative abundance: DNA carried over from the inoculum should track its source abundance, while growth should raise it. Detection (at least 2 reads) was a secondary outcome. The shortcut counted an ASV as having a cultured relative if its identity was at least 94.5%.

Media were classed by one ordered rule from recipe facts: selective (antibacterial agents), habitat extract (soil, sludge, rumen fluid or similar), defined, dilute complex (at most 4 g/L organic matter) and rich complex. Temperatures were banded as cold (20 °C or below), moderate, or warm (33 °C or above).

## Calibration model

The primary model was a mixed-effects logistic regression fitted with glmmTMB [@glmmtmb] in R: a natural cubic spline of identity (boundary knots 75 and 100%, interior knots 90, 94.5, 97 and 99%) with a no-hit indicator, log10 source relative abundance, log10 plate pool depth, medium class, temperature band, oxygen and log2 incubation days as fixed effects, with sum-to-zero coding so that an unknown or unseen level counts as the average of the known levels, and random intercepts for study, source sample, plating condition and phylum. Phylum was registered as a fixed effect; because two phyla had almost no enrichment events and their estimates diverged, it was moved to a random intercept before any held-out metric was computed (Deviation 1). Predictions for a new study average over the study, source and condition random effects by Gauss-Hermite quadrature. Companion models were binned rates, an identity-only model, a model without identity, and gradient boosting (LightGBM) with fixed settings.

## Validation

Every model was refitted with each study group held out and evaluated on that group. Calibration was measured by the observed-to-expected ratio and the calibration slope (registered ranges 0.67 to 1.5 and 0.8 to 1.25), discrimination by the AUC averaged over groups, with 95% intervals from 1,000 bootstrap resamples of study groups. The curve's intervals come from 1,000 draws of the fixed effects; the between-study spread is the 10th to 90th percentile of the study random effect. Recovery shares were bootstrapped over source samples within each study.

## Medium, curation and taxon analyses

Medium-by-taxon effects were estimated in the studies with two or more known medium classes by adding random family, family-by-study and family-by-medium-class intercepts to the primary model, tested by a likelihood-ratio test and by held-out log loss. For the curation test, each ASV's tied best hits were linked to BacDive strains, each strain's growth media were classed from MediaDive recipes [@mediadive] by the same rule, and the share of linked strains recorded on the plate's own class entered the model beside the shares recorded on each class. Under-recovered taxa were found by comparing observed enrichments with a model without taxonomy, shifted to each held-out study's own total, with a consistency rule across studies.

## Survey application

Survey studies were selected from MGnify [@mgnify] by biome lineage, in a seeded random order, with rules on title words, run numbers, read depth and a primer screen of the first reads, fixed before any sample was processed. Amplicon reads from ENA went through the training pipeline; for shotgun metagenomes, MGnify v5.0 reads encoding SSU rRNA were placed against SILVA cut to the V4 window, trimmed to the covered part, and compared with the cultured V4 reference with terminal gaps excluded. Predictions averaged over the training plating conditions with each study group weighted equally. The range used the 10th and 90th percentiles of the study effect; the shortcut was judged to overstate recovery in a habitat if the lower 2.5% bound of the mean difference between the shortcut share and the upper end was above 0 and the difference was positive in at least two thirds of its studies. Before the surveys, the procedure was checked on each training study with a model that had not seen it. Published estimates were copied from the papers by script.

## Software

Python {{ver_python}} (pandas {{ver_pandas}}, NumPy {{ver_numpy}}, SciPy {{ver_scipy}}, LightGBM {{ver_lightgbm}}, matplotlib {{ver_matplotlib}}), R {{ver_r}} with glmmTMB {{ver_glmmtmb}}, VSEARCH {{ver_vsearch}} and cutadapt {{ver_cutadapt}}.

# Data and code availability

All raw reads are public under the accessions in Supplementary Table S1. The code, configuration, analysis plan and processed tables are at https://github.com/oswan-1225/plate-truth; the full pipeline runs from the numbered scripts. The harmonized paired dataset (census, condition metadata, ASV tables and nearest-relative table) is packaged for a Zenodo release with a data dictionary.

# Supplementary Table S1: studies analysed

| Study | Habitat | Accession | Region | Sources | Conditions |
|---|---|---|---|---|---|
{{table_studies}}

# Figure legends

![](../figures/fig1_census.png)

**Figure 1. The studies analysed.** (a) Plating conditions analysed (log scale) and number of source samples per study, coloured by habitat. (b) Medium classes, temperature bands and oxygen regimes covered by each study's analysed rows (filled), conditions whose class could not be resolved (open), and the amplified 16S region. The two Zheng et al. datasets form one study group.

![](../figures/fig2_calibration.png)

**Figure 2. Identity to the nearest cultured relative ranks growth in a new study but does not set its level.** (a) Probability that an ASV is enriched on plates (at least 5 reads and at least twice its source relative abundance) against identity to the nearest cultured relative, for a study the model has not seen, averaged over the analysed communities and conditions. Line: estimate; band: 95% interval of the fixed effects; dashed: 10th and 90th percentile study. Orange: observed share enriched per identity bin, pooled over studies, each predicted only from the other studies (95% intervals); black ticks: mean held-out prediction. ASVs without a cultured relative at 75% or above are drawn at 77%. (b) Held-out AUC per study group. (c) Held-out observed / expected enrichments per study group; grey: registered calibration range.

![](../figures/fig3_recovery.png)

**Figure 3. The shortcut overstates recovery in every study.** (a) Per study group, the share of the source community's reads with a cultured relative at 94.5% or above (the shortcut), enriched in any plating condition, enriched in the median condition, and predicted for the median condition by the held-out model. (b) Shortcut share divided by the share enriched in any condition, with 95% intervals from resampling source samples; studies with one source sample have no interval (open).

![](../figures/fig4_shifts.png)

**Figure 4. What shifts the curve.** (a) The curve for environmental and host-associated habitats, from a model with an identity-by-habitat-type interaction; bands: 95% intervals. (b) Held-out calibration slope per study group with 95% intervals; above 1, identity separates growers more sharply in that study than the model predicts. (c) Odds ratios with 95% intervals for medium class, temperature band and oxygen relative to the average level of each factor, and for source abundance, plate depth and incubation time, from the full-data model.

![](../figures/fig5_medium_taxa.png)

**Figure 5. Medium-class preferences by taxon.** Family-by-medium-class effects (a) and phylum-by-medium-class effects (b) on the log-odds scale, in the five studies that plated on two or more medium classes, adjusted for relatedness, source abundance, depth and conditions. Blue: enriched more often on that class than the taxon's own average and the class average; red: less often. Dots: 95% interval excludes 0. Blank: no rows.

![](../figures/fig6_application.png)

**Figure 6. The shortcut against predicted recovery in public surveys.** Per habitat and survey type, predicted share of reads in taxa enriched on one plate (dark bar, from a low- to a high-recovery study) and on any of 8 plates (light bar), with 95% intervals for the ends, and the shortcut share (diamond, 95% interval over studies). Filled diamond: the shortcut exceeds the one-plate upper end by the registered rule; ringed: it also exceeds the any-of-8 upper end. Predictions are averaged over the plating conditions of the training studies. Habitats marked "extrapolation" have no training study.

![](../figures/fig_phase4_curation.png)

**Supplementary Figure S1. The curation test.** Odds ratio for the match term (all of a taxon's linked collection strains recorded on the plate's own medium class against none) with 95% intervals: the registered model and its sensitivity analyses, the matched-plate Mantel-Haenszel odds ratio, and (grey, exploratory) refits with recorded-class effects per study and within the two largest studies alone. Dashed line: no effect.

![](../figures/fig_phase4_under_recovered.png)

**Supplementary Figure S2. Families recovered less or more often than relatedness predicts.** Observed / expected enriched rows per family against a model of identity, source abundance and conditions without taxonomy, fitted with each study held out and shifted to that study's total, with 95% intervals from resampling studies. Red: systematically under-recovered; blue: over-recovered; grey: neither. Triangles: no enriched rows, drawn at 0.003.

![](../figures/fig_phase5_training_check.png)

**Supplementary Figure S3. The survey procedure applied to the training studies.** For each training study, with a model that had not seen it: predicted one-plate range (dark bar), predicted range on any of the study's own number of plates (light bar), observed per-plate recovery (filled dot), observed recovery on any plate (open dot) and the shortcut share (diamond). Symmetric log scale below 0.01.

![](../figures/fig_phase5_composition.png)

**Supplementary Figure S4. Survey composition by the family calls.** Share of each habitat's reads (mean over studies) in families called under-recovered, over-recovered or neither, and in families that could not be evaluated.

# References

{{references}}
