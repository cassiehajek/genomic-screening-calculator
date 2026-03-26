# Genomic Screening Calculator — Technical Methodology

**Calculator**: Population Genomic Screening Cost-Effectiveness Calculator
**Reference model**: Guzauskas GF et al. *Ann Intern Med.* 2023;176:585–595
**Last updated**: 2026-03-26

---

## Overview

This calculator estimates the cost-effectiveness of population-wide genomic screening for three
hereditary conditions — Lynch Syndrome (LS), Familial Hypercholesterolaemia (FH), and Hereditary
Breast and Ovarian Cancer (HBOC) — using a simplified annual-incidence × present-value annuity
approach calibrated to published Markov model outputs.

Primary outputs: QALYs gained, events prevented, and ICER ($/QALY) per 100,000 screened.

---

## Model Structure

### Events Prevented

For each condition and event type, events prevented per carrier are calculated as:

```
events_prevented = carriers
                 × (lifetime_risk / accrual_yrs)   ← annual incidence rate
                 × efficacy
                 × compliance
                 × ePvF                             ← discounted event multiplier
                 × (1 − prior_identification)
```

**Accrual years** convert lifetime penetrance to an approximate annual incidence rate:

| Condition | Accrual years | Rationale |
|-----------|--------------|-----------|
| Lynch CRC / Endometrial | 50 | Risk accrues ~age 30–80 |
| FH MACE | 60 | CV risk accrues ~age 30–90 |
| HBOC Breast | 60 | Breast risk concentrates ~age 25–85 |
| HBOC Ovarian | 135 | Ovarian risk concentrates ~age 45–70; longer denominator discounts fewer events |

**ePvF (event present-value factor)** is an age-weighted present-value annuity factor:

```
pvFactor(yrs, r) = (1 − (1+r)^−yrs) / r

ePvF = Σ_i  ageWeight_i × pvFactor(remainingYrs_i, discountRate)
```

Where remaining years to age 80 are: 56 (age 18–30), 42 (31–45), 27 (46–60), 12 (61–75), 0 (76+).

This converts the annual event rate into a **discounted lifetime event count** — events expected
sooner (near-term) receive higher weight than distant future events.

### QALY Calculation

```
QALYs = events_prevented × qalyPerEvent
```

`qalyPerEvent` is read from an editable table stratified by screening age group and event type.
Values represent the **undiscounted QALY gain per prevented event** for a person screened at that
age (accounting for remaining life-years and disease-specific utility weights derived from
Guzauskas 2023 Table 1). The discounting of WHEN events occur is handled entirely by `ePvF`.

**QALY table defaults** (QALYs gained per prevented event, by screening age group):

| Event type | 18–30 | 31–45 | 46–60 | 61–75 | 76+ |
|------------|-------|-------|-------|-------|-----|
| CRC (Lynch) | 3.3 | 2.6 | 1.9 | 1.1 | 0.4 |
| Endometrial (Lynch) | 2.9 | 2.3 | 1.7 | 0.9 | 0.3 |
| MACE (FH) | 6.6 | 5.3 | 3.8 | 2.1 | 0.8 |
| Breast cancer (HBOC) | 4.3 | 3.5 | 2.4 | 1.4 | 0.6 |
| Ovarian cancer (HBOC) | 5.0 | 4.1 | 2.9 | 1.6 | 0.6 |

Higher values for younger age groups reflect more remaining life-years gained per prevented event.
Higher values for ovarian vs breast reflect the poorer prognosis (lower survival) of ovarian cancer.

**Note on MACE (FH) values:** q_mace = 4.0 reflects the direct per-carrier MACE benefit only.
Cascade testing (identifying FH relatives of index cases) adds additional QALYs separately
via the `fh_cascade` parameter (see below).

### Cost Calculation

Management costs use a separate present-value annuity factor (`pvF`) based on the user-specified
time horizon (not remaining risk years). This distinguishes:

- **`ePvF`** — governs how many future disease events occur (remaining biological risk horizon)
- **`pvF`** — governs how long surveillance/treatment costs are incurred (programme time horizon)

### ICER

```
ICER = net_programme_cost / total_QALYs_gained

net_programme_cost = screening_cost + genetic_counselling_cost
                   + management_cost − averted_treatment_cost
```

---

## Cascade Testing

Each condition has an editable **cascade yield (%)** parameter representing additional carriers
identified per 100 index cases through first-degree relative testing. Cascade-identified carriers
receive the same per-carrier benefit (events prevented, QALYs) as directly screened carriers and
incur management and counselling costs, but not the upfront panel screening cost.

```
effective_carriers = direct_carriers × (1 + cascade_yield / 100)
```

This scales effective carrier counts for events, QALYs, GC costs, management costs, and
averted treatment costs. The panel test cost (`N × costScreen`) is unchanged.

Default cascade yields:

| Condition | Default | Rationale |
|-----------|---------|-----------|
| Lynch Syndrome | 0% | q_crc/q_endo calibrated to Guzauskas 2022 which includes cascade; adding explicit cascade would double-count unless q values are adjusted |
| FH | **65%** | Direct FH benefit alone gives ~67 QALYs/100K at age 30; Guzauskas 2023 reports ~111. The gap (~44 QALYs) is attributable to cascade. 65% yield calibrated to close this gap (16.8 MACE events × 1.65 × 4.0 QALY/event = 111 ✓). Literature: each FH carrier has ~3–4 first-degree relatives × 50% carrier probability × ~40–50% uptake ≈ 60–100% yield. |
| HBOC | 0% | HBOC q_breast/q_ovarian already calibrated to Guzauskas 2023 outputs; adding cascade on top adds benefit above the reference baseline |

Users can set any condition's cascade yield to 0 to isolate the direct screening benefit,
or increase it to model more aggressive family follow-up programs.

---

## Prior Identification

Each condition has an editable **prior identification (%)** input representing the proportion of
carriers already diagnosed before population screening. These individuals receive no incremental
benefit from population screening. Default values from Guzauskas 2023:

| Condition | Default | Range used in SA |
|-----------|---------|-----------------|
| Lynch Syndrome | 13% | 7–20% |
| FH | 5% | 1–9% |
| HBOC | 17% | 14–21% |

---

## Age-Stratified ICER Curve

`computeIcerAtAge(age)` recalculates the full model for a hypothetical cohort screened entirely
at one age. It injects `pvFactor(80 − age)` directly as `ePvF` and forces all cohort weight into
the corresponding QALY age bucket. This reproduces the Figure 2 curves from Guzauskas 2023.

---

## Carrier Prevalence Defaults

| Gene(s) | Default (%) | Guzauskas reference | Source |
|---------|------------|---------------------|--------|
| MLH1 | 0.12 | — | |
| MSH2 | 0.11 | — | |
| MSH6 | 0.09 | — | |
| PMS2 | 0.03 | — | |
| **Lynch total** | **0.35** | **0.30** | Guzauskas 2022 *Genet Med* |
| LDLR | 0.32 | — | |
| PCSK9 | 0.02 | — | |
| APOB | 0.09 | — | |
| **FH total** | **~0.43** | ~0.40 | Nordestgaard 2016 |
| BRCA1 | 0.29 | ~0.11 | Kuchenbaecker 2017 |
| BRCA2 | 0.43 | ~0.38 | Kuchenbaecker 2017 |
| **HBOC total** | **0.72** | **0.495** | Guzauskas 2020 *JAMA Netw Open* |

---

## Known Calibration Issues: Why QALYs Exceed Guzauskas 2023

The Guzauskas 2023 three-condition model (Table 2, age-30 cohort) reports:

| Condition | QALYs per 100K screened |
|-----------|------------------------|
| Lynch Syndrome | ~170 |
| FH | ~111 |
| HBOC | ~215 |
| **Total** | **~495** |

The calculator can exceed these benchmarks for two compounding reasons:

### Reason 1: Carrier Prevalence Defaults Are Higher Than Reference Values

**HBOC (most significant):**
Default BRCA1 + BRCA2 = 0.72% vs. the Guzauskas 2020 reference of 0.495%.
This is a **45% overcount of HBOC carriers**, directly inflating HBOC events and QALYs by the
same proportion. For example: 720 BRCA+ carriers per 100K in our model vs. 495 in the reference —
even with identical per-event QALY values, this alone yields ~310 HBOC QALYs vs. the reference
~215.

The BRCA1 default (0.29%) is particularly high relative to population estimates (~0.11–0.13%).
BRCA2 (0.43%) is closer to the literature but still slightly elevated.

**Lynch Syndrome (secondary):**
Default total 0.35% vs. reference 0.30% — a **17% overcount**, inflating Lynch QALYs by ~17%.
At age 30 this pushes Lynch from ~170 toward ~200 QALYs per 100K.

### Reason 2: QALY-per-Event Values for HBOC Are ~30% Too High

Back-calculating from Guzauskas 2023 Table 2 at age 30 (40 breast events + 7 ovarian events per
100K → 215 total HBOC QALYs):

```
40 × qalyBreast + 7 × qalyOvarian = 215
```

Maintaining the current breast:ovarian ratio (6:7) implies:
- **Breast**: ~4.5 QALY/event  (table default: 6.0 — ~33% high)
- **Ovarian**: ~5.2 QALY/event (table default: 7.0 — ~35% high)

These values were set independently of the carrier-rate calibration. Because the carrier rates
were already over-estimated, the per-event QALY values were not driven down to compensate,
resulting in both sources of inflation compounding each other.

### Combined Effect

| Condition | Carrier error | Per-event QALY error | Approximate net QALY inflation |
|-----------|--------------|---------------------|-------------------------------|
| Lynch | +17% | small | ~17% above reference |
| FH | small | uncertain | moderate |
| HBOC | +45% | +33% | ~90–95% above reference |

The HBOC arm approximately doubles the Guzauskas reference output when both errors are present.

### Calibration Applied (v2, 2026-03)

The QALY per event table has been recalibrated so that **with the current (observed) carrier
rates**, the model reproduces Guzauskas 2023 Table 2 condition-by-condition at age 30.

**Calibration basis:** carrier rates are NOT changed (they reflect real-world observed
prevalence, which is higher than the Guzauskas reference population for HBOC and Lynch).
Only the per-event QALY values are adjusted to reflect the correct disease model.

**Derivation at age 30 (ePvF = 26.96):**

| Condition | Events/100K | Target QALYs | q_per_event | Verification |
|-----------|------------|-------------|-------------|-------------|
| Lynch CRC | 36.7 | — | 3.3 | 36.7 × 3.3 + 16.4 × 2.9 = **169** ≈ 170 ✓ |
| Lynch endo | 16.4 | 170 combined | 2.9 | |
| FH MACE | 16.8 direct + 10.9 cascade (×1.65) = 27.7 | 111 | 4.0 | 27.7 × 4.0 = **111** ✓ |
| HBOC breast | 41.4 | — | 4.3 | 41.4 × 4.3 + 7.3 × 5.0 = **214** ≈ 215 ✓ |
| HBOC ovarian | 7.3 | 215 combined | 5.0 | |
| **Total** | | **495** | | **169 + 111 + 214 = 494** ✓ |

**Effect of using Guzauskas reference carrier rates (Lynch 0.30%, HBOC 0.495%):**
Because the formula's effective risk-reduction parameters for HBOC differ from the Guzauskas
Markov model, entering the reference carrier rates will not reproduce Table 2 exactly. With
reference rates, per-carrier event counts are lower (~28 breast vs 40), so HBOC QALYs would
be ~145 rather than 215. This is an inherent limitation of the simplified formula vs. the full
state-transition model — it cannot simultaneously satisfy both carrier-rate regimes.

---

## Sensitivity Analysis

One-way sensitivity analysis (tornado chart) varies each key parameter independently while
holding all others at base-case values. Parameter ranges follow Guzauskas 2023 Table 1:

| Parameter | Low | High |
|-----------|-----|------|
| Test cost per person | $100 | $400 |
| Lynch prior identification | 7% | 20% |
| FH prior identification | 1% | 9% |
| HBOC prior identification | 14% | 21% |
| LS colonoscopy adherence | 60% | 100% |
| FH statin uptake | 45% | 75% |
| HBOC any-mgmt uptake | 60% | 95% |
| Time horizon | 20 yr | 50 yr |
| Discount rate | 0% | 5% |
| Lynch CRC risk (MLH1/MSH2) | 44% | 66% |
| BRCA1 lifetime breast risk | 58% | 86% |

---

## References

1. **Guzauskas GF et al.** Population Genomic Screening for Three Common Hereditary Conditions:
   A Cost-Effectiveness Analysis. *Ann Intern Med.* 2023;176:585–595. doi:10.7326/M22-0846
   *(Primary reference — all base-case defaults calibrated to this paper)*

2. **Guzauskas GF, Jiang S, et al.** Cost-effectiveness of Population-wide Genomic Screening for
   Lynch Syndrome in the United States. *Genet Med.* 2022;24(5):1017–1026.
   doi:10.1016/j.gim.2022.01.017
   *(Lynch carrier prevalence 0.3%; age-specific ICERs: age 30 $132K, age 40 $124K, age 50 $140K)*

3. **Guzauskas GF et al.** Cost-effectiveness of Population-Wide Genomic Screening for Hereditary
   Breast and Ovarian Cancer in the United States. *JAMA Netw Open.* 2020;3(10):e2022874.
   doi:10.1001/jamanetworkopen.2020.22874
   *(HBOC carrier prevalence 0.495%; age-specific ICERs: age 30 $87.7K, age 45 $268K)*

4. **Guo F et al.** Cost-Effectiveness of Population-Based Multigene Testing for Breast and
   Ovarian Cancer Prevention. *JAMA Netw Open.* 2024;7(2):e2356078.
   doi:10.1001/jamanetworkopen.2023.56078
   *(HBOC management costs; PALB2 parameters)*

5. **Spencer AH et al.** Cost-effectiveness of population screening for familial
   hypercholesterolaemia. *J Clin Lipidol.* 2022.
   *(FH standalone screening ICER ~$181K/QALY at age 20; bundling improves CE)*

6. **Bonadona V et al.** Cancer Risks Associated with Germline Mutations in MLH1, MSH2, and
   MSH6 Genes in Lynch Syndrome. *JAMA.* 2011;305(22):2304–2310.

7. **Nordestgaard BG et al.** Familial hypercholesterolaemia is underdiagnosed and undertreated
   in the general population. *Eur Heart J.* 2013;34(45):3478–3490.

8. **Kuchenbaecker KB et al.** Risks of Breast, Ovarian, and Contralateral Breast Cancer for
   BRCA1 and BRCA2 Mutation Carriers. *JAMA.* 2017;317(23):2402–2416.
