# Regression Discontinuity: Causal Inference in Hospital Infection Control

**Project 3: Causal Inference: RD**

**Firm:** Sentinel Health Analytics  
**Team:** Muhammad Sawaiz Fatar, Victor Ostolaza, Vedaant Rath, Arpita Ram Samant, Sam Sheng

---

## The Question

Does mandatory isolation actually reduce hospital-acquired infections (HAI)?

---

## Background

**Metropolitan General Health System (MGHS)** is an 850-bed hospital network in Chicago serving 45,000 patients annually. Their current HAI rate is 4.8 per 1,000 patient-days, well above the benchmark of 3.5. Each infection costs roughly $40,000. They've invested $8M in isolation infrastructure with a $4M annual operating cost, and in March 2024 implemented a rule: patients with an Infection Risk Score of 70 or above are placed in isolation wards.

The problem with just comparing isolated vs. non-isolated patients is obvious. Sicker patients get higher risk scores and go to isolation. They might have worse outcomes simply because they're already sicker, not because isolation doesn't work. That's the selection bias problem and it's exactly why we need a smarter approach.

The RD solution: compare patients who scored just above 70 vs. just below 70. These patients are nearly identical in health status. The only difference is isolation assignment. That clean switch at 70 is what makes the comparison credible.

---

## Variables

| Variable | Definition |
|---|---|
| **Outcome (Y)** | Hospital-acquired infection — binary yes/no |
| **Running Variable** | Infection Risk Score (0–100, continuous) |
| **Cutoff** | Risk Score = 70 |
| **Treatment (D)** | Isolation ward placement (1 = isolated, 0 = regular ward) |
| **Unobserved (U)** | Latent severity — "how sick the patient really is," the stuff that's hard to capture perfectly in one neat variable |

---

## The Two Settings We Analyze

### Sharp RD : MGHS (Strict Enforcement)

At MGHS, the rule is followed exactly. Every patient above 70 goes to isolation, no exceptions. This is a **Sharp RD** = treatment is a deterministic function of the running variable.

### Fuzzy RD : Community Health Network (New Client)

**Community Health Network (CHN)** is a 300-bed hospital system serving rural areas outside Chicago with 25,000 patients annually. CHN faces real resource constraints: only 20 isolation beds, physician staffing varies by shift, and a Risk Score above 70 recommends isolation but doesn't guarantee it. High-risk patients may not get isolated when beds are full. Low-risk patients sometimes get isolated based on physician judgment.

The question CHN hires us to answer: is the isolation guideline still effective when it's not strictly enforced?

This is a **Fuzzy RD** — the threshold creates a strong nudge but not a deterministic rule, so we need to use the threshold as an instrument and run 2SLS.

---

## The Omitted Variable Problem (Why Fuzzy RD Needs IV)

At CHN, there's a confounding variable we can't observe: **Hospital Crowding**.

When the hospital is crowded, two things happen at once. First, isolation wards fill up, so high-risk patients who should be isolated end up in regular wards. Second, crowded conditions increase patient-to-patient contact, overwork staff, and raise infection risk across the board — regardless of isolation status. This means crowding affects both the isolation decision AND the outcome, which is the classic OVB problem. A naive comparison would make isolation look ineffective even if it genuinely works.

Our instrument is the 70-point threshold itself. The threshold strongly predicts isolation (below 70, only ~7% get isolated; above 70, ~86% do), but the threshold is not related to crowding. This is the variation we need for a fair estimate.

---

## Data

Both datasets are simulated with N = 1,000 patients. The data-generating process is built to match the real-world mechanics of each setting:

**Sharp (MGHS):** Latent severity (U) drives the risk score, the risk score strictly determines treatment at 70, and isolation reduces infection probability (true effect = negative). The analyst-facing dataset contains only `id`, `risk_score`, `D`, and `HAI` -> latent severity is hidden, as it would be in real data.

**Fuzzy (CHN):** Latent severity now also influences whether isolation actually happens (creating OVB). A fuzzy factors variable (bed availability, staffing, physician judgment) introduces non-compliance independently of U that's crucial for instrument validity. The analyst-facing dataset hides both U and the fuzzy factors.

---

## Methods

### Sharp RD

We first show the naive DiM estimate (raw comparison of isolated vs. non-isolated patients), which gives a misleading result because OVB hasn't been addressed. We then run regression controlling for the running variable, check for nonlinearity using a quadratic term and a joint F-test, and test for differing marginal effects by allowing different slopes on each side of the cutoff.

We run this analysis across multiple bandwidths (±3, ±5, ±10, ±15 points around 70) so we can report how stable the effect is as we vary the comparison window.

We also run the manipulation check (histogram of risk scores around the cutoff) to confirm there's no suspicious bunching just below 70, which would suggest doctors or patients are gaming the score to avoid isolation.

### Fuzzy RD

We use the threshold as an instrument and implement 2SLS (`ivreg`). The first stage regresses actual isolation on the threshold dummy and the risk score. We run a full spec search: nonlinearity check, interaction check (differing slopes), and a joint F-test to confirm the simpler linear specification is appropriate.

We then compare Naive OLS vs. 2SLS to show the direction of OVB, and run calipers (±3, ±5, ±10, ±15) to check robustness of the IV estimate near the threshold.

Diagnostic tests (Weak Instruments, Wu-Hausman) validate that the instrument is strong and that IV was actually necessary.

---

## Key Results

### Sharp RD : MGHS

When we zoom in to patients near the 70 cutoff, the sign on isolation is consistently negative, patients just above 70 have fewer infections than patients just below 70. The estimated effect ranges from about 4 to 13 percentage points depending on the bandwidth, and is statistically significant at ±15 (ATE ≈ −0.077). The simpler linear model is preferred — nonlinearity and interaction terms are both jointly unimportant.

For MGHS the takeaway is straightforward: the isolation rule looks like it's doing what you paid for. For every 100 patients near the cutoff, isolation prevents about 7–8 infections. At $40,000 per infection, that's roughly $280k–$320k saved per 100 borderline patients.

### Fuzzy RD : CHN

| Method | ATE Estimate | What it means |
|---|---|---|
| Naive OLS | -0.0112 | Isolation reduces HAI by 1.1% |
| 2SLS (IV) | +0.0418 | Isolation increases HAI by 4.2% |

These tell very different stories. OLS is biased because the hidden crowding factor is still contaminating the estimate. Once we correct for it using the threshold as an instrument, the picture flips. Neither estimate is strong enough to be conclusive (2SLS p = 0.107), but across all calipers the effect is consistently not clearly negative ranging from +3.5% to +8.1%.

Our recommendation to CHN: the flexibility introduced in the policy may be diluting its effectiveness. Either return to stricter enforcement of the rule, or collect more data before drawing firm conclusions.

---

## Repository Structure

```
.
├── Regression_Discontinuity_Report.qmd   # Main report
├── README.md
```

---

## Requirements

- R (>= 4.0)
- Packages: `tidyverse`, `ggdag`, `dagitty`, `ggplot2`, `AER`, `car`, `ivreg`

To render the report:

```r
quarto::quarto_render("Regression_Discontinuity_Report.qmd")
```
