# EPT
Model scripts for expedited partner therapy analysis

# Evaluating Expedited Partner Therapy among Men who have Sex with Men (MSM) in the United States

This repository holds the source to code to reproduce the analysis featured in our STI and HIV transmission model among men who have sex with men in the United States. This study investigated the possible implementation of EPT for MSM in the United States, evaluating different scale-up strategies along an EPT continuum -- from provision of EPT medication to index partners by health care providers, provision of EPT medication to non-index partners bys index patients, uptake of EPT medication by non-index partners, and treatment efficacy among non-index partners -- and the impact on STI incidence among MSM.

## Citation

> Weiss KM, Jones JS, Katz DA, Gift TL, Bernstein K, Workowski K, Rosenberg ES, Jenness SM. Epidemiological impact of expedited partner therapy for men who have sex with men: A modeling study. _Under Review._

<img src="https://github.com/EpiModel/EPT/raw/master/analysis/Fig1.png">

## Abstract
Expedited partner therapy (EPT) is an intervention for patients with gonorrhea (NG) or chlamydia (CT), providing index patients with prescriptions or medication to give to their partners. EPT is recommended for heterosexuals but not for men who have sex with men (MSM), partially due to concerns about overtreatment of uninfected partners and missed opportunities for HIV diagnosis. We extended our stochastic network-based mathematical model of HIV, NG, and CT among MSM to include EPT. EPT implementation was simulated for 10 years. Counterfactual scenarios varied EPT coverage, provision, uptake, and partnership window duration. We estimated sexually transmitted infection (STI) incidence, proportion of infections averted (PIA), and process outcomes under each scenario. Delivery of EPT to 20% of eligible MSM index patients (coverage) reduced cumulative STI incidence by 27% (interquartile range (IQR): 13%-39%) over 10 years compared to current estimated STI screening levels. A 20% increase in providing medication to non-index partners (provision) averted 32% (IQR: 20%-41%) of STI infections compared to estimated STI screening levels. When targeted by partnership type, EPT solely to casual partners maximized the population infections averted. The proportion of partners given medication who had no current STI varied from 52% to 63%, depending on coverage level. The proportion of partners given medication with undiagnosed HIV infection was 4% across scenarios. EPT could reduce bacterial STI incidence for MSM. However, this intervention could result in missed opportunities for HIV/STI prevention and a substantial increase in use of antimicrobials by STI-uninfected MSM, raising concerns about cost and antimicrobial resistance. 

<img src="https://github.com/EpiModel/EPT/raw/master/analysis/Fig2.png">

## Model Code

These models are written and executed in the R statistical software language. To run these files, it is necessary to first install our epidemic modeling software, [EpiModel](http://epimodel.org/), and our extension package specifically for modeling HIV and STI transmission dynamics among MSM, [EpiModelHIV](http://github.com/statnet/EpiModelHIV).

In R:
```
install.packages("EpiModel", dep = TRUE)

# install remotes if necessary, install.packages("remotes")
remotes::install_github("statnet/tergmLite")
remotes::install_github("statnet/EpiModelHPC")
remotes::install_github("statnet/EpiModelHIV", ref = "ept")
```

