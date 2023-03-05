# DorsalGallega_Wolves

Wolves data (Dorsal Gallega mountain ridge, Spain) and R + Nimble code to fit the SCR model used in the paper: Estimating wolf (_Canis lupus_) densities using video camera-traps and spatial capture-recapture analysis, from José Jiménez, Daniel Cara, Francisco García Dominguez and Jose Angel Barasona (2023). 

- _SCR_NegBin_RE.R_: negative binomial observation model with a trap-level random effect in the basal detection rate.
- _SCR_NegBin_RE_faster.R_: Faster version of the negative binomial model with a trap-level random effect  in the basal detection rate, with a custom Nimble function ("GetDetectionRate") and a vectorised negative binomial distribution ("dNBVector"), based on the work of Ben Augustine (https://github.com/benaug).
- _SCR_Poisson_RE.R_: Poisson observation model with a trap-level random effect in the basal detection rate.
- _SCR_Poisson.R_: Poisson observation model with constant basal detection rate.
