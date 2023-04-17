# DorsalGallega_Wolves

Wolves data (Dorsal Gallega mountain ridge, Spain) and R + Nimble codes to fit the SCR models used in the paper: Estimating wolf (_Canis lupus_) densities using video camera-traps and spatial capture-recapture analysis, from José Jiménez, Daniel Cara, Francisco García Dominguez and Jose Angel Barasona (2023). 

- _SCR_NegBin.R_: negative binomial observation model with constant baseline detection rate.
- _SCR_NegBin_eps[j].R_: negative binomial observation model with a trap-level random effect in the baseline detection rate.
- _SCR_Poisson.R_: Poisson observation model with constant baseline detection rate.
- _SCR_Poisson_eps[i].R_: Poisson observation model with a individual-level random effect in the baseline detection rate.
- _SCR_Poisson_eps[j].R_: Poisson observation model with a trap-level random effect in the baseline detection rate.
- _SCR_Poisson_eps[i,j].R_: Poisson observation model with a individual-by-trap random effect in the baseline detection rate.
- _SCR_Poisson_eps[i,j,k].R_: Poisson observation model with a individual-by-trap-by occasion random effect in the baseline detection rate.

