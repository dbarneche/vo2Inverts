**SUPPLEMENTARY MATERIAL**  

*Estimating metabolic rates at $C_{CO_2}$*  
We use the species-specific parameters estimated and calculated in the main text ($C_{50\%\dot{V}\textrm{O}_{\textrm{2}max}}$, $\dot{V}\textrm{O}_\textrm{2}_{max}$ and $C_{CO_2}$) to back calculate what the $\dot{V}\textrm{O}_\textrm{2}$ of each individual will be at their respective $C_{CO_2}$. To do so, for each $i$ individual, we use the following the formula:  
  
$$\dot{V}\textrm{O}_\textrm{2}_{at C_{CO_2}~i} = \frac{\dot{V}\textrm{O}_\textrm{2}_{max~s} \times Y_{max~i} \times C_{CO_2~s}}{C_{50\%\dot{V}\textrm{O}_{\textrm{2}max~s}} + C_{CO_2~s}}$$ (S1)  
  
where $\dot{V}\textrm{O}_\textrm{2}_{max~s}$ and $C_{CO_2~s}$ are species-specific parameters, with $s$ representing a species, and $Y_{max~i}$ representing the highest measurement of individual $i$ (i.e. same value used to standardise $\dot{V}\textrm{O}_\textrm{2}$ between 0 and 1 in the first place).  
  
*Estimating mass-independent metabolic rates at $C_{CO_2}$*  
Using the individual-specific calculated metabolic rates above ($\dot{V}\textrm{O}_\textrm{2}_{at C_{CO_2}}$ in ml O~2~ h^-1^), we estimated species-specific mass-independent (i.e normalized) metabolic rates, $V_o$ (ml O~2~ g$~^{\alpha}$ h^-1^), following the equation:  
  
$$\dot{V}\textrm{O}_\textrm{2}_{at C_{CO_2}} = V_o m^{\alpha}$$ (S2),  
  
where $\alpha$ is the dimensionless mass-scaling exponent of metabolic rates and $m$ is body mass. We take the natural logarithm at both sides of equation S2 to fit the model:  
  
$$\textrm{ln}\dot{V}\textrm{O}_\textrm{2}_{at C_{CO_2}} = [\textrm{ln}V_o + {\Delta}_s \textrm{ln}V_o] + [\alpha + {\Delta}_s \alpha] \textrm{ln}m$$ (S3),  
  
where ${\Delta}_s$ represents log-additive species-specific deviations in $\textrm{ln}V_o$ and $\alpha$.  
  
*Model fitting*  
We follow the exact same procedure already described in the main text. We fit equation S3 above in Bayesian framework by calling *JAGS* version 4.2.0 from the R package *R2jags* version 0.05-6 [@r2jags] in order to derive posterior distributions and associated 95\% credible intervals (CIs) for the fitted parameters, $\textrm{ln}V_o$, ${\Delta}_s \textrm{ln}V_o$, $\alpha$ and ${\Delta}_s \alpha$. Random effects, ${\Delta}_s \textrm{ln}V_o$ and ${\Delta}_s \alpha$, were assumed to be normally distributed, with means of 0. Fitted parameters were assigned priors that were vague (i.e. locally uniform over the region supported by the likelihood) [@kruschke2014book]. The posterior distributions of model parameters were estimated using Markov chain Monte Carlo (MCMC) methods by constructing three chains of 1.5 $\times$ 10^6^ steps each, including 7.5 $\times$ 10^5^-step burn-in periods. Chains were thinned using a 375-step interval, so a total of 6,000 steps were retained to estimate posterior distributions (i.e. 3 $\times$ (1.5 $\times$ 10^6^ - 7.5 $\times$ 10^5^)/375 = 6,000).  
  
Finally, using the species-specific parameters estimated using equation S3 above, we estimate the oxygen consumption of native and invasive species at their respective $C_{CO_2}$ values estimated using equation 1. Means and associated 95% credible intervals were obtained by using the full posterior MCMC sample matrix from equation S3.  
