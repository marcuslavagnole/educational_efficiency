In this repository, one can find the R routines used in the article [Implications of Brazilian institutional guidelines on educational efficiency](https://doi.org/10.1353/eco.2020.0009). It is a joint work with Kalinca L. Becker and Mario Jorge Mendonça, published in _Economía - LACEA journal_.

The paper investigates the relationship between inefficiency in the Brazilian education system and municipal wealth, discussing how the current legislation possibly influences it. To that end, we apply a stochastic frontier model that accommodates covariates in the asymmetric error component to analyze the impact of per capita GDP on inefficiency. This methodology is applied to a data set on the Rio Grande do Sul municipalities for the years 2007 and 2017. The results indicate a positive effect, suggesting wealthier municipalities are less efficient in allocating resources.

This repo includes:

- **MCMC_Code.R** : Main file containing the MCMC routine; 
- **Full_Conditionals.R** : Auxiliary file with all full conditional distributions;
- **dados.RData** : real data set;
- **W_RS.RData** : adjacency matrix for the municipalities of Rio Grande do Sul state;
- **W_RS_imediata.RData** : adjacency matrix for the immediate regions of Rio Grande do Sul;
- **W_RS_intermediaria.RData** : adjacency matrix for the intermediate regions of Rio Grande do Sul.
