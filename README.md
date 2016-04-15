## finemapped_eqtls

## installation

```R
library("devtools")
devtools::install_github('zaczap/finemapped_eqtls')
```

## examples

```R
library(finemappedEqtls)
eqtls = fetch_eqtl_data("BRCA1")
plot_multipopulation_eqtls(eqtls, populations=c('CEU','JPT'), start = 41000000, end = 41500000)
```

```R
eqtls = fetch_eqtl_data("ORMDL3")
plot_finemapped_eqtls(eqtls, 'JPT', start = 38000000, 38050000)
```

```R
eqtls = fetch_eqtl_data('ORMDL3')
plot_eqtl_locus(eqtls, populations=c('CHB','JPT'), discovery='JPT')
```
