# evee-tools
Tools for comparative RNA-Seq analysis. Data coming soon!
<br><br>
`makeExpressionMatrix.py` - Given a file of ortholog groups and list of rsem output files, creates a matrix file (gene x sample) of TPM values, and an index file mapping each column to the corresponding species.
<br><br> *Note that the following R scripts depend on edgeR and ouch libraries.*
<br><br>`residuals.R` - Given outputs of `makeExpressionMatrix.py`, calculates the expression residuals between each species and a reference species.

`fitOUModel.R`- Fits a Brownian motion and Ornstein-Uhlenbeck model, and the performs likelihood ratio test between the two nested models.

`scoreGenes.R` - Calculates a z-score for expression levels from new data set.

`ouRegimes.R`- Fits hypotheses of multivariate Ornstein-Uhlenbeck models and calculates AIC and BIC scores.
