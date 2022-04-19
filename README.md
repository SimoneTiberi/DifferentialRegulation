# DifferentialRegulation
*DifferentialRegulation* is a method for detecting differentially regulated genes between two groups of samples (e.g., healthy vs. disease, or treated vs. untreated samples), by targeting differences in the balance of spliced and unspliced mRNA abundances, obtained from single-cell RNA-sequencing (scRNA-seq) data.
*DifferentialRegulation* accounts for the sample-to-sample variability, and embeds multiple samples in a Bayesian hierarchical model.
In particular, when providing equivaelence classes data (via `EC_list`), reads that are compatible with multiple genes, or multiple splicing versions of a gene (unspliced spliced or ambiguous), are allocated to the gene of origin and their splicing version.
Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques (Metropolis-within-Gibbs).

## Alignment and quantification with *alevin-fry*
*DifferentialRegulation* inputs scRNA-seq data, aligned via *alevin-fry*.

NOTE: when using *alevin-fry*, set options `--d` (or `--dump-eqclasses`), to obtain the equivalence classes, and `--use-mtx`, to store counts in a `quants_mat.mtx` file (as expected by our `load_USA` function).

We also recommend using the `--CR-like-EM` option, which also allows equivalence classes of reads that map to multiple genes.
