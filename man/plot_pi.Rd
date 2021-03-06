% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pi.R
\name{plot_pi}
\alias{plot_pi}
\title{Plot the estimated proportions of US or USA counts in each group}
\usage{
plot_pi(results, gene_id, cluster_id, type = "USA", CI = TRUE, CI_level = 0.95)
}
\arguments{
\item{results}{a \code{list} of 3 \code{\linkS4class{data.frame}} objects, computed via \code{\link{DifferentialRegulation}}.}

\item{gene_id}{a character, indicating the gene to plot.}

\item{cluster_id}{a character, indicating the cell cluster to plot.}

\item{type}{a character (either 'SU' or 'SUA').
If "SU", it plots the proportion of Spliced and Unspliced counts 
(Ambiguous counts are assigned 50:50 to Spliced and Unspliced counts).
If "SUA" (default), it plots the proportion of Spliced, Unspliced and Ambiguous counts 
(Ambiguous counts are kept separately).
Note that, although US reads are easier to interpret, 
USA reads more closely reflect what is being compared between conditions.}

\item{CI}{a logical ('TRUE' by default), indicating whether to plot a 
profile Wald type confidence interval around the estimated proportions.}

\item{CI_level}{a numeric between 0 and 1, indicating the level of the confidence interval.}
}
\value{
A \code{ggplot} object.
}
\description{
\code{plot_pi} plots the posterior means of the proportions
of US (if 'type' = 'US') or USA (if 'type' = 'USA') counts, in each group.
If 'CI' is TRUE, a profile Wald type confidence interval will also be added;
the level of the confidence interval is specified by 'CI_level'.
}
\examples{
# load internal data to the package:
data_dir = system.file("extdata", package = "DifferentialRegulation")

# specify samples ids:
sample_ids = paste0("organoid", c(1:3, 16:18))
# set directories of each sample input data (obtained via alevin-fry):
base_dir = file.path(data_dir, "alevin-fry", sample_ids)
file.exists(base_dir)

# set paths to USA counts, cell id and gene id:
# Note that alevin-fry needs to be run with '--use-mtx' option
# to store counts in a 'quants_mat.mtx' file.
path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")

# load USA counts:
sce = load_USA(path_to_counts,
               path_to_cell_id,
               path_to_gene_id,
               sample_ids)
 
# define the design of the study:
design = data.frame(sample = sample_ids,
                    group = c( rep("3 mon", 3), rep("6 mon", 3) ))
design

# cell types should be assigned to each cell;
# here we load pre-computed cell types:
path_to_DF = file.path(data_dir,"DF_cell_types.txt")
DF_cell_types = read.csv(path_to_DF, sep = "\t", header = TRUE)
matches = match(colnames(sce), DF_cell_types$cell_id)
sce$cell_type = DF_cell_types$cell_type[matches]
               
PB_counts = compute_PB_counts(sce = sce,
                              EC_list = NULL,
                              design =  design,
                              sample_col_name = "sample",
                              group_col_name = "group",
                              sce_cluster_name = "cell_type",
                              min_cells_per_cluster = 100, 
                              min_counts_per_gene_per_group = 20)
                              
# Differential regulation test based on estimated USA (unspliced, spliced, ambiguous) counts
set.seed(169612)
results_USA = DifferentialRegulation(PB_counts, EC = FALSE)

# DifferentialRegulation returns of a list of 3 data.frames:
# "Differential_results" contains results from differential testing only;
# "US_results" has estimates and standard deviation (SD) for pi_S and pi_U (proportion of Spliced and Unspliced counts);
# "USA_results" has estimates and standard deviation (SD) for pi_S, pi_U and pi_A (proportion of Spliced, Unspliced and Ambiguous counts).
names(results_USA)

# We visualize differential results:
head(results_USA$Differential_results)

# For improved performance, at a higher computational cost,
# we recommend using equivalence classes (EC) (here not run for computational reasons)
# see help(DifferentialRegulation) examples.

# plot top (i.e., most significant) result:
# plot USA proportions:
plot_pi(results_USA,
        type = "USA",
        gene_id = results_USA$Differential_results$Gene_id[1],
        cluster_id = results_USA$Differential_results$Cluster_id[1])

# plot US proportions:
plot_pi(results_USA,
        type = "US",
        gene_id = results_USA$Differential_results$Gene_id[1],
        cluster_id = results_USA$Differential_results$Cluster_id[1])

}
\seealso{
\code{\link{DifferentialRegulation}}
}
\author{
Simone Tiberi \email{simone.tiberi@uzh.ch}
}
