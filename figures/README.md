# Figures

For each figure, the filename tells what have been studied.

The only studied taxonomic divisions are ``BCT`` (Bacteria|Archea), ``HUMAN`` (_Homo sapiens_|_Homonidae_), ``MICE`` (_Mus_|_Muridae_) and ``PLN`` (Plant)

The prefix is either:
* metag [Only metagenomic indexes]
* no_metag [Only non-metagenomic indexes]

The following subword means what is the figure:
* accuracy [Measure of the number of accession that were well classified PER span and tax. division]
* confusion_matrix [Confusion matrix of predicted tax. divisions PER tax. division]
* timeout_proportion [Proportion of timed out NCBI metadata queries PER span and tax. division]
* unresolved_proportion [Proportion of failed NCBI metadata queries PER span and tax. division]
* unknown_proportion [Proportion of probably deleted NCBI SRA accession PER span and tax. division]

Then, another subwords:
* with_errors [Means that errors like timeout, unresolved and unknown are counted as a misclassification]
* without_errors [Means that only the accessions which metadatas were successfully retrieved are counted]
