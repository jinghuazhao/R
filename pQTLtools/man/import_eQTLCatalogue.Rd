\name{import_eQTLCatalogue}
\alias{import_eQTLCatalogue}
\title{Import eQTL Catalogue}
\usage{import_eQTLCatalogue(ftp_path, region, selected_gene_id, column_names, verbose = TRUE)}
\description{
An adopted function which imports eQTL Catalogue.
}
\arguments{
\item{ftp_path}{URL}
\item{region}{chr:start-end}
\item{selected_gene_id}{An ensGene ID}
\item{column_names}{Column names of the dataset}
\item{verbose}{Extra information}
}
\details{
This function is based on the eQTL-Catalogue-resources.
}

\value{
A summary statistic object
}

\section{References}{
See https://www.ebi.ac.uk/about/news/press-releases/eQTL-catalogue
}

\examples{
\dontrun{
run_coloc <- function(eqtl_sumstats, gwas_sumstats)
{
  eQTL_dataset = list(beta = eqtl_sumstats$beta,
                      varbeta = eqtl_sumstats$se^2,
                      N = (eqtl_sumstats$an)[1]/2,
                      MAF = eqtl_sumstats$maf, 
                      type = "quant", 
                      snp = eqtl_sumstats$id)
  gwas_dataset = list(beta = gwas_sumstats$ES,
                      varbeta = gwas_sumstats$SE^2, 
                      type = "quant", 
                      snp = gwas_sumstats$id,
                      MAF = gwas_sumstats$MAF, 
                      N = gwas_sumstats$SS)
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}

lapply(c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer"), require, character.only = T)
tfp <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths = read.delim(tfp, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
tfpi <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_imported.tsv")
imported_tabix_paths = read.delim(tfpi, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# MPV association at the ARHGEF3 locus
region = "3:56615721-57015721"
platelet_df = dplyr::filter(tabix_paths, study == "CEDAR", tissue_label == "platelet")
column_names = colnames(readr::read_tsv("CEDAR_microarray_platelet.all.tsv.gz", n_max = 1))
summary_stats = import_eQTLCatalogue("CEDAR_microarray_platelet.all.tsv.gz", region, selected_gene_id = "ENSG00000163947", column_names)
summary_stats
ggplot(summary_stats, aes(x = position, y = -log(pvalue, 10))) + geom_point()
# GWAS sumstat from the same region
gwas_stats = gwasvcf::query_gwas("ebi-a-GCST004599.vcf.gz", chrompos = "3:56649749-57049749")
gwas_stats = gwasvcf::vcf_to_granges(gwas_stats) %>% keepSeqlevels("3") %>% renameSeqlevels("chr3")
chain = rtracklayer::import.chain("hg19ToHg38.over.chain")
gwas_stats_hg38 = rtracklayer::liftOver(gwas_stats, chain) %>% 
  unlist() %>% 
  renameSeqlevels("3") %>%
  dplyr::as_tibble() %>%
  dplyr::transmute(chromosome = seqnames, position = start, AF, ES, SE, LP, SS) %>%
  dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
  dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
  dplyr::group_by(id) %>% #Keep bi-alleilic variants
  dplyr::mutate(row_count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(row_count == 1) 
ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
# Colocalisation
res = run_coloc(summary_stats, gwas_stats_hg38)

# a. all other eQTL datasets
microarray_df = dplyr::filter(tabix_paths, quant_method == "microarray") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list = setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id)
column_names = colnames(readr::read_tsv("CEDAR_microarray_T-cell_CD8.all.tsv.gz", n_max = 1))
summary_list = purrr::map("~/rds/public_databases/GTEx/csv/Adipose_Subcutaneous.tsv.gz", ~import_eQTLCatalogue(., region,
                          selected_gene_id = "ENSG00000163947", column_names))
coloc_df_microarray = purrr::map_df(summary_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

# b. Uniformly processed RNA-seq datasets
rnaseq_df = dplyr::filter(tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list = setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
column_names = colnames(readr::read_tsv("Alasoo_2018_ge_macrophage_naive.all.tsv.gz", n_max = 1))
safe_import = purrr::safely(import_eQTLCatalogue)
summary_list = purrr::map("Alasoo_2018_ge_macrophage_naive.all.tsv.gz", ~safe_import(., region, selected_gene_id = "ENSG00000163947", column_names))
result_list = purrr::map(summary_list, ~.$result)
result_list = result_list[!unlist(purrr::map(result_list, is.null))]
coloc_df_rnaseq = purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

# c. GTEx_v8 imported eQTL datasets
rnaseq_df = dplyr::filter(imported_tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list = setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
column_names = colnames(readr::read_tsv("~/rds/public_databases/GTEx/csv/Adipose_Subcutaneous.tsv.gz", n_max = 1))
safe_import = purrr::safely(import_eQTLCatalogue)
summary_list = purrr::map("~/rds/public_databases/GTEx/csv/Adipose_Subcutaneous.tsv.gz", ~safe_import(., region,
                          selected_gene_id = "ENSG00000163947", column_names))
result_list = purrr::map(summary_list, ~.$result)
result_list = result_list[!unlist(purrr::map(result_list, is.null))]
result_filtered = purrr::map(result_list, ~dplyr::filter(., !is.na(se)))
coloc_df_imported = purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

coloc_df = dplyr::bind_rows(coloc_df_microarray, coloc_df_rnaseq, coloc_df_imported)
dplyr::arrange(coloc_df, -PP.H4.abf)
ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
}
}
\note{Adapted function}
\keyword{utilities}
