
###process_rna.R

###Make TPM/FPKM matrix with symbols

###4/09/2020 - from `make_rna_mat_R.ipynb` 12/19/2019

library(tidyverse)
library(reshape2)
library( org.Hs.eg.db )
library(AnnotationDbi)

THRES=1

######################### READ FILES ##########################
# test if there is at least one argument: if not, return an error
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
# default output file
dirname = args[1]
print(dirname)

tpm_df = read.csv(paste0(dirname, '/tissue_tpm.csv'),row.names=1)
fpkm_df = read.csv(paste0(dirname, '/tissue_fpkm.csv'),row.names=1)

# ######################### MAPPING SYMBOLS ##########################
#
# tpm_df$symbol <- mapIds(org.Hs.eg.db, keys=row.names(tpm_df), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
# fpkm_df$symbol <- mapIds(org.Hs.eg.db, keys=row.names(fpkm_df), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
#
# tpm_df_filt = tpm_df[!is.na(tpm_df$symbol),]
# tpm_df_filt = tpm_df_filt %>%
#     group_by(symbol) %>%
#     summarise_all(funs(mean(., na.rm=TRUE)))%>%
#     column_to_rownames("symbol")
#
# fpkm_df_filt = fpkm_df[!is.na(fpkm_df$symbol),]
# fpkm_df_filt = fpkm_df_filt %>%
#     group_by(symbol) %>%
#     summarise_all(funs(mean(., na.rm=TRUE)))%>%
#     column_to_rownames("symbol")
#
# write.csv(tpm_df_filt, paste0(dirname, '/tissue_tpm_sym.csv'))
# write.csv(fpkm_df_filt, paste0(dirname, '/tissue_fpkm_sym.csv'))
#
# ######################### MAPPING SYMBOLS ##########################
#
# tpm_df_filt_long = tpm_df_filt%>%
#     rownames_to_column()%>%
#     melt(id.vars = c("rowname"))
#
# tpm_df_filt_long_thres = filter(tpm_df_filt_long, value > THRES)
# tpm_df_thres = tpm_df_filt_long_thres%>%
#     dcast(rowname~variable)%>%
#     column_to_rownames("rowname")
#
# tpm_df_thres_bin = tpm_df_thres%>%
#     replace(!is.na(.), 1)%>%
#     replace(is.na(.), 0)
#
# fpkm_df_filt_long = fpkm_df_filt%>%
#     rownames_to_column()%>%
#     melt(id.vars = c("rowname"))
#
# fpkm_df_filt_long_thres = filter(fpkm_df_filt_long, value > THRES)
# fpkm_df_thres = fpkm_df_filt_long_thres%>%
#     dcast(rowname~variable)%>%
#     column_to_rownames("rowname")
#
# fpkm_df_thres_bin = fpkm_df_thres %>%
#     replace(!is.na(.), 1)%>%
#     replace(is.na(.), 0)
#
# write.csv(tpm_df_filt, paste0(dirname, '/tissue_tpm_sym.csv'))
# write.csv(fpkm_df_filt, paste0(dirname, '/tissue_fpkm_sym.csv'))
#
# write.csv(tpm_df_thres, paste0(dirname, '/tissue_tpm_sym_thres.csv'))
# write.csv(fpkm_df_thres, paste0(dirname, '/tissue_fpkm_sym_thres.csv'))
#
# write.csv(tpm_df_thres_bin, paste0(dirname, '/tissue_tpm_sym_bin.csv'))
# write.csv(fpkm_df_thres_bin, paste0(dirname, '/tissue_fpkm_sym_bin.csv'))
