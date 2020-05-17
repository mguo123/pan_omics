
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

######################### MAPPING SYMBOLS ##########################
mapSym = function (rna_path, save_prefix){
  rna_df = read.csv(rna_path,row.names=1,check.names = FALSE)
  rna_df$symbol <- mapIds(org.Hs.eg.db, keys=row.names(rna_df), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  rna_df_filt = rna_df[!is.na(rna_df$symbol),]
  rna_df_filt = rna_df_filt %>%
      group_by(symbol) %>%
      summarise_all(funs(mean(., na.rm=TRUE)))%>%
      column_to_rownames("symbol")

  rna_df_filt_long = rna_df_filt%>%
      rownames_to_column()%>%
      melt(id.vars = c("rowname"))

  rna_df_filt_long_thres = filter(rna_df_filt_long, value > THRES)
  rna_df_thres = rna_df_filt_long_thres%>%
      dcast(rowname~variable)%>%
      column_to_rownames("rowname")

  rna_df_thres_bin = rna_df_thres%>%
      replace(!is.na(.), 1)%>%
      replace(is.na(.), 0)

  write.csv(rna_df_filt, paste0(save_prefix,'_sym.csv'))

  write.csv(rna_df_thres, paste0(save_prefix,'_sym_thres.csv'))

  write.csv(rna_df_thres_bin, paste0(save_prefix,'_sym_bin.csv'))
}

# tpm_df = read.csv(paste0(dirname, '/tissue_tpm.csv'),row.names=1,check.names = FALSE)
# fpkm_df = read.csv(paste0(dirname, '/tissue_fpkm.csv'),row.names=1,check.names = FALSE)
# tpm_df_sample = read.csv(paste0(dirname, '/sample_tpm.csv'),row.names=1,check.names = FALSE)
# fpkm_df_sample = read.csv(paste0(dirname, '/sample_fpkm.csv'),row.names=1,check.names = FALSE)

print('mapping tissue tpm')
mapSym(paste0(dirname, '/tissue_tpm.csv'), paste0(dirname, '/tissue_tpm'))
print('mapping tissue fpkm')
mapSym(paste0(dirname, '/tissue_fpkm.csv'), paste0(dirname, '/tissue_fpkm'))
print('mapping sample tpm')
mapSym(paste0(dirname, '/sample_tpm.csv'), paste0(dirname, '/sample_tpm'))
print('mapping sample fpkm')
mapSym(paste0(dirname, '/sample_fpkm.csv'), paste0(dirname, '/sample_fpkm'))

# tpm_df$symbol <- mapIds(org.Hs.eg.db, keys=row.names(tpm_df), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
# fpkm_df$symbol <- mapIds(org.Hs.eg.db, keys=row.names(fpkm_df), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

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
