# Makes reference/norms_cnorm.rds, the cNORM values engine.R compares with when cNORM is not
# installed. Run from inst/validation with cNORM installed:  Rscript reference/make_norms_cnorm.R
source("norms_sim.R")
saveRDS(cnorm_ref(norm_sim()), "reference/norms_cnorm.rds", version = 2)
