library(INLA)
library(brinla)
library(data.table)
analysis_list = readRDS("data/derived_data/analysis_list.rds")



prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"

# Run suite of models


#Random intercept (different intercept at multiple levels) - a model which assumes trends are similar across space/species

for(a in names(analysis_list)[1:10]){ 
  
  tmp_df = analysis_list[[a]][[1]]
  tmp_df$log_abundance_trim = tmp_df$log_abundance
  tmp_df$cent_abundance_trim = tmp_df$cent_abundance
  
  uni_sites = unique(tmp_df$site_spec_code)
  uni_sites_sample = sample(uni_sites,(length(uni_sites)/2))
  remove_abundances = list()
  for(c in c(1:(length(uni_sites)))){
    trim_tmp_df = subset(tmp_df, site_spec_code == uni_sites[c])
    if(uni_sites[c] %in% uni_sites_sample){
      trim_tmp_df$log_abundance_trim[c(nrow(trim_tmp_df))]  = NA
      trim_tmp_df$cent_abundance_trim[c(nrow(trim_tmp_df))]  = NA
    } else {
      
    }
    remove_abundances[[c]] = trim_tmp_df
  }
  remove_abundances = rbindlist(remove_abundances)
  analysis_list[[a]][[1]] = remove_abundances
  
  print(a)
  message("   Model 1 - Intercept model")
  m1 = inla(log_abundance_trim ~
              year_centre +
              f(site_spec_code,  model = "iid", constr = F, hyper = prior_prec) +
              f(tips_code, model = "iid", constr = F, hyper = prior_prec) +
              f(genus_code, model = "iid", constr = F, hyper = prior_prec) +
              f(site_code, model = "iid", constr = F, hyper = prior_prec) +
              f(region_code, model = "iid", constr = F, hyper = prior_prec),
            data = analysis_list[[a]][[1]],  family = "gaussian", 
            control.predictor=list(compute=TRUE),
            control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                 mean = 0, prec = 1),
            num.threads = "1:1")
    m1 = inla.rerun(m1)
    m1 = inla.rerun(m1)
  
  #Random slope (different slopes at multiple levels) - a model which captures dissilmar trends across space/species, but assumes independence instead of covariance
  
  
  message("   Model 2 - Slope model")
  m2 = inla(cent_abundance_trim ~ 
              year_centre +
              f(site_spec_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(tips_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(genus_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(site_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(region_code, year_centre, model = "iid", constr = F, hyper = prior_prec),
            data = analysis_list[[a]][[1]],  family = "gaussian", 
            control.predictor=list(compute=TRUE),
            control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                 mean = 0, prec = 1),
            num.threads = "1:1")
  m2 = inla.rerun(m2)
  m2 = inla.rerun(m2)
  m2_sig = bri.hyperpar.summary(m2)[1,4]

  
  ar_prior = list(theta1 = list(prior="pc.prec", param=c(m2_sig*3, 0.01)),
                  theta2 = list(prior="pc.cor1", param=c(0, 0.9), initial = 0))
  
  
  #Correlated effect - a model which captures covariances and similarties in trends between neighbouring sites and species
  message("   Model 3 - Correlation model")
  m3 = inla(cent_abundance_trim ~ 
              year_centre +
              f(year3, model = "ar1", replicate = site_spec_code2, hyper = ar_prior) +
              f(site_spec_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(tips_code2, year_centre, model = "generic0", 
                constr = F, Cmatrix = analysis_list[[a]][[2]], hyper = prior_prec) +
              f(tips_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(genus_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(site_code2, year_centre, model = "generic0", 
                constr = F, Cmatrix = analysis_list[[a]][[3]], hyper = prior_prec) +
              f(site_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(region_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec),
            data = analysis_list[[a]][[1]], family = "gaussian", 
            control.predictor=list(compute=TRUE),
            control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                 mean = 0, prec = 1),
            num.threads = "1:1")
  m3 = inla.rerun(m3)
  m3 = inla.rerun(m3)

  
  results = list(m1,m2,m3,remove_abundances)
  saveRDS(results, paste0("outputs/model_output_miss_abundance_",a,".rds"))
  rm(m1,m2,m3,remove_abundances)
}



