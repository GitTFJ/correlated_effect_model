library(INLA)
library(brinla)
analysis_list = readRDS("data/derived_data/analysis_list_predict2.rds")



prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"


# Run suite of models


#Random intercept (different intercept at multiple levels) - a model which assumes trends are similar across space/species

cmb_df = NULL
for(a in names(analysis_list)){ 
  print(a)
  
  m1 = inla(log_abundance ~
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
  
  #Random slope (different slopes at multiple levels) - a model which captures dissilmar trends across space/species, but assumes independence instead of covariance
  
  
  message("   Model 2 - Slope model")
  m2 = inla(cent_abundance ~ 
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
  
  m2_sig = bri.hyperpar.summary(m2)[1,4]
  ar_prior = list(theta1 = list(prior="pc.prec", param=c(m2_sig*3, 0.01)),
                  theta2 = list(prior="pc.cor1", param=c(0, 0.9), initial = 0))
  #Correlated effect - a model which captures covariances and similarties in trends between neighbouring sites and species
  message("   Model 3 - Correlation model")
  m3 = inla(cent_abundance ~ 
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
  
  
  
  results = list(m1,m2,m3)
  saveRDS(results, paste0("outputs/model_output_predict2_",a,".rds"))
  #rm(m1,m2,m3)
}
summary(m3)
