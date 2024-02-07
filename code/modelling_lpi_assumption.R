library(INLA)
library(brinla)



prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"

  m1 = inla(log_abundance ~
              year_centre +
              f(site_spec_code,  model = "iid", constr = F, hyper = prior_prec) +
              f(tips_code, model = "iid", constr = F, hyper = prior_prec) +
              f(genus_code, model = "iid", constr = F, hyper = prior_prec) +
              f(site_code, model = "iid", constr = F, hyper = prior_prec) +
              f(region_code, model = "iid", constr = F, hyper = prior_prec),
            data = analysis_list_lpi[[1]],  family = "gaussian", 
            control.predictor=list(compute=TRUE),
            control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                 mean = 0, prec = 1),
            num.threads = 4)

  
  
  #Random slope (different slopes at multiple levels) - a model which captures dissilmar trends across space/species, but assumes independence instead of covariance
  
  
  message("   Model 2 - Slope model")
  m2 = inla(cent_abundance ~ 
              year_centre +
              f(site_spec_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(tips_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(genus_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(site_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(region_code, year_centre, model = "iid", constr = F, hyper = prior_prec),
            data = analysis_list_lpi[[1]],  family = "gaussian", 
            control.predictor=list(compute=TRUE),
            control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                 mean = 0, prec = 1),
            num.threads = 4)
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
                constr = F, Cmatrix = analysis_list_lpi[[2]], hyper = prior_prec) +
              f(tips_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(genus_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(site_code2, year_centre, model = "generic0", 
                constr = F, Cmatrix = analysis_list_lpi[[3]], hyper = prior_prec) +
              f(site_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec) +
              f(region_code, year_centre, model = "iid", 
                constr = F, hyper = prior_prec),
            data = analysis_list_lpi[[1]], family = "gaussian", 
            control.predictor=list(compute=TRUE),
            control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                 mean = 0, prec = 1),
            num.threads = 4)

  results = list(m1,m2,m3)
  saveRDS(results, "outputs/model_output_lpi_assumption.rds")
