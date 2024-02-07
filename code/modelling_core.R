library(INLA)
library(brinla)
analysis_list = readRDS("data/derived_data/analysis_list.rds")



prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"



# Run suite of models


#Random intercept (different intercept at multiple levels) - a model which assumes trends are similar across space/species

cmb_df = NULL
cmb_df2 = NULL
for(a in names(analysis_list)[c(1:10)]){ 
  for(b in c(1:10)){
  print(a)
  print(b)
  message("   Model 1 - Intercept model")
  if(b == 1){
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
    m1_coef_mn = m1$summary.fixed[2,4]
    m1_coef_lc5 = inla.hpdmarginal(0.5, m1$marginals.fixed[[2]])[1]
    m1_coef_uc5 = inla.hpdmarginal(0.5, m1$marginals.fixed[[2]])[2]
    
    
    
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
    m2_coef_mn = m2$summary.fixed[2,4]
    m2_coef_lc5 = inla.hpdmarginal(0.5, m2$marginals.fixed[[2]])[1]
    m2_coef_uc5 = inla.hpdmarginal(0.5, m2$marginals.fixed[[2]])[2]
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
    m3_coef_mn = m3$summary.fixed[2,4]
    m3_coef_lc5 = inla.hpdmarginal(0.5, m3$marginals.fixed[[2]])[1]
    m3_coef_uc5 = inla.hpdmarginal(0.5, m3$marginals.fixed[[2]])[2]
  } else {
    m1 = inla.rerun(m1)
    m1_coef_mn = m1$summary.fixed[2,4]
    m1_coef_lc5 = inla.hpdmarginal(0.5, m1$marginals.fixed[[2]])[1]
    m1_coef_uc5 = inla.hpdmarginal(0.5, m1$marginals.fixed[[2]])[2]
    
    
    
    #Random slope (different slopes at multiple levels) - a model which captures dissilmar trends across space/species, but assumes independence instead of covariance
    
    
    message("   Model 2 - Slope model")
    m2 = inla.rerun(m2)
    m2_coef_mn = m2$summary.fixed[2,4]
    m2_coef_lc5 = inla.hpdmarginal(0.5, m2$marginals.fixed[[2]])[1]
    m2_coef_uc5 = inla.hpdmarginal(0.5, m2$marginals.fixed[[2]])[2]
    m2_sig = bri.hyperpar.summary(m2)[1,4]
    
    ar_prior = list(theta1 = list(prior="pc.prec", param=c(m2_sig*3, 0.01)),
                    theta2 = list(prior="pc.cor1", param=c(0, 0.9), initial = 0))
    
    #Correlated effect - a model which captures covariances and similarties in trends between neighbouring sites and species
    message("   Model 3 - Correlation model")
    m3 = inla.rerun(m3)
    m3_coef_mn = m3$summary.fixed[2,4]
    m3_coef_lc5 = inla.hpdmarginal(0.5, m3$marginals.fixed[[2]])[1]
    m3_coef_uc5 = inla.hpdmarginal(0.5, m3$marginals.fixed[[2]])[2]
  }
  
  if(b==10){
      m1_coef_mn = m1$summary.fixed[2,4]
      m1_coef_sd = m1$summary.fixed[2,2]
      m1_coef_lc5 = inla.hpdmarginal(0.5, m1$marginals.fixed[[2]])[1]
      m1_coef_uc5 = inla.hpdmarginal(0.5, m1$marginals.fixed[[2]])[2]
      m1_coef_lc = m1$summary.fixed[2,3]
      m1_coef_uc = m1$summary.fixed[2,5]
      m1_fit = m1$summary.fitted.values[,4]
      m1_res = analysis_list[[a]][[1]]$log_abundance - m1_fit
      m1_pred = sd(m1_fit)^2
      m1_sig = bri.hyperpar.summary(m1)[1,4]
      m1_obv = bri.hyperpar.summary(m1)[2,4]
      m1_tip = bri.hyperpar.summary(m1)[3,4]
      m1_gen = bri.hyperpar.summary(m1)[4,4]
      m1_sit = bri.hyperpar.summary(m1)[5,4]
      m1_squ = bri.hyperpar.summary(m1)[6,4]
      
      m2_coef_mn = m2$summary.fixed[2,4]
      m2_coef_sd = m2$summary.fixed[2,2]
      m2_coef_lc5 = inla.hpdmarginal(0.5, m2$marginals.fixed[[2]])[1]
      m2_coef_uc5 = inla.hpdmarginal(0.5, m2$marginals.fixed[[2]])[2]
      m2_coef_lc = m2$summary.fixed[2,3]
      m2_coef_uc = m2$summary.fixed[2,5]
      m2_fit = m2$summary.fitted.values[,4]
      m2_res = analysis_list[[a]][[1]]$log_abundance - m2_fit
      m2_pred = sd(m2_fit)^2
      m2_sig = bri.hyperpar.summary(m2)[1,4]
      m2_obv = bri.hyperpar.summary(m2)[2,4]
      m2_tip = bri.hyperpar.summary(m2)[3,4]
      m2_gen = bri.hyperpar.summary(m2)[4,4]
      m2_sit = bri.hyperpar.summary(m2)[5,4]
      m2_squ = bri.hyperpar.summary(m2)[6,4]
      
      m3_coef_mn = m3$summary.fixed[2,4]
      m3_coef_sd = m3$summary.fixed[2,2]
      m3_coef_lc5 = inla.hpdmarginal(0.5, m3$marginals.fixed[[2]])[1]
      m3_coef_uc5 = inla.hpdmarginal(0.5, m3$marginals.fixed[[2]])[2]
      m3_coef_lc = m3$summary.fixed[2,3]
      m3_coef_uc = m3$summary.fixed[2,5]
      m3_fit = m3$summary.fitted.values[,4]
      m3_res = analysis_list[[a]][[1]]$cent_abundance - m3_fit
      m3_pred = sd(m3_fit)^2
      m3_sig = bri.hyperpar.summary(m3)[1,4]
      m3_obv_auto = bri.hyperpar.summary(m3)[2,4]
      m3_phi = bri.hyperpar.summary(m3)[3,4]
      m3_obv = bri.hyperpar.summary(m3)[4,4]
      m3_tip_h = bri.hyperpar.summary(m3)[5,4]
      m3_tip = bri.hyperpar.summary(m3)[6,4]
      m3_gen = bri.hyperpar.summary(m3)[7,4]
      m3_sit_h = bri.hyperpar.summary(m3)[8,4]
      m3_sit = bri.hyperpar.summary(m3)[9,4]
      m3_squ = bri.hyperpar.summary(m3)[10,4]
      
      tmp_df2 = data.frame(
        code = a,
        model = c(1:3),
        coef = c(m1_coef_mn,m2_coef_mn,m3_coef_mn),
        coef_sd = c(m1_coef_sd,m2_coef_sd,m3_coef_sd),
        coef_lc5 = c(m1_coef_lc5,m2_coef_lc5,m3_coef_lc5),
        coef_uc5 = c(m1_coef_uc5,m2_coef_uc5,m3_coef_uc5),
        coef_lc = c(m1_coef_lc,m2_coef_lc,m3_coef_lc),
        coef_uc = c(m1_coef_uc,m2_coef_uc,m3_coef_uc),
        fix = c(m1_pred,m2_pred,m3_pred),
        obv = c(m1_obv,m2_obv,m3_obv),
        obv_auto = c(NA,NA,m3_obv_auto),
        tip = c(m1_tip,m2_tip,m3_tip),
        gen = c(m1_gen,m2_gen,m3_gen),
        tip_h = c(NA,NA,m3_tip_h),
        sit = c(m1_sit,m2_sit,m3_sit),
        squ = c(m1_squ,m2_squ,m3_squ),
        sit_h = c(NA,NA,m3_sit_h),
        sig = c(m1_sig,m2_sig,m3_sig),
        phi = c(NA,NA,m3_phi)
      )
      
      cmb_df2 = rbind(cmb_df2, tmp_df2)
  } else {}

  tmp_df = data.frame(
    dataset = a,
    run = b,
    model = c(1,2,3),
    coef = c(m1_coef_mn, m2_coef_mn, m3_coef_mn),
    low = c(m1_coef_lc5, m2_coef_lc5, m3_coef_lc5),
    upp = c(m1_coef_uc5, m2_coef_uc5, m3_coef_uc5),
    status = c(m1$mode$mode.status, m3$mode$mode.status, m3$mode$mode.status))
  print(tmp_df)
  cmb_df = rbind(cmb_df, tmp_df)
  results = list(m1,m2,m3)
  saveRDS(results, paste0("outputs/model_output_convergence",b,"_",a,".rds"))
  }
  rm(m1,m2,m3)
}
write.csv(cmb_df, "outputs/model_summary_convergence.csv")
write.csv(cmb_df2, "outputs/model_summary.csv")
