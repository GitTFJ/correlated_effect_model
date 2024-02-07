library(INLA)
library(brinla)
analysis_list = readRDS("data/derived_data/analysis_list_phylo.rds")




prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"


# Run suite of models


#Random intercept (different intercept at multiple levels) - a model which assumes trends are similar across space/species

cmb_df = NULL
for(a in names(analysis_list)[c(1:10)]){
  print(a)
  
  
  m2 = inla(cent_abundance ~ 
              year_centre +
              f(site_spec_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(tips_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
              f(genus_code1, year_centre, model = "iid", constr = F, hyper = prior_prec) +
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
  
  message("   Model 3 - Correlation model with the trimmed rotl data")
  
  m3a = inla(cent_abundance ~ 
               year_centre +
               f(year3, model = "ar1", replicate = site_spec_code2, hyper = ar_prior) +
               f(site_spec_code, year_centre, model = "iid", 
                 constr = F, hyper = prior_prec) +
               f(tips_code2, year_centre, model = "generic0",
                 constr = F, Cmatrix = analysis_list[[a]][[2]], hyper = prior_prec) +
               f(tips_code, year_centre, model = "iid", 
                 constr = F, hyper = prior_prec) +
               f(genus_code1, year_centre, model = "iid",
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
  m3a = inla.rerun(m3a)
  m3a = inla.rerun(m3a)
  
  m3a_coef_mn = m3a$summary.fixed[2,4]
  m3a_coef_sd = m3a$summary.fixed[2,2]
  m3a_coef_lc5 = inla.hpdmarginal(0.5, m3a$marginals.fixed[[2]])[1]
  m3a_coef_uc5 = inla.hpdmarginal(0.5, m3a$marginals.fixed[[2]])[2]
  m3a_coef_lc = m3a$summary.fixed[2,3]
  m3a_coef_uc = m3a$summary.fixed[2,5]
  m3a_fit = m3a$summary.fitted.values[,4]
  m3a_res = analysis_list[[a]][[1]]$cent_abundance - m3a_fit
  m3a_pred = sd(m3a_fit)^2
  m3a_sig = bri.hyperpar.summary(m3a)[1,4]
  m3a_obv_auto = bri.hyperpar.summary(m3a)[2,4]
  m3a_phi = bri.hyperpar.summary(m3a)[3,4]
  m3a_obv = bri.hyperpar.summary(m3a)[4,4]
  m3a_tip_h = bri.hyperpar.summary(m3a)[5,4]
  m3a_tip = bri.hyperpar.summary(m3a)[6,4]
  m3a_gen = bri.hyperpar.summary(m3a)[7,4]
  m3a_sit_h = bri.hyperpar.summary(m3a)[8,4]
  m3a_sit = bri.hyperpar.summary(m3a)[9,4]
  m3a_squ = bri.hyperpar.summary(m3a)[10,4]
  
  message("   Model 3 - Correlation model with the timetree data")

  m3b = inla(cent_abundance ~ 
               year_centre +
               f(year3, model = "ar1", replicate = site_spec_code2, hyper = ar_prior) +
               f(site_spec_code, year_centre, model = "iid", 
                 constr = F, hyper = prior_prec) +
               f(tips_code2, year_centre, model = "generic0", 
                 constr = F, Cmatrix = analysis_list[[a]][[4]], hyper = prior_prec) +
               f(tips_code, year_centre, model = "iid", 
                 constr = F, hyper = prior_prec) +
               f(genus_code2, year_centre, model = "iid", 
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
  m3b = inla.rerun(m3b)
  m3b = inla.rerun(m3b)
  
  m3b_coef_mn = m3b$summary.fixed[2,4]
  m3b_coef_sd = m3b$summary.fixed[2,2]
  m3b_coef_lc5 = inla.hpdmarginal(0.5, m3b$marginals.fixed[[2]])[1]
  m3b_coef_uc5 = inla.hpdmarginal(0.5, m3b$marginals.fixed[[2]])[2]
  m3b_coef_lc = m3b$summary.fixed[2,3]
  m3b_coef_uc = m3b$summary.fixed[2,5]
  m3b_fit = m3b$summary.fitted.values[,4]
  m3b_res = analysis_list[[a]][[1]]$cent_abundance - m3b_fit
  m3b_pred = sd(m3b_fit)^2
  m3b_sig = bri.hyperpar.summary(m3b)[1,4]
  m3b_obv_auto = bri.hyperpar.summary(m3b)[2,4]
  m3b_phi = bri.hyperpar.summary(m3b)[3,4]
  m3b_obv = bri.hyperpar.summary(m3b)[4,4]
  m3b_tip_h = bri.hyperpar.summary(m3b)[5,4]
  m3b_tip = bri.hyperpar.summary(m3b)[6,4]
  m3b_gen = bri.hyperpar.summary(m3b)[7,4]
  m3b_sit_h = bri.hyperpar.summary(m3b)[8,4]
  m3b_sit = bri.hyperpar.summary(m3b)[9,4]
  m3b_squ = bri.hyperpar.summary(m3b)[10,4]
  
  
  
  tmp_df = data.frame(
    code = a,
    model = c("3_rotl_trim", "3_timetree"),
    coef = c(m3a_coef_mn, m3b_coef_mn),
    coef_sd = c(m3a_coef_sd, m3b_coef_sd),
    coef_lc5 = c(m3a_coef_lc5, m3b_coef_lc5),
    coef_uc5 = c(m3a_coef_uc5, m3b_coef_uc5),
    coef_lc = c(m3a_coef_lc, m3b_coef_lc),
    coef_uc = c(m3a_coef_uc, m3b_coef_uc),
    fix = c(m3a_pred, m3b_pred),
    obv = c(m3a_obv, m3b_obv),
    obv_auto = c(m3a_obv_auto, m3b_obv_auto),
    tip = c(m3a_tip, m3b_tip),
    gen = c(m3a_gen, m3b_gen),
    tip_h = c(m3a_tip_h, m3b_tip_h),
    sit = c(m3a_sit, m3b_sit),
    squ = c(m3a_squ, m3b_squ),
    sit_h = c(m3a_sit_h, m3b_sit_h),
    sig = c(m3a_sig, m3b_sig),
    phi = c(m3a_phi, m3b_phi)
  )  
  
  cmb_df = rbind(cmb_df, tmp_df)
  results = list(m3a, m3b)
  saveRDS(results, paste0("outputs/model_output_sensitivity_to_phylo",a,".rds"))
}
write.csv(cmb_df, "outputs/model_summary_phylo.csv")

