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
for(a in names(analysis_list)[1:10]){
  print(a)
  
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
            num.threads = 4)
  m2_sig = bri.hyperpar.summary(m2)[1,4]
  ar_prior = list(theta1 = list(prior="pc.prec", param=c(m2_sig*3, 0.01)),
                  theta2 = list(prior="pc.cor1", param=c(0, 0.9), initial = 0))
  
  print("   Model 3 - Correlation model without time")
  m3_a = inla(cent_abundance ~ 
                year_centre +
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
              data = analysis_list[[a]][[1]], family = "gaussian", num.threads = 4)
  m3_a_coef_mn = m3_a$summary.fixed[2,4]
  m3_a_coef_sd = m3_a$summary.fixed[2,2]
  m3_a_coef_lc = m3_a$summary.fixed[2,3]
  m3_a_coef_uc = m3_a$summary.fixed[2,5]
  
  print("   Model 3 - Correlation model without space")
  m3_b = inla(cent_abundance ~ 
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
                f(site_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(region_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec),
              data = analysis_list[[a]][[1]], family = "gaussian", num.threads = 4)
  m3_b_coef_mn = m3_b$summary.fixed[2,4]
  m3_b_coef_sd = m3_b$summary.fixed[2,2]
  m3_b_coef_lc = m3_b$summary.fixed[2,3]
  m3_b_coef_uc = m3_b$summary.fixed[2,5]
  
  print("   Model 3 - Correlation model without phylogeny")
  m3_c = inla(cent_abundance ~ 
                year_centre +
                f(year3, model = "ar1", replicate = site_spec_code2, hyper = ar_prior) +
                f(site_spec_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
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
              data = analysis_list[[a]][[1]], family = "gaussian", num.threads = 4)
  m3_c_coef_mn = m3_c$summary.fixed[2,4]
  m3_c_coef_sd = m3_c$summary.fixed[2,2]
  m3_c_coef_lc = m3_c$summary.fixed[2,3]
  m3_c_coef_uc = m3_c$summary.fixed[2,5]
  
  
  print("   Model 3 - Correlation model without time and space")
  m3_d = inla(cent_abundance ~ 
                year_centre +
                f(site_spec_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(tips_code2, year_centre, model = "generic0",
                  constr = F, Cmatrix = analysis_list[[a]][[2]], hyper = prior_prec) +
                f(tips_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(genus_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(site_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(region_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec),
              data = analysis_list[[a]][[1]], family = "gaussian", num.threads = 4)
  m3_d_coef_mn = m3_d$summary.fixed[2,4]
  m3_d_coef_sd = m3_d$summary.fixed[2,2]
  m3_d_coef_lc = m3_d$summary.fixed[2,3]
  m3_d_coef_uc = m3_d$summary.fixed[2,5]
  
  print("   Model 3 - Correlation model without time and phylogeny")
  m3_e = inla(cent_abundance ~ 
                year_centre +
                f(site_spec_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
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
              data = analysis_list[[a]][[1]], family = "gaussian", num.threads = 4)
  m3_e_coef_mn = m3_e$summary.fixed[2,4]
  m3_e_coef_sd = m3_e$summary.fixed[2,2]
  m3_e_coef_lc = m3_e$summary.fixed[2,3]
  m3_e_coef_uc = m3_e$summary.fixed[2,5]
  
  print("   Model 3 - Correlation model without space and phylogeny")
  m3_f = inla(cent_abundance ~ 
                year_centre +
                f(year3, model = "ar1", replicate = site_spec_code2, hyper = ar_prior) +
                f(site_spec_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(tips_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(genus_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(site_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec) +
                f(region_code, year_centre, model = "iid", 
                  constr = F, hyper = prior_prec),
              data = analysis_list[[a]][[1]], family = "gaussian", num.threads = 4)
  m3_f_coef_mn = m3_f$summary.fixed[2,4]
  m3_f_coef_sd = m3_f$summary.fixed[2,2]
  m3_f_coef_lc = m3_f$summary.fixed[2,3]
  m3_f_coef_uc = m3_f$summary.fixed[2,5]
  
  tmp_df = data.frame(
    code = a,
    model_type = c(
      "Time absent", 
      "Space absent", 
      "Phylogeny absent",
      "Time and space absent",
      "Time and phylogeny absent",
      "Space and phylogeny absent"),
    coef = c(m3_a_coef_mn, m3_b_coef_mn, m3_c_coef_mn, m3_d_coef_mn, m3_e_coef_mn, m3_f_coef_mn),
    coef_sd = c(m3_a_coef_sd, m3_b_coef_sd, m3_c_coef_sd, m3_d_coef_sd, m3_e_coef_sd, m3_f_coef_sd),
    coef_lc = c(m3_a_coef_lc, m3_b_coef_lc, m3_c_coef_lc, m3_d_coef_lc, m3_e_coef_lc, m3_f_coef_lc),
    coef_uc = c(m3_a_coef_uc, m3_b_coef_uc, m3_c_coef_uc, m3_d_coef_uc, m3_e_coef_uc, m3_f_coef_uc)
  )
  
  
  cmb_df = rbind(cmb_df, tmp_df)
  results = list(m3_a,m3_b,m3_c,m3_d,m3_e,m3_f)
  saveRDS(results, paste0("outputs/model_output_structure_",a,".rds"))
  
}
write.csv(cmb_df, "outputs/model_summary_structure.csv")
