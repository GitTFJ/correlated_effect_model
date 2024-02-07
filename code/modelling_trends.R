library(INLA)
library(brinla)
library(data.table)
library(dplyr)
analysis_list = readRDS("data/derived_data/analysis_list.rds")



prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"


# Run suite of models
cmb_df = NULL
for(a in names(analysis_list)[c(1:10)]){ 
  
  m_tmp1 = inla(cent_abundance ~
                  f(site_spec_code, year_centre, model = "iid", constr = F, hyper = prior_prec), data = analysis_list[[a]][[1]])
  sig = bri.hyperpar.summary(m_tmp1)[1,4]
  ar_prior = list(theta1 = list(prior="pc.prec", param=c(sig*3, 0.01)),
                  theta2 = list(prior="pc.cor1", param=c(0, 0.9), initial = 0))
  m_tmp2 = inla(cent_abundance ~
                  f(year3, model = "ar1", replicate = site_spec_code2, hyper = ar_prior) +
                  f(site_spec_code, year_centre, model = "iid", 
                    constr = F, hyper = prior_prec),
                data = analysis_list[[a]][[1]])
  sites = m_tmp2$summary.random$site_spec_code
  sites$sig = ifelse(sites$`0.025quant` > 0 | sites$`0.975quant` < 0, "sig", "ns")
  
  cmb_trends = NULL
  for(b in c(1:50)){
    print(b)
    tmp_df = analysis_list[[a]][[1]]
    tmp_df$log_abundance_trim = tmp_df$log_abundance
    tmp_df$cent_abundance_trim = tmp_df$cent_abundance
    
    uni_sites = sites[which(sites$sig == "sig"),]$ID
    uni_sites_sample = sample(uni_sites,1)
    remove_trends = list()
    for(c in uni_sites){
      trim_tmp_df = subset(tmp_df, site_spec_code == c)
      if(uni_sites[c] %in% uni_sites_sample){
        trim_tmp_df$log_abundance_trim = NA
        trim_tmp_df$cent_abundance_trim = NA
      } else {
        
      }
      remove_trends[[c]] = trim_tmp_df
    }
    remove_trends = rbindlist(remove_trends)
    
    
    #Random slope (different slopes at multiple levels) - a model which captures dissilmar trends across space/species, but assumes independence instead of covariance
    
    
    message("   Model 2 - Slope model")
    m2 = inla(cent_abundance_trim ~ 
                year_centre +
                f(site_spec_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
                f(tips_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
                f(genus_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
                f(site_code, year_centre, model = "iid", constr = F, hyper = prior_prec) +
                f(region_code, year_centre, model = "iid", constr = F, hyper = prior_prec),
              data = remove_trends,  family = "gaussian", 
              control.predictor=list(compute=TRUE),
              control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                   mean = 0, prec = 1),
              num.threads = "1:1")
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
              data = remove_trends, family = "gaussian", 
              control.predictor=list(compute=TRUE),
              control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                   mean = 0, prec = 1),
              num.threads = "1:1")
    
    
    link_df = unique(tmp_df[,c("site_spec_code", "tips_code", "genus_code", "site_code", "region_code")])
    tip_rarity = link_df %>%
      group_by(tips_code) %>%
      summarise(N_tip = n())
    genus_rarity = link_df %>%
      group_by(genus_code) %>%
      summarise(N_gen = n())
    site_rarity = link_df %>%
      group_by(site_code) %>%
      summarise(N_sit = n())
    reg_rarity = link_df %>%
      group_by(region_code) %>%
      summarise(N_reg = n())
    
    link_df = left_join(link_df, tip_rarity)
    link_df = left_join(link_df, genus_rarity)
    link_df = left_join(link_df, site_rarity)
    link_df = left_join(link_df, reg_rarity)
    
    link_df$sum = link_df$N_tip + link_df$N_gen + link_df$N_sit + link_df$N_reg
    
    link_df = link_df[link_df$site_spec_code %in% uni_sites_sample, ]
    
    site_spec_code_med = m2$summary.random$site_spec_code[,c(1,2,4,6)]
    colnames(site_spec_code_med) = c("ID", "m2_site_spec_mn","m2_site_spec_lc", "m2_site_spec_uc")
    link_df = left_join(link_df, site_spec_code_med, by = c("site_spec_code" = "ID"))
    
    tips_code_med = m2$summary.random$tips_code[,c(1,2)]
    colnames(tips_code_med) = c("ID", "m2_tip_mn")
    link_df = left_join(link_df, tips_code_med, by = c("tips_code" = "ID"))
    
    genus_code_med = m2$summary.random$genus_code[,c(1,2)]
    colnames(genus_code_med) = c("ID", "m2_gen_mn")
    link_df = left_join(link_df, genus_code_med, by = c("genus_code" = "ID"))
    
    site_code_med = m2$summary.random$site_code[,c(1,2)]
    colnames(site_code_med) = c("ID", "m2_spec_mn")
    link_df = left_join(link_df, site_code_med, by = c("site_code" = "ID"))
    
    region_code_med = m2$summary.random$region_code[,c(1,2)]
    colnames(region_code_med) = c("ID", "m2_reg_mn")
    link_df = left_join(link_df, region_code_med, by = c("region_code" = "ID"))
    
    link_df$m2_coef = m2$summary.fixed$mean[2]
    
    site_spec_code_med = m3$summary.random$site_spec_code[,c(1,2,4,6)]
    colnames(site_spec_code_med) = c("ID", "m3_site_spec_mn","m3_site_spec_lc", "m3_site_spec_uc")
    link_df = left_join(link_df, site_spec_code_med, by = c("site_spec_code" = "ID"))
    
    tips_code_med = m3$summary.random$tips_code[,c(1,2)]
    colnames(tips_code_med) = c("ID", "m3_tip_mn")
    link_df = left_join(link_df, tips_code_med, by = c("tips_code" = "ID"))
    
    tips_code_med = m3$summary.random$tips_code2[,c(1,2)]
    colnames(tips_code_med) = c("ID", "m3_tip2_mn")
    link_df = left_join(link_df, tips_code_med, by = c("tips_code" = "ID"))
    
    genus_code_med = m3$summary.random$genus_code[,c(1,2)]
    colnames(genus_code_med) = c("ID", "m3_gen_mn")
    link_df = left_join(link_df, genus_code_med, by = c("genus_code" = "ID"))
    
    site_code_med = m3$summary.random$site_code[,c(1,2)]
    colnames(site_code_med) = c("ID", "m3_spec_mn")
    link_df = left_join(link_df, site_code_med, by = c("site_code" = "ID"))
    
    site_code_med = m3$summary.random$site_code2[,c(1,2)]
    colnames(site_code_med) = c("ID", "m3_spec2_mn")
    link_df = left_join(link_df, site_code_med, by = c("site_code" = "ID"))
    
    region_code_med = m3$summary.random$region_code[,c(1,2)]
    colnames(region_code_med) = c("ID", "m3_reg_mn")
    link_df = left_join(link_df, region_code_med, by = c("region_code" = "ID"))
    
    link_df$m3_coef = m3$summary.fixed$mean[2]
    
    link_df$pop_trend2 = 
      link_df$m2_site_spec_mn + 
      link_df$m2_tip_mn +
      link_df$m2_gen_mn +
      link_df$m2_spec_mn +
      link_df$m2_reg_mn +
      link_df$m2_coef
    
    link_df$pop_trend2_lc = 
      link_df$m2_site_spec_lc + 
      link_df$m2_tip_mn +
      link_df$m2_gen_mn +
      link_df$m2_spec_mn +
      link_df$m2_reg_mn +
      link_df$m2_coef
    
    link_df$pop_trend2_uc = 
      link_df$m2_site_spec_uc + 
      link_df$m2_tip_mn +
      link_df$m2_gen_mn +
      link_df$m2_spec_mn +
      link_df$m2_reg_mn +
      link_df$m2_coef
    
    link_df$pop_trend3 = 
      link_df$m3_site_spec_mn + 
      link_df$m3_tip_mn +
      link_df$m3_tip2_mn +
      link_df$m3_gen_mn +
      link_df$m3_spec_mn +
      link_df$m3_spec2_mn +
      link_df$m3_reg_mn +
      link_df$m3_coef
    
    link_df$pop_trend3_lc = 
      link_df$m3_site_spec_lc + 
      link_df$m3_tip_mn +
      link_df$m3_tip2_mn +
      link_df$m3_gen_mn +
      link_df$m3_spec_mn +
      link_df$m3_spec2_mn +
      link_df$m3_reg_mn +
      link_df$m3_coef
    
    link_df$pop_trend3_uc = 
      link_df$m3_site_spec_uc + 
      link_df$m3_tip_mn +
      link_df$m3_tip2_mn +
      link_df$m3_gen_mn +
      link_df$m3_spec_mn +
      link_df$m3_spec2_mn +
      link_df$m3_reg_mn +
      link_df$m3_coef
    
    
    tmp_trends = data.frame(
      code = a,
      run = b,
      site = uni_sites_sample,
      true = sites[which(sites$ID == uni_sites_sample),]$mean,
      true_lc = sites[which(sites$ID == uni_sites_sample),]$`0.025quant`,
      true_uc = sites[which(sites$ID == uni_sites_sample),]$`0.975quant`,
      m2 = link_df$pop_trend2,
      m2_lc = link_df$pop_trend2_lc,
      m2_uc = link_df$pop_trend2_uc,
      m3 = link_df$pop_trend3,
      m3_lc = link_df$pop_trend3_lc,
      m3_uc = link_df$pop_trend3_uc,
      rarity = link_df$sum
    )
    cmb_trends = rbind(cmb_trends, tmp_trends)
  }
  cmb_df = rbind(cmb_df, cmb_trends)
}
write.csv(cmb_trends, "outputs/model_summary_miss_trends.csv")
