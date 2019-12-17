# Sensitivity analysis:
# Model's behavior under changes in the value of Beta (see SM4):

value_of_beta = 0.1 # or other values
nbsimu_scenario_1 = 300
t_end = 416

free_simus = matrix(NA, nsimu, t_end)
cont_simus = matrix(NA, nsimu, t_end)
decont_simus = matrix(NA, nsimu, t_end)

for (simu in 1:nsimu){
  print(paste0(simu,"/",nsimu))
  load(paste0("~/",value_of_beta,"/Results_",value_of_beta,"_simu_",simu,".rdata"))
  apop_prev = results[[1]]
  
  was_farm_cont = rep(F, dim(apop_prev)[3])
  for (t in 1:t_end){
    is_farm_cont = (colSums(apop_prev[,"PR",,t]) !=0 )
    cont_simus[simu, t] = sum(is_farm_cont)
    
    decont_simus[simu, t] = length(which((is_farm_cont==F) & (was_farm_cont==T)))
    
    free_simus[simu, t] = 10542 - cont_simus[simu, t] - decont_simus[simu, t]
    
    was_farm_cont[is_farm_cont] = T
  }
}

cont_ext = as.data.frame(matrix(NA, t_end, 10))
colnames(cont_ext) = c("week","free","free_inf","free_sup","cont","cont_inf","cont_sup","decont","decont_inf","decont_sup")
cont_ext$week = 1:t_end
for (t in 1:t_end){
  cont_ext$free[t] = mean(free_simus[,t])
  cont_ext$free_inf[t] = quantile(free_simus[,t], probs = 0.025)
  cont_ext$free_sup[t] = quantile(free_simus[,t], probs = 0.975)
  
  cont_ext$cont[t] = mean(cont_simus[,t])
  cont_ext$cont_inf[t] = quantile(cont_simus[,t], probs = 0.025)
  cont_ext$cont_sup[t] = quantile(cont_simus[,t], probs = 0.975)
  
  cont_ext$decont[t] = mean(decont_simus[,t])
  cont_ext$decont_inf[t] = quantile(decont_simus[,t], probs = 0.025)
  cont_ext$decont_sup[t] = quantile(decont_simus[,t], probs = 0.975)
}

cont_ext[,-1] = cont_ext[,-1]*100/10542 # %

# Plot

library(ggplot2)
p = ggplot(data = cont_ext, aes(x=week))
p = p + xlab("Time (weeks)") + ylab("Percentage of farms (%)") + ylim(0,100) + ggtitle(paste0("Beta = ",value_of_beta))

p = p+ geom_ribbon(aes(ymin=free_inf, ymax=free_sup), color="black", alpha=1/4, fill="#E69F00")
p = p + geom_line(aes(y=free, col="1"), size=1)

p = p+ geom_ribbon(aes(ymin=cont_inf, ymax=cont_sup), color="black", alpha=1/4, fill="#56B4E9")
p = p + geom_line(aes(y=cont, col="2"), size=1)

p = p+ geom_ribbon(aes(ymin=decont_inf, ymax=decont_sup), color="black", alpha=1/4, fill="#CC79A7")
p = p + geom_line(aes(y=decont, col="3"), size=1)

p = p + scale_colour_manual(values = c("1"="#E69F00", "2"="#56B4E9", "3"="#CC79A7"),
                            name = "Farms:",
                            labels = c("Never contaminated", "Contaminated", "Decontaminated"))
p

######################################################################
# For the sensitivity analysis assessing the impact of our initial hypothesis on the value of the carriage duration, previously defined functions (including "Calibration" and "Control measures" have to run with a different value for VPR)
