# Control measures (scenario 1):

nbsimu_scenario_1=300
t_end=416

f_high_betw = carac_f[order(carac_f$betweeness,decreasing=T),"lieu"][1:100]
f_high_ind = carac_f[order(carac_f$indegree,decreasing=T),"lieu"][1:100]
f_high_outd = carac_f[order(carac_f$outdegree,decreasing=T),"lieu"][1:100]
f_high_outf = carac_f[order(carac_f$outflux,decreasing=T),"lieu"][1:100]

f_cont_meas = f_high_betw # or f_high_ind, or f_high_outd, or f_high_outf

simulation_scenario_1_control_measures = function(beta_modif){

  for (simu in 1:nbsimu_scenario_1){
    
    f_cont=sample(farms$lieu[which(Reduce("+",apop_initial)["PX",]>0)],round(0.05*nrow(farms))) # Farms initially contaminated
    # f_cont_meas = sample(x=farms$lieu, size=100) # In case farms where control measures are applied are chosen randomly
    
    # 4D matrix with pig populations for each week t (prevalence)
    apop_prev=array(dim=c(2,2,I,as.numeric(t_end-t_start+1)),data=0)
    dimnames(apop_prev)=list(c("Breeding herd","Production herd"),c("PX","PR"),as.character(farms$lieu),t_start:t_end)
    # PX: pigs that don't carry MRSA
    # PR: pigs that carry MRSA
    
    # 3D matrix with prevalence in the batch about to leave for the slaughterhouse
    to_abat_prev=array(dim=c(2,I,as.numeric(t_end-t_start+1)),data=NA)
    dimnames(to_abat_prev)=list(c("PX","PR"),as.character(farms$lieu),t_start:t_end)
    
    # Matrice apop reinitialisee
    apop=apop_initial
    
    # Farm age vector (randomly shifted):
    agep=sample(x=1:4,size=I,replace=T)
    names(agep)=as.character(farms$lieu)
    
    # For farms initially contaminated:
    #####
    pigs=array(0,dim=c(length(f_cont),12))
    dimnames(pigs)=list(f_cont,1:12)
    for(i in 1:12){pigs[f_cont,i]=apop[[i]]["PX",f_cont]}
    for (f in f_cont){
      nporcs=sum(pigs[f,])
      draw=pigs[f,]+1
      while(any((pigs[f,]-draw)<0)){
        draw=t(rmultinom(n=1,size=floor(nporcs*0.16),prob=pigs[f,]))
      }
      for(i in 1:12){
        apop[[i]]["PX",f]=apop[[i]]["PX",f]-draw[i]
        apop[[i]]["PR",f]=apop[[i]]["PR",f]+draw[i]
      }
    }
    for(i in 1:12){if(any(apop[[i]]<0)){print("Erreur dans l'etat initial des elevages contamines")}}
    rm(draw)
    rm(nporcs)
    rm(pigs)
    #####
    
    for (t in t_start:t_end){
      
      if(t<=207){
        beta_in_farms_with_contmeas = betaPPR
      }else{
        beta_in_farms_with_contmeas = beta_modif
      }
      
      # Saves population on week t
      apop_prev["Breeding herd",,,t]=apop$G1+apop$G2+apop$G3+apop$G4+apop$LA
      apop_prev["Production herd",,,t]=apop$FA+apop$PW1+apop$PW2+apop$F1+apop$F2+apop$F3+apop$F4
      
      # If no colonized pig in the whole network, computation is stopped
      if(sum(Reduce("+",apop)["PR",])==0){
        for(i in t:t_end){
          # apop_prev identical the followink weeks
          apop_prev["Breeding herd",,,i]=apop$G1+apop$G2+apop$G3+apop$G4+apop$LA
          apop_prev["Production herd",,,i]=apop$FA+apop$PW1+apop$PW2+apop$F1+apop$F2+apop$F3+apop$F4
          
          to_abat_prev[,names(which(agep%%4==0)),i]=apop$F4[,names(which(agep%%4==0))]
          
          agep=agep+1
        }
        break
      }
      
      # 1) Demographic dynamics:
      
      if(any(agep_init%%4==0)){
        f_chang=names(which(agep_init%%4==0))
        
        apop2=apop # save
        
        # 1.1) Sows demographics:
        
        #####
        apop$G1[,f_chang]=apop2$LA[,f_chang]
        apop$G2[,f_chang]=apop2$G1[,f_chang]
        apop$G3[,f_chang]=apop2$G2[,f_chang]
        apop$G4[,f_chang]=apop2$G3[,f_chang]
        apop$LA[,f_chang]=apop2$G4[,f_chang]
        
        # Gilts imports
        f_repro_chang=intersect(f_chang,as.character(farms_repro$lieu))
        
        f_exp_gilts_t=f_exp_gilts[which(colSums(apop2$F4[,f_exp_gilts])!=0)]
        
        prop_X_F4=apop2$F4["PX",f_exp_gilts_t]/(colSums(apop2$F4[,f_exp_gilts_t]))
        prop_R_F4=apop2$F4["PR",f_exp_gilts_t]/(colSums(apop2$F4[,f_exp_gilts_t]))
        
        p_gilts_X=colSums(mov_gilts[f_exp_gilts_t,f_repro_chang]*prop_X_F4)
        p_gilts_R=colSums(mov_gilts[f_exp_gilts_t,f_repro_chang]*prop_R_F4)
        
        # p_gilts_X+p_gilts_R must be >0
        f_imp_gilts=intersect(f_repro_chang,names(which(p_gilts_X+p_gilts_R!=0)))
        
        # Draw number of gilts X and Y imported
        gilts_draw=rmultinomial(n=length(f_imp_gilts),size=n_imp_gilts[f_imp_gilts],prob=cbind(p_gilts_X[f_imp_gilts],p_gilts_R[f_imp_gilts]))
        gilts_X=gilts_draw[,1]
        gilts_R=gilts_draw[,2]
        
        names(gilts_X)=names(n_imp_gilts[f_imp_gilts])
        names(gilts_R)=names(n_imp_gilts[f_imp_gilts])
        
        # Culled sows:
        culling_draw=rmultinomial(n=length(f_imp_gilts),size=n_imp_gilts[f_imp_gilts],prob=cbind(apop2$LA["PX",f_imp_gilts],apop2$LA["PR",f_imp_gilts]))
        
        apop$G1["PX",f_imp_gilts]=apop2$LA["PX",f_imp_gilts]-culling_draw[,1]+gilts_X
        apop$G1["PR",f_imp_gilts]=apop2$LA["PR",f_imp_gilts]-culling_draw[,2]+gilts_R
        if(any(apop$G1<0)){
          print("Populations G1 <0:")
          print(summary(apop$PW1[which(apop$G1<0)]))
          apop$G1[which(apop$G1<0)]=0
        }
        # /!\ Beware of rounded proportions    
        #####
        
        # 1.2) Fattening pigs demographics:
        
        # Piglets born have same status than mother:
        apop$FA[,f_chang]=round(11.8*apop2$G4[,f_chang])
        
        apop$PW1[,f_chang]=apop2$FA[,f_chang]
        apop$PW2[,f_chang]=apop2$PW1[,f_chang]
        apop$F1[,f_chang]=apop2$PW2[,f_chang]
        apop$F2[,f_chang]=apop2$F1[,f_chang]
        apop$F3[,f_chang]=apop2$F2[,f_chang]
        apop$F4[,f_chang]=apop2$F3[,f_chang]
        
        
        # 1.2.1) Exports AB
        
        #####
        f_chang_exportAB=Reduce(intersect,list(f_chang,f_exportAB,names(which(colSums(apop2$FA)!=0))))
        
        n_exportAB=round(prop_export[f_chang_exportAB]*colSums(apop2$FA[,f_chang_exportAB]))
        
        draw_leave_farm=rmultinomial(n=length(f_chang_exportAB),size=n_exportAB,prob=cbind(apop2$FA["PX",f_chang_exportAB],apop2$FA["PR",f_chang_exportAB]))
        apop$PW1[,f_chang_exportAB]=apop2$FA[,f_chang_exportAB]-t(draw_leave_farm)
        if(any(apop$PW1<0)){
          print("Populations PW1 <0:")
          print(summary(apop$PW1[which(apop$PW1<0)]))
          apop$PW1[which(apop$PW1<0)]=0
        }
        # /!\ Beware of multinomial draws    
        
        movAB_number=floor(movAB[f_chang_exportAB,f_imp_AB]*n_exportAB)
        
        prop_X_AB=draw_leave_farm[,1]/(rowSums(draw_leave_farm))
        prop_X_AB[which(!(is.finite(prop_X_AB)))]=0 # Null origin population ==> NaN = 0
        prop_R_AB=draw_leave_farm[,2]/(rowSums(draw_leave_farm))
        prop_R_AB[which(!(is.finite(prop_R_AB)))]=0
        
        p_AB_X=colSums(movAB_number*prop_X_AB) # Column by column multiplication then sum
        p_AB_R=colSums(movAB_number*prop_R_AB)
        
        f_imp_AB_t=names(which((p_AB_X+p_AB_R)!=0))
        
        AB_draw=rmultinomial(n=length(f_imp_AB_t),size=colSums(movAB_number[,f_imp_AB_t]),prob=cbind(p_AB_X[f_imp_AB_t],p_AB_R[f_imp_AB_t]))
        
        apop$PW1[,f_imp_AB_t]=apop$PW1[,f_imp_AB_t]+t(AB_draw)
        
        #####
        
        # 1.2.2) Exports BC
        
        #####
        f_chang_exportBC=Reduce(intersect,list(f_chang,f_exportBC,names(which(colSums(apop2$PW2)!=0))))
        
        n_exportBC=floor(prop_export[f_chang_exportBC]*colSums(apop2$PW2[,f_chang_exportBC]))
        
        draw_leave_farm=rmultinomial(n=length(f_chang_exportBC),size=n_exportBC,prob=cbind(apop2$PW2["PX",f_chang_exportBC],apop2$PW2["PR",f_chang_exportBC]))
        apop$F1[,f_chang_exportBC]=apop2$PW2[,f_chang_exportBC]-t(draw_leave_farm)
        if(any(apop$F1<0)){
          print("Populations F1 <0:")
          print(summary(apop$F1[which(apop$F1<0)]))
          apop$F1[which(apop$F1<0)]=0
        }
        # /!\ Beware of multinomial draws    
        
        movBC_number=floor(movBC[f_chang_exportBC,f_imp_BC]*n_exportBC)
        
        prop_X_BC=draw_leave_farm[,1]/(rowSums(draw_leave_farm))
        prop_X_BC[which(!(is.finite(prop_X_BC)))]=0
        prop_R_BC=draw_leave_farm[,2]/(rowSums(draw_leave_farm))
        prop_R_BC[which(!(is.finite(prop_R_BC)))]=0
        
        p_BC_X=colSums(movBC_number*prop_X_BC)
        p_BC_R=colSums(movBC_number*prop_R_BC)
        
        f_imp_BC_t=names(which((p_BC_X+p_BC_R)!=0))
        
        BC_draw=rmultinomial(n=length(f_imp_BC_t),size=colSums(movBC_number[,f_imp_BC_t]),prob=cbind(p_BC_X[f_imp_BC_t],p_BC_R[f_imp_BC_t]))
        
        apop$F1[,f_imp_BC_t]=apop$F1[,f_imp_BC_t]+t(BC_draw)
        
        #####
        
        for(i in 1:12){if(any(apop[[i]]<0)){print(paste("Population negative dans le compartiment:",names(apop)[i]))}}
        
      }
      
      # 2) Epidemic dynamics:
      
      #####
      
      apop2=apop # save
      
      # Post-Weaning compartments:
      
      epi_output=epidemio_cont_meas(compart=apop2$PW1,contact_compart=list(apop2$PW1,apop2$PW2),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$PW1=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$PW2,contact_compart=list(apop2$PW1,apop2$PW2),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$PW2=epi_output[[1]]
      
      # Finishing compartments:
      
      epi_output=epidemio_cont_meas(compart=apop2$F1,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$F1=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$F2,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$F2=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$F3,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$F3=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$F4,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$F4=epi_output[[1]]
      
      # Gestation compartments:
      
      epi_output=epidemio_cont_meas(compart=apop2$G1,contact_compart=list(apop2$G1,apop2$G2,apop2$G3,apop2$G4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$G1[,farms_repro$lieu]=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$G2,contact_compart=list(apop2$G1,apop2$G2,apop2$G3,apop2$G4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$G2[,farms_repro$lieu]=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$G3,contact_compart=list(apop2$G1,apop2$G2,apop2$G3,apop2$G4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$G3[,farms_repro$lieu]=epi_output[[1]]
      
      epi_output=epidemio_cont_meas(compart=apop2$G4,contact_compart=list(apop2$G1,apop2$G2,apop2$G3,apop2$G4),beta=betaPPR,beta_in_farms_with_contmeas)
      apop$G4[,farms_repro$lieu]=epi_output[[1]]
      
      #####
      
      t=t+1
      agep=agep+1
    }
    
    # Saves total pig population:
    pop_tot=rbind(colSums(apop_prev["Breeding herd",,,],dims=2),colSums(apop_prev["Production herd",,,],dims=2))
    dimnames(pop_tot)=list(c("Breeding herd","Production herd"),t_start:t_end)
    
    # Saves all farms on first repetition of the model, and then only farms with at least one PR pig at some point, and the 2 first in the list
    if(simu>=2){
      f_pigsR_nonemp=union(farms$lieu[1:2],farms$lieu[which(rowSums(colSums(apop_prev[,"PR",,]))!=0)])
      apop_prev=apop_prev[,,f_pigsR_nonemp,]
    }
    
    # Output matrices:
    results=list(apop_prev,1,pop_tot,to_abat_prev)
    # After each repetition, save:
    save(results,file=paste("~/",beta_modif,"/Results_",beta_modif,"_simu_",simu,".rdata",sep=""))
  }
}

######################################################################

# Analysis and figure 6:

library(mratios)

# Function:
# "output" = "n_f_cont" or "apop_prev"
eff_beta = function(vbeta, output, v_t_val){
  print(paste0("Analyse pour beta=",vbeta,", output: ", output," et valeurs de temps:"))
  print(v_t_val)
  
  sim_outp_rand = array(data=0, dim=c(length(v_t_val), nbsimu_scenario_1))
  sim_outp_outdeg = array(data=0, dim=c(length(v_t_val), nbsimu_scenario_1))
  sim_outp_betw = array(data=0, dim=c(length(v_t_val), nbsimu_scenario_1))
  sim_outp_outflux = array(data=0, dim=c(length(v_t_val), nbsimu_scenario_1))
  sim_outp_indegree = array(data=0, dim=c(length(v_t_val), nbsimu_scenario_1))
  sim_outp_baseline = array(data=0, dim=c(length(v_t_val), nbsimu_scenario_1))
  
  rownames(sim_outp_rand) = as.character(v_t_val)
  rownames(sim_outp_outdeg) = as.character(v_t_val)
  rownames(sim_outp_betw) = as.character(v_t_val)
  rownames(sim_outp_outflux) = as.character(v_t_val)
  rownames(sim_outp_indegree) = as.character(v_t_val)
  rownames(sim_outp_baseline) = as.character(v_t_val)
  
  for (simu in 1:nbsimu_scenario_1){

    # Baseline:
    load(paste0("~/Random farms/0.262/Results_0.262_simu_",simu,".rdata"))
    for (t_val in v_t_val){
      sim_outp_baseline[as.character(t_val),simu] = val_output(output,results,t_val)
    }
    
    # Scenarios "random" and based on criteria:
    load(paste0("~/Random farms/",as.character(vbeta),"/Results_",as.character(vbeta),"_simu_",simu,".rdata"))
    for (t_val in v_t_val){
      sim_outp_rand[as.character(t_val),simu] = val_output(output,results,t_val)
    }
    
    load(paste0("~/Highest outdegree farms/",as.character(vbeta),"/Results_",as.character(vbeta),"_simu_",simu,".rdata"))
    for (t_val in v_t_val){
      sim_outp_outdeg[as.character(t_val),simu] = val_output(output,results,t_val)
    }
    
    load(paste0("~/Highest betweeness farms/",as.character(vbeta),"/Results_",as.character(vbeta),"_simu_",simu,".rdata"))
    for (t_val in v_t_val){
      sim_outp_betw[as.character(t_val),simu] = val_output(output,results,t_val)
    }
    
    load(paste0("~/Highest outflux farms/",as.character(vbeta),"/Results_",as.character(vbeta),"_simu_",simu,".rdata"))
    for (t_val in v_t_val){
      sim_outp_outflux[as.character(t_val),simu] = val_output(output,results,t_val)
    }
    
    load(paste0("~/Highest indegree farms/",as.character(vbeta),"/Results_",as.character(vbeta),"_simu_",simu,".rdata"))
    for (t_val in v_t_val){
      sim_outp_indegree[as.character(t_val),simu] = val_output(output,results,t_val)
    }
  }
  
  return(list(sim_outp_baseline, sim_outp_rand, sim_outp_outdeg, sim_outp_betw, sim_outp_outflux, sim_outp_indegree))
}

which_output = "n_f_cont" # or "apop_prev"
all_indicators = T  # or "F"

if(which_output == "apop_prev"){
  val_beta_50 = eff_beta(vbeta = 0.131, output = "apop_prev", v_t_val = c(208,260,312,364,416))
  val_beta_25 = eff_beta(vbeta = 0.1965, output = "apop_prev", v_t_val = c(208,260,312,364,416))
}else if(which_output == "n_f_cont"){
  val_beta_50 = eff_beta(vbeta = 0.131, output = "n_f_cont", v_t_val = c(208,260,312,364,416))
  val_beta_25 = eff_beta(vbeta = 0.1965, output = "n_f_cont", v_t_val = c(208,260,312,364,416))
}

val_beta_baseline = val_beta_50[[1]]

val_beta_50_rand = val_beta_50[[2]]
val_beta_50_outdeg = val_beta_50[[3]]
val_beta_50_betw = val_beta_50[[4]]
val_beta_50_outflux = val_beta_50[[5]]
val_beta_50_indegree = val_beta_50[[6]]

val_beta_25_rand = val_beta_25[[2]]
val_beta_25_outdeg = val_beta_25[[3]]
val_beta_25_betw = val_beta_25[[4]]
val_beta_25_outflux = val_beta_25[[5]]
val_beta_25_indegree = val_beta_25[[6]]

cont_meas = as.data.frame(matrix(NA,30,7))
colnames(cont_meas) = c("reduction", "time", "targeted", "ratio", "inf", "sup", "pvalue")

# With "1 year", "2 years" and "4 years"
#####
cont_meas$reduction = c(rep("Beta reduction of 50%",15), rep("Beta reduction of 25%",15))
cont_meas$time = c(rep("1 year",5), rep("2 years",5), rep("4 years",5), rep("1 year",5), rep("2 years",5), rep("4 years",5))
cont_meas$targeted = rep(c("Random","Highest outdegree","Highest betweeness","Highest outflux", "Highest indegree"), 6)

val_test = ttestratio(x = val_beta_50_rand["260",], y = val_beta_baseline["260",])
cont_meas$ratio[1] = val_test$estimate["x/y"]
cont_meas[1, c(5,6)] = val_test$conf.int
cont_meas$pvalue[1] = val_test$p.value

val_test = ttestratio(val_beta_50_outdeg["260",], val_beta_baseline["260",])
cont_meas$ratio[2] = val_test$estimate["x/y"]
cont_meas[2, c(5,6)] = val_test$conf.int
cont_meas$pvalue[2] = val_test$p.value

val_test = ttestratio(val_beta_50_betw["260",], val_beta_baseline["260",])
cont_meas$ratio[3] = val_test$estimate["x/y"]
cont_meas[3, c(5,6)] = val_test$conf.int
cont_meas$pvalue[3] = val_test$p.value

val_test = ttestratio(val_beta_50_outflux["260",], val_beta_baseline["260",])
cont_meas$ratio[4] = val_test$estimate["x/y"]
cont_meas[4, c(5,6)] = val_test$conf.int
cont_meas$pvalue[4] = val_test$p.value

val_test = ttestratio(val_beta_50_indegree["260",], val_beta_baseline["260",])
cont_meas$ratio[5] = val_test$estimate["x/y"]
cont_meas[5, c(5,6)] = val_test$conf.int
cont_meas$pvalue[5] = val_test$p.value

val_test = ttestratio(val_beta_50_rand["312",], val_beta_baseline["312",])
cont_meas$ratio[6] = val_test$estimate["x/y"]
cont_meas[6, c(5,6)] = val_test$conf.int
cont_meas$pvalue[6] = val_test$p.value

val_test = ttestratio(val_beta_50_outdeg["312",], val_beta_baseline["312",])
cont_meas$ratio[7] = val_test$estimate["x/y"]
cont_meas[7, c(5,6)] = val_test$conf.int
cont_meas$pvalue[7] = val_test$p.value

val_test = ttestratio(val_beta_50_betw["312",], val_beta_baseline["312",])
cont_meas$ratio[8] = val_test$estimate["x/y"]
cont_meas[8, c(5,6)] = val_test$conf.int
cont_meas$pvalue[8] = val_test$p.value

val_test = ttestratio(val_beta_50_outflux["312",], val_beta_baseline["312",])
cont_meas$ratio[9] = val_test$estimate["x/y"]
cont_meas[9, c(5,6)] = val_test$conf.int
cont_meas$pvalue[9] = val_test$p.value

val_test = ttestratio(val_beta_50_indegree["312",], val_beta_baseline["312",])
cont_meas$ratio[10] = val_test$estimate["x/y"]
cont_meas[10, c(5,6)] = val_test$conf.int
cont_meas$pvalue[10] = val_test$p.value

val_test = ttestratio(val_beta_50_rand["416",], val_beta_baseline["416",])
cont_meas$ratio[11] = val_test$estimate["x/y"]
cont_meas[11, c(5,6)] = val_test$conf.int
cont_meas$pvalue[11] = val_test$p.value

val_test = ttestratio(val_beta_50_outdeg["416",], val_beta_baseline["416",])
cont_meas$ratio[12] = val_test$estimate["x/y"]
cont_meas[12, c(5,6)] = val_test$conf.int
cont_meas$pvalue[12] = val_test$p.value

val_test = ttestratio(val_beta_50_betw["416",], val_beta_baseline["416",])
cont_meas$ratio[13] = val_test$estimate["x/y"]
cont_meas[13, c(5,6)] = val_test$conf.int
cont_meas$pvalue[13] = val_test$p.value

val_test = ttestratio(val_beta_50_outflux["416",], val_beta_baseline["416",])
cont_meas$ratio[14] = val_test$estimate["x/y"]
cont_meas[14, c(5,6)] = val_test$conf.int
cont_meas$pvalue[14] = val_test$p.value

val_test = ttestratio(val_beta_50_indegree["416",], val_beta_baseline["416",])
cont_meas$ratio[15] = val_test$estimate["x/y"]
cont_meas[15, c(5,6)] = val_test$conf.int
cont_meas$pvalue[15] = val_test$p.value

val_test = ttestratio(val_beta_25_rand["260",], val_beta_baseline["260",])
cont_meas$ratio[16] = val_test$estimate["x/y"]
cont_meas[16, c(5,6)] = val_test$conf.int
cont_meas$pvalue[16] = val_test$p.value

val_test = ttestratio(val_beta_25_outdeg["260",], val_beta_baseline["260",])
cont_meas$ratio[17] = val_test$estimate["x/y"]
cont_meas[17, c(5,6)] = val_test$conf.int
cont_meas$pvalue[17] = val_test$p.value

val_test = ttestratio(val_beta_25_betw["260",], val_beta_baseline["260",])
cont_meas$ratio[18] = val_test$estimate["x/y"]
cont_meas[18, c(5,6)] = val_test$conf.int
cont_meas$pvalue[18] = val_test$p.value

val_test = ttestratio(val_beta_25_outflux["260",], val_beta_baseline["260",])
cont_meas$ratio[19] = val_test$estimate["x/y"]
cont_meas[19, c(5,6)] = val_test$conf.int
cont_meas$pvalue[19] = val_test$p.value

val_test = ttestratio(val_beta_25_indegree["260",], val_beta_baseline["260",])
cont_meas$ratio[20] = val_test$estimate["x/y"]
cont_meas[20, c(5,6)] = val_test$conf.int
cont_meas$pvalue[20] = val_test$p.value

val_test = ttestratio(val_beta_25_rand["312",], val_beta_baseline["312",])
cont_meas$ratio[21] = val_test$estimate["x/y"]
cont_meas[21, c(5,6)] = val_test$conf.int
cont_meas$pvalue[21] = val_test$p.value

val_test = ttestratio(val_beta_25_outdeg["312",], val_beta_baseline["312",])
cont_meas$ratio[22] = val_test$estimate["x/y"]
cont_meas[22, c(5,6)] = val_test$conf.int
cont_meas$pvalue[22] = val_test$p.value

val_test = ttestratio(val_beta_25_betw["312",], val_beta_baseline["312",])
cont_meas$ratio[23] = val_test$estimate["x/y"]
cont_meas[23, c(5,6)] = val_test$conf.int
cont_meas$pvalue[23] = val_test$p.value

val_test = ttestratio(val_beta_25_outflux["312",], val_beta_baseline["312",])
cont_meas$ratio[24] = val_test$estimate["x/y"]
cont_meas[24, c(5,6)] = val_test$conf.int
cont_meas$pvalue[24] = val_test$p.value

val_test = ttestratio(val_beta_25_indegree["312",], val_beta_baseline["312",])
cont_meas$ratio[25] = val_test$estimate["x/y"]
cont_meas[25, c(5,6)] = val_test$conf.int
cont_meas$pvalue[25] = val_test$p.value

val_test = ttestratio(val_beta_25_rand["416",], val_beta_baseline["416",])
cont_meas$ratio[26] = val_test$estimate["x/y"]
cont_meas[26, c(5,6)] = val_test$conf.int
cont_meas$pvalue[26] = val_test$p.value

val_test = ttestratio(val_beta_25_outdeg["416",], val_beta_baseline["416",])
cont_meas$ratio[27] = val_test$estimate["x/y"]
cont_meas[27, c(5,6)] = val_test$conf.int
cont_meas$pvalue[27] = val_test$p.value

val_test = ttestratio(val_beta_25_betw["416",], val_beta_baseline["416",])
cont_meas$ratio[28] = val_test$estimate["x/y"]
cont_meas[28, c(5,6)] = val_test$conf.int
cont_meas$pvalue[28] = val_test$p.value

val_test = ttestratio(val_beta_25_outflux["416",], val_beta_baseline["416",])
cont_meas$ratio[29] = val_test$estimate["x/y"]
cont_meas[29, c(5,6)] = val_test$conf.int
cont_meas$pvalue[29] = val_test$p.value

val_test = ttestratio(val_beta_25_indegree["416",], val_beta_baseline["416",])
cont_meas$ratio[30] = val_test$estimate["x/y"]
cont_meas[30, c(5,6)] = val_test$conf.int
cont_meas$pvalue[30] = val_test$p.value
#####

cont_meas$targeted = factor(cont_meas$targeted, levels = c("Random", "Highest outdegree", "Highest indegree", "Highest betweeness", "Highest outflux"))

library(ggplot2)

if(all_indicators){
  p = ggplot(data = cont_meas, aes(x = time, colour = targeted))
}else{
  p = ggplot(data = cont_meas[which(cont_meas$targeted%in%c("Random","Highest outdegree","Highest indegree")),], aes(x = time, colour = targeted))
}

if(which_output == "apop_prev"){
  p = p + ggtitle("Effect of control measures on the total prevalence")
  p = p + ylab("Prevalence ratio (with measures / without measures)")
}else if(which_output == "n_f_cont"){
  p = p + ggtitle("Effect of control measures on the PFC")
  p = p + ylab("PFC ratio (with measures / without measures)")
}

p = p + xlab("Time after control measures introduction")
p = p + geom_errorbar(aes(ymin = inf, ymax = sup), width = 0.3, size = 0.8, position = position_dodge(0.9))
p = p + geom_point(aes(y = ratio, shape = targeted), stat="identity", size = 2, position = position_dodge(0.9))
p = p + labs(colour = "Selection of the farms:", shape = "Selection of the farms:")
p = p + scale_shape_manual(values = c(15, 16, 17, 4, 18))
p = p + facet_wrap(~reduction)
p = p + geom_hline(yintercept = 1, linetype = "dashed", size = 1)
p = p + geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed", size = 0.5, col = "grey60")
p = p + theme(legend.position="bottom")
p

