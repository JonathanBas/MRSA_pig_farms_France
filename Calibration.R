# Calibration (scenario 1):

nbsimu_scenario_1=300
t_end=416

simulation_scenario_1_calibration = function(betap){
  
  betaPPR=betap
  
  for (simu in 1:nbsimu_scenario_1){
    
    f_cont=sample(farms$lieu[which(Reduce("+",apop_initial)["PX",]>0)],round(0.05*nrow(farms))) # Farms initially contaminated
    
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
      
      epi_output=epidemio(compart=apop2$PW1,contact_compart=list(apop2$PW1,apop2$PW2),beta=betaPPR)
      apop$PW1=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$PW2,contact_compart=list(apop2$PW1,apop2$PW2),beta=betaPPR)
      apop$PW2=epi_output[[1]]
      
      # Finishing compartments:
      
      epi_output=epidemio(compart=apop2$F1,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR)
      apop$F1=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$F2,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR)
      apop$F2=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$F3,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR)
      apop$F3=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$F4,contact_compart=list(apop2$F1,apop2$F2,apop2$F3,apop2$F4),beta=betaPPR)
      apop$F4=epi_output[[1]]
      
      # Gestation compartments:
      
      epi_output=epidemio(compart=apop2$G1[,farms_repro$lieu],contact_compart=list(apop2$G1[,farms_repro$lieu],apop2$G2[,farms_repro$lieu],apop2$G3[,farms_repro$lieu],apop2$G4[,farms_repro$lieu]),beta=betaPPR)
      apop$G1[,farms_repro$lieu]=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$G2[,farms_repro$lieu],contact_compart=list(apop2$G1[,farms_repro$lieu],apop2$G2[,farms_repro$lieu],apop2$G3[,farms_repro$lieu],apop2$G4[,farms_repro$lieu]),beta=betaPPR)
      apop$G2[,farms_repro$lieu]=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$G3[,farms_repro$lieu],contact_compart=list(apop2$G1[,farms_repro$lieu],apop2$G2[,farms_repro$lieu],apop2$G3[,farms_repro$lieu],apop2$G4[,farms_repro$lieu]),beta=betaPPR)
      apop$G3[,farms_repro$lieu]=epi_output[[1]]
      
      epi_output=epidemio(compart=apop2$G4[,farms_repro$lieu],contact_compart=list(apop2$G1[,farms_repro$lieu],apop2$G2[,farms_repro$lieu],apop2$G3[,farms_repro$lieu],apop2$G4[,farms_repro$lieu]),beta=betaPPR)
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
    save(results,file=paste("~/",betaPPR,"/Results_",betaPPR,"_simu_",simu,".rdata",sep=""))
  }
}

curve_calib_beta = function(betap){

  calib_beta = array(data=0, dim=c(t_end,nbsimu_scenario_1))

  for (simu in 1:nbsimu_scenario_1){
    load(paste0("~/",betap,"/Results_",betap,"_simu_",simu,".rdata",sep=""))
    calib_beta[,simu] = val_output_2("apop_prev",results)
  }
  
  beta_prev = as.data.frame(matrix(0,t_end,4))
  colnames(beta_prev) = c("week","mean","inf","sup")
  
  beta_prev$week = 1:t_end
  
  beta_prev$mean = rowMeans(calib_beta)
  for (t in 1:t_end){
    beta_prev$inf[t] = quantile(calib_beta[t,], probs=0.025)
    beta_prev$sup[t] = quantile(calib_beta[t,], probs=0.975)
  }

  ecart = sum((beta_prev$mean[(208-52):208]-0.8)^2)
  print(paste("Residuals square sum:", ecart))
  
  return(beta_prev)
}
