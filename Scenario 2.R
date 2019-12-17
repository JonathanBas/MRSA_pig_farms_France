# Scenario 2:

nbsimu_scenario_2 = 50

simulation_scenario_2 = function(X){
  
  f_cont=farms_repro$lieu[X] # Xth seed farm in the list of farms with sows
  
  dir.create(paste0("~/",f_cont))
  
  Categ=farms$site[which(farms$lieu==f_cont)]
  
  t_end=as.numeric(t_equilibre[f_cont])

  for (simu in 1:nbsimu_scenario_2){
    
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
    
    agep[f_cont]=4
    
    apop$G1["PX",f_cont]=apop$G1["PX",f_cont]-n_imp_gilts[f_cont]
    apop$G1["PR",f_cont]=n_imp_gilts[f_cont]
    
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
    results=list(apop_prev,0,pop_tot,to_abat_prev)
    # After each repetition, save:
    save(results,file=paste("~/",f_cont,"/Results_",Categ,"_f_",f_cont,"_simu_",simu,".rdata",sep=""))
  }
}

######################################################################

mat_analyse=as.data.frame(matrix(data=0,nrow=length(farms_repro$lieu),ncol=16,dimnames=list(farms_repro$lieu,c("lieu","site","size_bre","size_pro","outdegree","indegree","betweeness","closeness","coreness","eigcen","outflux","influx","numb_cont_farms","prev_slaugh","prev_all_ages","prev_farmer"))))
carac_f=as.data.frame(matrix(data=0,nrow=length(farms$lieu),ncol=19,dimnames=list(farms$lieu,c("lieu","site","size_bre","size_pro","outdegree","indegree","betweeness","closeness","coreness","eigcen","outflux","influx","times_cont","t_of_cont","l_of_cont","importance_of_cont","nb_cont_f_when_cont","Bree_first","Prod_first"))))

# For getting the proportion of each farm category in the contaminated farms:
mat_eff_seed_carac = as.data.frame(matrix(data=0,nrow=length(farms_repro$lieu),ncol=13,dimnames=list(farms_repro$lieu,c("lieu","site","outdegree","numb_cont_farms","SEL","MU","FF","FA","FPW","PW","PWF","FI","prev_slaugh"))))

# Computation of "mat_analyse", and then "carac_f":
for(f in farms_repro$lieu){
  iter=which(farms_repro$lieu==f)
  Categ=farms_repro$site[which(farms_repro$lieu==f)]
  
  mat_analyse$lieu[iter]=f
  mat_analyse$site[iter]=Categ
  mat_analyse$size_bre[iter]=apop_initial$G1["PX",f]+apop_initial$G2["PX",f]+apop_initial$G3["PX",f]+apop_initial$G4["PX",f]+apop_initial$LA["PX",f]
  mat_analyse$size_pro[iter]=apop_initial$FA["PX",f]+apop_initial$PW1["PX",f]+apop_initial$PW2["PX",f]+apop_initial$F1["PX",f]+apop_initial$F2["PX",f]+apop_initial$F3["PX",f]+apop_initial$F4["PX",f]

  mat_eff_seed_carac$lieu[iter]=f
  mat_eff_seed_carac$site[iter]=Categ
  
  for(simu in 1:nbsimu_scenario_2){
    load(paste0("~/",f,"/Results_",Categ,"_f_",f,"_simu_",simu,".rdata"))
    
    mat_analyse$numb_cont_farms[iter]=mat_analyse$numb_cont_farms[iter]+val_output(output="n_f_cont",results,iter,simu)/nbsimu_scenario_2
    mat_analyse$prev_slaugh[iter]=mat_analyse$prev_slaugh[iter]+val_output(output="to_abat_prev",results,iter,simu)/nbsimu_scenario_2
    mat_analyse$prev_all_ages[iter]=mat_analyse$prev_all_ages[iter]+val_output(output="apop_prev",results,iter,simu)/nbsimu_scenario_2

    m_output=results[[1]]
    farms_nonnuls=names(which(rowSums(colSums(m_output[,"PR",,]))!=0))
    farms_nonnuls=farms_nonnuls[which(farms_nonnuls!=f)]

    # mat_eff_seed_carac: at t = 52 weeks
    num_far = which(farms$lieu %in% names(which(colSums(m_output[,"PR",,52])!=0)))
    for (cate in c("SEL","MU","FF","FA","FPW","PW","PWF","FI")){
      mat_eff_seed_carac[iter,cate] = mat_eff_seed_carac[iter,cate] + sum(farms$site[num_far] == cate)/nbsimu_scenario_2
    }

    for(far in farms_nonnuls){
      nfarm=which(farms$lieu==far)

      carac_f$times_cont[nfarm] = carac_f$times_cont[nfarm] + 1 # =1 if farm was contaminated

      weeks_cont=which(colSums(m_output[,"PR",far,])!=0) # Weeks when the farm was contaminated

      if(m_output["Breeding herd","PR",far,min(weeks_cont)] != 0){
        carac_f$Bree_first[nfarm] = carac_f$Bree_first[nfarm] + 1
      }

      if(m_output["Production herd","PR",far,min(weeks_cont)] != 0){
        carac_f$Prod_first[nfarm] = carac_f$Prod_first[nfarm] + 1
      }

      carac_f$t_of_cont[nfarm] = carac_f$t_of_cont[nfarm] + min(weeks_cont) # Semaine ou la ferme a ete contaminee
      # carac_f$l_of_cont[nfarm] = carac_f$l_of_cont[nfarm] + (length(weeks_cont)) #/ dim(m_output)[4]) # Number of weeks the farm was contaminated (out of the total duration of the simulation)

      # Importance of contamination: Mean prevalence of MRSA during each week the farm is contaminated
      if(length(weeks_cont)>1){
        carac_f$importance_of_cont[nfarm] = carac_f$importance_of_cont[nfarm] + mean(colSums(m_output[, "PR", far, weeks_cont])/colSums(m_output[, , far, weeks_cont],dims=2))
      }else{
        carac_f$importance_of_cont[nfarm] = carac_f$importance_of_cont[nfarm] + ( sum(m_output[, "PR", far, weeks_cont]) / sum(m_output[, , far, weeks_cont]) )
      }

      etat_fermes_avant_cont_far=colSums(m_output[,"PR",,1:min(weeks_cont)]) # Number of PR in all farms between week 1 and the contamination week
      carac_f$nb_cont_f_when_cont[nfarm] = carac_f$nb_cont_f_when_cont[nfarm] + length(which(rowSums(etat_fermes_avant_cont_far)!=0)) # Number of farms contaminated when "far" becomes contamined

    }
  }
  
}

for(f in farms$lieu){
  iter=which(farms$lieu==f)
  Categ=farms$site[which(farms$lieu==f)]
  
  carac_f$lieu[iter]=f
  carac_f$site[iter]=Categ
  carac_f$size_bre[iter]=apop_initial$G1["PX",f]+apop_initial$G2["PX",f]+apop_initial$G3["PX",f]+apop_initial$G4["PX",f]+apop_initial$LA["PX",f]
  carac_f$size_pro[iter]=apop_initial$FA["PX",f]+apop_initial$PW1["PX",f]+apop_initial$PW2["PX",f]+apop_initial$F1["PX",f]+apop_initial$F2["PX",f]+apop_initial$F3["PX",f]+apop_initial$F4["PX",f]
}

# Using package "igraph" for centrality indicators:
#####
movtot = movAB+movBC+mov_gilts

tab_graph = data.frame(matrix(0, length(which(movtot!=0)), 3))
colnames(tab_graph)=c('Source','Target','Weight')

iter = 0
for(i in 1:nrow(farms)){
  rec_non_nuls = which(movtot[farms$lieu[i],]!=0)
  
  if(length(rec_non_nuls)!=0){
    tab_graph$Source[ (iter+1) : (iter+length(rec_non_nuls)) ] = rep(farms$lieu[i],length(rec_non_nuls))
    tab_graph$Target[ (iter+1) : (iter+length(rec_non_nuls)) ] = names(rec_non_nuls)
  }
  iter = iter+length(rec_non_nuls)
}

library(igraph)
graph_network = graph_from_edgelist(as.matrix(tab_graph[,1:2]), directed=TRUE)

# Vectors with centrality indicators: betweenness, eigen centrality, closeness, coreness
betw = betweenness(graph_network)
eigen_centr = eigen_centrality(graph=graph_network, directed=T)
closen = closeness(graph = graph_network, mode = "all") # We consider non-directional network (mode=all). Also, non consensual definition of closeness for not connected graphes
coren = coreness(graph = graph_network, mode = "all") # We consider non-directional network (mode=all)

# Idem for degrees
# CAUTION: Degrees do NOT include loops (movements in same farm)
outd = degree(graph_network, mode="out", loops=F)
ind = degree(graph_network, mode="in", loops=F)

for (i in 1:nrow(mat_analyse)){
  if(any(names(betw) == mat_analyse$lieu[i])){
    
    mat_analyse$betweeness[i] = betw[ which(names(betw)==mat_analyse$lieu[i]) ]
    mat_analyse$eigcen[i] = eigen_centr$vector[ which(names(eigen_centr$vector)==mat_analyse$lieu[i]) ]
    mat_analyse$closeness[i] = closen[ which(names(closen)==mat_analyse$lieu[i]) ]
    mat_analyse$coreness[i] = coren[ which(names(coren)==mat_analyse$lieu[i]) ]
    mat_analyse$outdegree[i] = outd[ which(names(outd)==mat_analyse$lieu[i]) ]
    mat_analyse$indegree[i] = ind[ which(names(ind)==mat_analyse$lieu[i]) ]
  }
}

for (i in 1:nrow(carac_f)){
  if(any(names(betw) == carac_f$lieu[i])){
    
    carac_f$betweeness[i] = betw[ which(names(betw)==carac_f$lieu[i]) ]
    carac_f$eigcen[i] = eigen_centr$vector[ which(names(eigen_centr$vector)==carac_f$lieu[i]) ]
    carac_f$closeness[i] = closen[ which(names(closen)==carac_f$lieu[i]) ]
    carac_f$coreness[i] = coren[ which(names(coren)==carac_f$lieu[i]) ]
    carac_f$outdegree[i] = outd[ which(names(outd)==carac_f$lieu[i]) ]
    carac_f$indegree[i] = ind[ which(names(ind)==carac_f$lieu[i]) ]
  }
}
#####

# Similarly, we use package "igraph" for calculating the influx/outflux (= Weighted indegree/outdegree, number of pigs imported/exported by each farm)
# CAUTION: loops NOT included
# /!\ We use non proportional matrices because we need weighted edges
#####
movtot = movAB_non_prop + movBC_non_prop + mov_gilts_non_prop
tab_graph = data.frame(matrix(0, length(which(movtot!=0)), 3))
colnames(tab_graph)=c('Source','Target','Weight')
iter = 0
for(i in 1:nrow(farms)){
  rec_non_nuls = which(movtot[farms$lieu[i],]!=0)
  if(length(rec_non_nuls)!=0){
    tab_graph$Source[ (iter+1) : (iter+length(rec_non_nuls)) ] = rep(farms$lieu[i],length(rec_non_nuls))
    tab_graph$Target[ (iter+1) : (iter+length(rec_non_nuls)) ] = names(rec_non_nuls)
    tab_graph$Weight[ (iter+1) : (iter+length(rec_non_nuls)) ] = movtot[farms$lieu[i],names(rec_non_nuls)]
  }
  iter = iter+length(rec_non_nuls)
}
library(igraph)
graph_network = graph_from_edgelist(as.matrix(tab_graph[,1:2]), directed=TRUE)
E(graph_network)$weight = as.numeric(tab_graph[,3])
infl = strength(graph = graph_network, mode = "in", loops = F)
outfl = strength(graph = graph_network, mode = "out", loops = F)

for (i in 1:nrow(mat_analyse)){
  if(any(names(infl)==mat_analyse$lieu[i])){
    mat_analyse$influx[i] = infl[ which(names(infl)==mat_analyse$lieu[i]) ]
    mat_analyse$outflux[i] = outfl[ which(names(outfl)==mat_analyse$lieu[i]) ]
  }
}

for (i in 1:nrow(carac_f)){
  if(any(names(infl)==carac_f$lieu[i])){
    carac_f$influx[i] = infl[ which(names(infl)==carac_f$lieu[i]) ]
    carac_f$outflux[i] = outfl[ which(names(outfl)==carac_f$lieu[i]) ]
  }
}
#####

######################################################################
# Figure 5:

load(file="C:/Users/Jonathan/Desktop/Documents resultats/mat_eff_seed_carac.rdata")
par(mfrow = c(1, 2))

# For seed farm category:
seed_cat_analyse = as.data.frame(matrix(NA, nrow=5*9, ncol=8))
colnames(seed_cat_analyse) = c("seed_carac","contaminated_category","PFC","PFC_inf","PFC_sup","prev_slaugh","prev_slaugh_inf","prev_slaugh_sup")
seed_cat_analyse$seed_carac <- factor(seed_cat_analyse$seed_carac, levels = c("SEL","MU","FF","FA","FPW"))
seed_cat_analyse$contaminated_category <- factor(seed_cat_analyse$contaminated_category, levels = c("SEL","MU","FF","FA","FPW","PW","PWF","FI","total"))
iter = 1
for(seed_cat in c("SEL","MU","FF","FA","FPW")){
  for(cont_cat in c("total","SEL","MU","FF","FA","FPW","PW","PWF","FI")){
    seed_cat_analyse$seed_carac[iter] = seed_cat
    seed_cat_analyse$contaminated_category[iter] = cont_cat
    
    if(cont_cat == "total"){
      res_test_PFC = t.test(mat_eff_seed_carac$numb_cont_farms[which(mat_eff_seed_carac$site == seed_cat)])
      seed_cat_analyse$PFC_inf[iter] = res_test_PFC$conf.int[1]
      seed_cat_analyse$PFC_sup[iter] = res_test_PFC$conf.int[2]
      
      res_test_prevsl = t.test(mat_eff_seed_carac$prev_slaugh[which(mat_eff_seed_carac$site == seed_cat)])
      seed_cat_analyse$prev_slaugh[iter] = res_test_prevsl$estimate
      seed_cat_analyse$prev_slaugh_inf[iter] = res_test_prevsl$conf.int[1]
      seed_cat_analyse$prev_slaugh_sup[iter] = res_test_prevsl$conf.int[2]
    }else{
      seed_cat_analyse$PFC[iter] = mean(mat_eff_seed_carac[which(mat_eff_seed_carac$site == seed_cat),cont_cat])
    }
    
    iter = iter + 1
  }
}

# For seed farm outdegree:

mat_eff_seed_carac$cat_outdegree[which(mat_eff_seed_carac$cat_outdegree %in% c("[11;25]", "> 25"))] = ">= 11"

seed_outd_analyse = as.data.frame(matrix(NA, nrow=2*9, ncol=8))
colnames(seed_outd_analyse) = c("seed_carac","contaminated_category","PFC","PFC_inf","PFC_sup","prev_slaugh","prev_slaugh_inf","prev_slaugh_sup")
seed_outd_analyse$seed_carac <- factor(seed_outd_analyse$seed_carac, levels = c("< 11",">= 11"))
seed_outd_analyse$contaminated_category <- factor(seed_outd_analyse$contaminated_category, levels = c("SEL","MU","FF","FA","FPW","PW","PWF","FI","total"))
iter = 1
for(seed_outdeg in c("< 11",">= 11")){
  for(cont_cat in c("total","SEL","MU","FF","FA","FPW","PW","PWF","FI")){
    seed_outd_analyse$seed_carac[iter] = seed_outdeg
    seed_outd_analyse$contaminated_category[iter] = cont_cat
    
    if(cont_cat == "total"){
      res_test_PFC = t.test(mat_eff_seed_carac$numb_cont_farms[which(mat_eff_seed_carac$cat_outdegree == seed_outdeg)])
      seed_outd_analyse$PFC_inf[iter] = res_test_PFC$conf.int[1]
      seed_outd_analyse$PFC_sup[iter] = res_test_PFC$conf.int[2]
      
      res_test_prevsl = t.test(mat_eff_seed_carac$prev_slaugh[which(mat_eff_seed_carac$cat_outdegree == seed_outdeg)])
      seed_outd_analyse$prev_slaugh[iter] = res_test_prevsl$estimate
      seed_outd_analyse$prev_slaugh_inf[iter] = res_test_prevsl$conf.int[1]
      seed_outd_analyse$prev_slaugh_sup[iter] = res_test_prevsl$conf.int[2]
    }else{
      seed_outd_analyse$PFC[iter] = mean(mat_eff_seed_carac[which(mat_eff_seed_carac$cat_outdegree == seed_outdeg),cont_cat])
    }
    
    iter = iter + 1
  }
}

library(ggplot2)
library(gridExtra)

p1 = ggplot(data = seed_cat_analyse, aes(x = seed_carac))
p1 = p1 + xlab("Seed farm's category")
p1 = p1 + ggtitle("Effect of the seed farm's category")
p1 = p1 + ylab("Percentage of farms contaminated after 1 year (%)")
p1 = p1 + geom_bar(aes(y = PFC, fill = contaminated_category), stat="identity", position="stack")
p1 = p1 + geom_errorbar(aes(ymin = PFC_inf, ymax = PFC_sup), width = 0.1, size = 0.8)
p1 = p1 + geom_errorbar(aes(x = as.numeric(seed_carac), ymin = prev_slaugh_inf*130, ymax = prev_slaugh_sup*130), width = 0.1, color = "red", size = 0.8)
p1 = p1 + geom_point(aes(x = as.numeric(seed_carac), y = prev_slaugh*130), shape = 23,  size = 3, fill = "red", color = "black")
p1 = p1 + scale_y_continuous(sec.axis = sec_axis(~./130, name = "Prevalence in pigs heading for slaughterhouse after 1 year (%)"))
p1 = p1 + theme(axis.text.y.right = element_text(color="red"), axis.title.y.right = element_text(color="red"))
p1 = p1 + scale_fill_discrete(name = "Category of the contaminated farms:",
                              labels = c("Nucleus (SEL)","Multiplier (MU)","Farrowing-to-Finishing (FF)","Farrowing (FA)","Farrowing-Post-Weaning (FPW)","Post-Weaning (PW)","Post-Weaning-Finishing (PWF)","Finishing (FI)"))
p1 = p1 + theme(legend.position="bottom")

p2 = ggplot(data = seed_outd_analyse, aes(x = seed_carac))
p2 = p2 + xlab("Seed farm's outdegree")
p2 = p2 + ggtitle("Effect of the seed farm's outdegree")
p2 = p2 + ylab("Percentage of farms contaminated after 1 year (%)")
p2 = p2 + geom_bar(aes(y = PFC, fill = contaminated_category), stat="identity", position="stack")
p2 = p2 + geom_errorbar(aes(ymin = PFC_inf, ymax = PFC_sup), width = 0.1, size = 0.8)
p2 = p2 + geom_errorbar(aes(x = as.numeric(seed_carac), ymin = prev_slaugh_inf*130, ymax = prev_slaugh_sup*130), width = 0.1, color = "red", size = 0.8)
p2 = p2 + geom_point(aes(x = as.numeric(seed_carac), y = prev_slaugh*130), shape = 23,  size = 3, fill = "red", color = "black")
p2 = p2 + scale_y_continuous(sec.axis = sec_axis(~./130, name = "Prevalence in pigs heading for slaughterhouse after 1 year (%)"))
p2 = p2 + theme(axis.text.y.right = element_text(color="red"), axis.title.y.right = element_text(color="red"))
p2 = p2 + scale_fill_discrete(name = "Category of the contaminated farms:",
                              labels = c("Nucleus (SEL)","Multiplier (MU)","Farrowing-to-Finishing (FF)","Farrowing (FA)","Farrowing-Post-Weaning (FPW)","Post-Weaning (PW)","Post-Weaning-Finishing (PWF)","Finishing (FI)"))

count_cat = as.data.frame(t(table(factor(farms$site))))[,-1]
colnames(count_cat) = c("site", "number")
count_cat$site <- factor(count_cat$site, levels = c("SEL","MU","FF","FA","FPW","PW","PWF","FI"))
count_cat$perc = count_cat$number *100 / sum(count_cat$number)

p3 = ggplot(data = count_cat)
p3 = p3 + geom_bar(aes(x = "",y = number, fill = site), stat="identity", width=1, colour="black", size=0.01)
p3 = p3 + coord_polar("y")
p3 = p3+ scale_fill_discrete(name = "Category of the contaminated farms:",
                             labels = c("Nucleus (SEL)","Multiplier (MU)","Farrowing-to-Finishing (FF)","Farrowing (FA)","Farrowing-Post-Weaning (FPW)","Post-Weaning (PW)","Post-Weaning-Finishing (PWF)","Finishing (FI)"))

gl = list(ggplotGrob(p1 + theme(legend.position="none")), ggplotGrob(p2 + theme(legend.position="none")), ggplotGrob(p3 + theme(legend.position="right")))
grid.arrange(
  grobs = gl,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2),
                        c(1, 2),
                        c(3,3))
)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

p4 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               p3 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))


# To compare proportions of contaminated categories VS proportions in database:
#####
# Seed farm outdegree
tab1 = seed_outd_analyse[which(seed_outd_analyse$contaminated_category != "total"), c("seed_carac","contaminated_category","PFC")]
rownames(tab1) = 1:nrow(tab1)
tab1$perc_PFC = NA
tab1$perc_real = NA
for (outd in levels(tab1$seed_carac)){
  for (cat_contam in c("SEL","MU","FF","FA","FPW","PW","PWF","FI")){
    i = which((tab1$seed_carac == outd) & (tab1$contaminated_category == cat_contam))
    tab1$perc_PFC[i] = tab1$PFC[i] *100 / sum(tab1$PFC[tab1$seed_carac == outd], na.rm = T)
    
    tab1$perc_real[i] = count_cat$perc[which(count_cat$site == cat_contam)]
  }
}
tab1$ratio = tab1$perc_PFC / tab1$perc_real
tab1[which(tab1$ratio>2 | tab1$ratio<0.5),]

# Seed farm category
tab2 = seed_cat_analyse[which(seed_cat_analyse$contaminated_category != "total"), c("seed_carac","contaminated_category","PFC")]
rownames(tab2) = 1:nrow(tab2)
tab2$perc_PFC = NA
tab2$perc_real = NA
for (seed_cat in levels(tab2$seed_carac)){
  for (cat_contam in c("SEL","MU","FF","FA","FPW","PW","PWF","FI")){
    i = which((tab2$seed_carac == seed_cat) & (tab2$contaminated_category == cat_contam))
    tab2$perc_PFC[i] = tab2$PFC[i] *100 / sum(tab2$PFC[tab2$seed_carac == seed_cat], na.rm = T)
    
    tab2$perc_real[i] = count_cat$perc[which(count_cat$site == cat_contam)]
  }
}
tab2$ratio = tab2$perc_PFC / tab2$perc_real
tab2[which(tab2$ratio>2 | tab2$ratio<0.5),]

#####

######################################################################
# Analysis:

str(mat_analyse)

mat_analyse$categ=F
mat_analyse$categ[which(mat_analyse$site%in%c("SEL","MU"))]=T

full = lm(numb_cont_farms ~ ., data = mat_analyse[,c((3:13),17)])
null = lm(numb_cont_farms ~ 1, data = mat_analyse[,c((3:13),17)])
final_mod = step(null,  scope=list(lower=null, upper=full), direction="both", criterion = "AIC")
summary(final_mod)

# We assess if any variable's VIF is > 10
library(car)
vif(final_mod)

# ==> "size_bre" is removed and stepwise done again:
full_2 = lm(numb_cont_farms ~ size_pro+outdegree+indegree+betweeness+closeness+coreness+eigcen+outflux+influx+categ, data = mat_analyse[,c((3:13),17)])
final_mod_2 = step(null,  scope=list(lower=null, upper=full_2), direction="both", criterion = "AIC")
summary(final_mod_2)
vif(final_mod_2)




