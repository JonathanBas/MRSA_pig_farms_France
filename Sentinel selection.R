# Sentinel selection (using results from Scenario 2):

library(igraph)

# Priority lists (methods for selecting sentinels):
#####
senti_priority_outdegree = carac_f$lieu[order(carac_f$outdegree, decreasing=T)]
senti_priority_outdegree = senti_priority_outdegree[1:121]

senti_priority_indegree = carac_f$lieu[order(carac_f$indegree, decreasing=T)]
senti_priority_indegree = senti_priority_indegree[1:121]

senti_priority_betweeness = carac_f$lieu[order(carac_f$betweeness, decreasing=T)]
senti_priority_betweeness = senti_priority_betweeness[1:121]

senti_priority_eigcen = carac_f$lieu[order(carac_f$eigcen, decreasing=T)]
senti_priority_eigcen = senti_priority_eigcen[1:121]

senti_priority_closeness = carac_f$lieu[order(carac_f$closeness, decreasing=T)]
senti_priority_closeness = senti_priority_closeness[1:121]

senti_priority_coreness = carac_f$lieu[order(carac_f$coreness, decreasing=T)]
senti_priority_coreness = senti_priority_coreness[1:121]

senti_priority_influx = carac_f$lieu[order(carac_f$influx, decreasing=T)]
senti_priority_influx = senti_priority_influx[1:121]

senti_priority_outflux = carac_f$lieu[order(carac_f$outflux, decreasing=T)]
senti_priority_outflux = senti_priority_outflux[1:121]

# Priority list based on this paper from Holme: "Objective measures for sentinel surveillance in network epidemiology":
# Alternates farms that are most often contaminated and farms that are contaminated as early as possible
perc_detec_first = carac_f$lieu[order(carac_f$times_cont, decreasing=T)]
t_of_detec_first = carac_f$lieu[order(carac_f$t_of_cont/carac_f$times_cont, decreasing=F)]
senti_priority_altern = c(rbind(perc_detec_first, t_of_detec_first))
senti_priority_altern = senti_priority_altern[1:121]
rm(perc_detec_first)
rm(t_of_detec_first)
#####

# "Invasion paths" method priority list based on Bajardi et al. (2012) and Schirdewahn et al. (2017)
# From the ICSN matrix (based on Jaccard index), function "algo_cluster" gives a list of k sentinel farms
# Details are in SM3
#####
mat_analyse_2 = cbind(mat_analyse, as.data.frame(matrix(0,nrow=nrow(mat_analyse),ncol=nrow(carac_f))))
colnames(mat_analyse_2[,12:10553]) = as.character(carac_f$lieu)

inf_path = matrix(0, nrow=nrow(mat_analyse), ncol=nrow(carac_f))
rownames(inf_path) = mat_analyse$lieu
colnames(inf_path) = carac_f$lieu

for (f_start in mat_analyse_2$lieu){
  Categ=farms_repro$site[which(farms_repro$lieu==f_start)]
  
  for(simu in 1:nbsimu_scenario_2){
    load(paste0("~/",f_start,"/Results_",Categ,"_f_",f_start,"_simu_",simu,".rdata"))
    
    inf_net = colSums(results[[1]][,"PR",,])
    farms_nonnuls = names(which(rowSums(inf_net)!=0))
    farms_nonnuls = farms_nonnuls[which(farms_nonnuls!=f_start)]
    
    inf_path[f_start, farms_nonnuls] = inf_path[f_start, farms_nonnuls] + 1
    
    for(far in farms_nonnuls){
      nfarm=which(farms$lieu==far)
      mat_analyse_2[f_start, 11+nfarm] = mat_analyse_2[f_start, 11+nfarm] + (as.numeric(inf_net[far,20] != 0)/50)
    }
  }
}

inf_path_2 = inf_path[-which(mat_analyse$numb_cont_farms<5),]

jac_ind = matrix(0,nrow(inf_path_2),nrow(inf_path_2))
colnames(jac_ind) = rownames(jac_ind) = rownames(inf_path_2)

for (i in 1:nrow(inf_path_2)){
  print(i)
  
  for (j in 1:(i-1)){
    inf_i = colnames(inf_path_2)[which(inf_path_2[i,] != 0)]
    inf_j = colnames(inf_path_2)[which(inf_path_2[j,] != 0)]
    
    if ( (length(inf_i)!=0) & (length(inf_j)!=0) ){
      jac_ind[i,j] = length(intersect(inf_i,inf_j)) / length(union(inf_i,inf_j))
    }
  }
}

size_cluster = data.frame(filter=seq(0,0.1050,0.0025), size=0)
for (i in seq(0,0.1050,0.0005)){
  jac_ind2 = jac_ind
  jac_ind2[which(jac_ind < i)] = 0
  icsn = graph_from_adjacency_matrix(adjmatrix=jac_ind2, mode="max", weighted=T)
  compon = components(graph = icsn, mode="weak")
  size_cluster$size [which(size_cluster$filter == i)] = compon$no
}
rm(mat_analyse_2)

compon_all = list()

# For 30 sentinels
filt = min(size_cluster$filter[which(size_cluster$size == 30)])
jac_ind2 = jac_ind
jac_ind2[which(jac_ind < filt)] = 0
icsn = graph_from_adjacency_matrix(adjmatrix=jac_ind2, mode="max", weighted=T)
compon_all[["30"]] = components(graph = icsn, mode="weak")

# For 60 sentinels
filt = min(size_cluster$filter[which(size_cluster$size == 60)])
jac_ind2 = jac_ind
jac_ind2[which(jac_ind < filt)] = 0
icsn = graph_from_adjacency_matrix(adjmatrix=jac_ind2, mode="max", weighted=T)
compon_all[["60"]] = components(graph = icsn, mode="weak")

# For 120 sentinels
filt = min(size_cluster$filter[which(size_cluster$size == 120)])
jac_ind2 = jac_ind
jac_ind2[which(jac_ind < filt)] = 0
icsn = graph_from_adjacency_matrix(adjmatrix=jac_ind2, mode="max", weighted=T)
compon_all[["120"]] = components(graph = icsn, mode="weak")

# /!\ k needs to be in the possible values in 'size_cluster$size'
algo_cluster = function(k){
  compon = compon_all[[as.character(k)]]
  
  senti_priority_infpath = rep(0,k)
  
  clust_bigger_1 = which(compon$csize > 1)
  clust_equal_1 = (1:k)[-clust_bigger_1]
  
  for(elem in clust_bigger_1){
    possibilities = names(compon$membership)[which(compon$membership == elem)]
    senti_priority_infpath[elem] = sample(x = possibilities, size = 1)
  }
  for(elem in clust_equal_1){
    senti_priority_infpath[elem] = names(compon$membership)[which(compon$membership == elem)]
  }
  
  return(senti_priority_infpath)
}

#####

# We want to compare different priority lists and random for sentinel selection, based on 3 efficiency criteria
# We consider epidemics starting from all seed farms
n_rand_samp = 100 # Number of random samples drawn
n_infpath_samp = 100 # Number of samples drawn for the "Invasion paths" method

# Fct that takes a sample of sentinel farms and returns efficiency indicators, for 1 repetition and 1 seed
fct_detec = function(samp, m_output, farms_nonnuls){
  epi_detect_iter = 0
  t_detect_iter = NA
  nb_f_detect_iter = NA
  
  if(any(samp%in%farms_nonnuls)){
    epi_detect_iter = 1 # One farm from sample was contaminated
    farms_nonnuls_samp = intersect(samp,farms_nonnuls) # Farms from sample that were contaminated
    
    pop_weeks_samp=colSums(m_output[,"PR",farms_nonnuls_samp,], dims = min(2,length(farms_nonnuls_samp)) )

    t_detect_iter = min(which(pop_weeks_samp!=0)) # Minimum week at which farm from sample was contaminated
    
    pop_weeks=colSums(m_output[,"PR",farms_nonnuls,]) # All contaminated farms (not only sample)
    
    if(length(farms_nonnuls)==1){
      nb_f_detect_iter = as.numeric(pop_weeks[t_detect_iter]!=0)
    }else{
      nb_f_detect_iter = length(which(pop_weeks[,t_detect_iter]!=0))
    }
    nb_f_detect_iter = nb_f_detect_iter*100 / nrow(carac_f)
  }
  return(c(epi_detect_iter, t_detect_iter, nb_f_detect_iter))
}

# Function that returns efficiency indicators for 1 seed farm (all repetitions) for all priority lists
scen_intro_analysis = function(n_f){
  val_num_senti = c("30", "60", "120")
  
  f = farms_repro$lieu[n_f]
  Categ=farms_repro$site[which(farms_repro$lieu==f)]
  
  #####
  # Matrix for all methods (priority lists) excepted "Invasion paths" method
  detect_score = array(data=NA, dim=c(nbsimu_scenario_2, 9, 3))
  dimnames(detect_score)[[2]] = c("outdegree", "indegree", "betweeness", "closeness", "coreness", "eigcen", "outflux", "influx", "altern")
  dimnames(detect_score)[[3]] = c("epi_detec", "t_detec", "n_f_detec")
  detect_score = rep(list(detect_score), length(val_num_senti))
  names(detect_score) = val_num_senti
  
  # Matrix for random selection
  detect_rand = as.data.frame(matrix(data=NA, nrow=n_rand_samp*nbsimu_scenario_2, ncol=3))
  colnames(detect_rand) = c("epi_detec", "t_detec", "n_f_detec")
  detect_rand = rep(list(detect_rand), length(val_num_senti))
  names(detect_rand) = val_num_senti
  
  # Matrix for "Invasion paths" method
  detect_infpath = as.data.frame(matrix(data=NA, nrow=n_infpath_samp*nbsimu_scenario_2, ncol=3))
  colnames(detect_infpath) = c("epi_detec", "t_detec", "n_f_detec")
  detect_infpath = rep(list(detect_infpath), length(val_num_senti))
  names(detect_infpath) = val_num_senti
  #####
  
  for(simu in 1:nbsimu_scenario_2){
    load(paste0("~/",f,"/Results_",Categ,"_f_",f,"_simu_",simu,".rdata"))
    m_output=results[[1]]
    farms_nonnuls=names(which(rowSums(colSums(m_output[,"PR",,]))!=0))
    farms_nonnuls=farms_nonnuls[which(farms_nonnuls!=f)]
    
    for (k in val_num_senti){
      detect_score[[k]][simu,"outdegree",] = fct_detec(samp = senti_priority_outdegree[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"indegree",] = fct_detec(samp = senti_priority_indegree[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"betweeness",] = fct_detec(samp = senti_priority_betweeness[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"eigcen",] = fct_detec(samp = senti_priority_eigcen[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"closeness",] = fct_detec(samp = senti_priority_closeness[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"coreness",] = fct_detec(samp = senti_priority_coreness[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"outflux",] = fct_detec(samp = senti_priority_outflux[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"influx",] = fct_detec(samp = senti_priority_influx[1:as.numeric(k)], m_output, farms_nonnuls)
      detect_score[[k]][simu,"altern",] = fct_detec(samp = senti_priority_altern[1:as.numeric(k)], m_output, farms_nonnuls)
      
      for(rand_samp in 1:n_rand_samp){
        detect_rand[[k]][(simu-1)*n_rand_samp + rand_samp, ] = fct_detec(samp = sample(sample(carac_f$lieu, as.numeric(k))), m_output, farms_nonnuls)
      }
      
      for(infpath_samp in 1:n_infpath_samp){
        senti_priority_infpath = algo_cluster(as.numeric(k))
        detect_infpath[[k]][(simu-1)*n_infpath_samp + infpath_samp, ] = fct_detec(samp = senti_priority_infpath, m_output, farms_nonnuls)
      }
    }
  }
  save(detect_score, file=paste0("~/detect_score_",f,".rdata"))
  save(detect_rand, file=paste0("~/detect_rand_",f,".rdata"))
  save(detect_infpath, file=paste0("~/detect_infpath_",f,".rdata"))
}

######################################################################

# Parallelization:
library(parallel)
cl<-makeCluster(7)
clusterExport(cl,c("carac_f","fct_detec","farms_repro","nbsimu_scenario_2","senti_priority_outdegree","senti_priority_indegree","senti_priority_betweeness","senti_priority_eigcen","senti_priority_closeness","senti_priority_coreness","senti_priority_outflux","senti_priority_influx","senti_priority_altern","algo_cluster","compon_all"))
parLapply(cl,X=1:nrow(mat_analyse),fun=scen_intro_analysis)
stopCluster(cl)

######################################################################
# Analysis:

library(Hmisc)

val_num_senti = c("30", "60", "120")

# "indic_eff_surv" = final matrix for surveillance results
#####
indic_eff_surv = data.frame(matrix(NA,3*11*3,6))
colnames(indic_eff_surv) = c("indic", "numb_sent", "score", "val", "inf", "sup")
indic_eff_surv$numb_sent <- factor(indic_eff_surv$numb_sent, levels = c("30","60","120"))

iter = 1
for(k in val_num_senti){
  for(indi in c("epi_detec", "t_detec", "n_f_detec")){
    for(sco in c("outdegree", "indegree", "betweeness", "closeness", "coreness", "eigcen", "outflux", "influx", "altern", "rand", "infpath")){
      indic_eff_surv$numb_sent[iter] = k
      indic_eff_surv$indic[iter] = indi
      indic_eff_surv$score[iter] = sco
      iter = iter + 1
    }
  }
}
#####

# "detect_score_tot", "detect_rand_tot" and "detect_infpath_tot" are similar to "detect_score", "detect_rand" and "detect_infpath", but add for every seed
#####
detect_score_tot = array(data=NA, dim=c(nbsimu_scenario_2 * nrow(farms_repro), 9, 3))
dimnames(detect_score_tot)[[2]] = c("outdegree", "indegree", "betweeness", "closeness", "coreness", "eigcen", "outflux", "influx", "altern")
dimnames(detect_score_tot)[[3]] = c("epi_detec", "t_detec", "n_f_detec")
detect_score_tot = rep(list(detect_score_tot), length(val_num_senti))
names(detect_score_tot) = val_num_senti

detect_rand_tot = as.data.frame(matrix(data=NA, nrow=n_rand_samp*nbsimu_scenario_2*nrow(farms_repro), ncol=3))
colnames(detect_rand_tot) = c("epi_detec", "t_detec", "n_f_detec")
detect_rand_tot = rep(list(detect_rand_tot), length(val_num_senti))
names(detect_rand_tot) = val_num_senti

detect_infpath_tot = as.data.frame(matrix(data=NA, nrow=n_infpath_samp*nbsimu_scenario_2*nrow(farms_repro), ncol=3))
colnames(detect_infpath_tot) = c("epi_detec", "t_detec", "n_f_detec")
detect_infpath_tot = rep(list(detect_infpath_tot), length(val_num_senti))
names(detect_infpath_tot) = val_num_senti
#####

for(f in farms_repro$lieu){
  iter = which(farms_repro$lieu==f)
  print(paste0("Analyse de la ferme: ",f," (",iter,"/",length(farms_repro$lieu),")"))
  
  load(file=paste0("~/detect_score_",f,".rdata"))
  load(file=paste0("~/detect_rand_",f,".rdata"))
  load(file=paste0("~/detect_infpath_",f,".rdata"))
  
  for(k in val_num_senti){
    detect_score_tot[[k]][(((iter-1)*nbsimu_scenario_2 + 1):((iter-1)*nbsimu_scenario_2 + nbsimu_scenario_2)),,] = detect_score[[k]]
    detect_rand_tot[[k]][(((iter-1)*n_rand_samp*nbsimu_scenario_2 + 1):((iter-1)*n_rand_samp*nbsimu_scenario_2 + n_rand_samp*nbsimu_scenario_2)),] = detect_rand[[k]]
    detect_infpath_tot[[k]][(((iter-1)*n_infpath_samp*nbsimu_scenario_2 + 1):((iter-1)*n_infpath_samp*nbsimu_scenario_2 + n_infpath_samp*nbsimu_scenario_2)),] = detect_infpath[[k]]
  }
}

# For all methods (priority lists) excepted "Invasion paths" method:
for(k in val_num_senti){
  for(sco in c("outdegree", "indegree", "betweeness", "closeness", "coreness", "eigcen", "outflux", "influx", "altern")){
    epi_d = detect_score_tot[[k]][, sco, "epi_detec"]
    indic_eff_surv[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "epi_detec") & (indic_eff_surv$score == sco)), c("val","inf","sup")] = 100 * binconf(sum(epi_d, na.rm=T), length(epi_d))
    
    t_d = detect_score_tot[[k]][, sco, "t_detec"]
    val_t_d = t.test(t_d)
    indic_eff_surv$val[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == sco))] = mean(t_d, na.rm=T)
    indic_eff_surv$inf[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == sco))] = quantile(t_d, 0.025, na.rm=T)
    indic_eff_surv$sup[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == sco))] = quantile(t_d, 0.975, na.rm=T)
    
    n_f_d = detect_score_tot[[k]][, sco, "n_f_detec"]
    val_n_f_d = t.test(n_f_d)
    indic_eff_surv$val[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == sco))] = mean(n_f_d, na.rm=T)
    indic_eff_surv$inf[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == sco))] = quantile(n_f_d, 0.025, na.rm=T)
    indic_eff_surv$sup[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == sco))] = quantile(n_f_d, 0.975, na.rm=T)
  }
}

# Random:
for(k in val_num_senti){
  epi_d = detect_rand_tot[[k]][, "epi_detec"]
  indic_eff_surv[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "epi_detec") & (indic_eff_surv$score == "rand")), c("val","inf","sup")] = 100 * binconf(sum(epi_d, na.rm=T), length(epi_d))
  
  t_d = detect_rand_tot[[k]][, "t_detec"]
  val_t_d = t.test(t_d)
  indic_eff_surv$val[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == "rand"))] = mean(t_d, na.rm=T)
  indic_eff_surv$inf[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == "rand"))] = quantile(t_d, 0.025, na.rm=T)
  indic_eff_surv$sup[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == "rand"))] = quantile(t_d, 0.975, na.rm=T)
  
  n_f_d = detect_rand_tot[[k]][, "n_f_detec"]
  val_n_f_d = t.test(n_f_d)
  indic_eff_surv$val[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == "rand"))] = mean(n_f_d, na.rm=T)
  indic_eff_surv$inf[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == "rand"))] = quantile(n_f_d, 0.025, na.rm=T)
  indic_eff_surv$sup[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == "rand"))] = quantile(n_f_d, 0.975, na.rm=T)
}

# "Invasion paths" method:
for(k in val_num_senti){
  epi_d = detect_infpath_tot[[k]][, "epi_detec"]
  indic_eff_surv[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "epi_detec") & (indic_eff_surv$score == "infpath")), c("val","inf","sup")] = 100 * binconf(sum(epi_d, na.rm=T), length(epi_d))
  
  t_d = detect_infpath_tot[[k]][, "t_detec"]
  val_t_d = t.test(t_d)
  indic_eff_surv$val[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == "infpath"))] = mean(t_d, na.rm=T)
  indic_eff_surv$inf[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == "infpath"))] = quantile(t_d, 0.025, na.rm=T)
  indic_eff_surv$sup[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "t_detec") & (indic_eff_surv$score == "infpath"))] = quantile(t_d, 0.975, na.rm=T)
  
  n_f_d = detect_infpath_tot[[k]][, "n_f_detec"]
  val_n_f_d = t.test(n_f_d)
  indic_eff_surv$val[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == "infpath"))] = mean(n_f_d, na.rm=T)
  indic_eff_surv$inf[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == "infpath"))] = quantile(n_f_d, 0.025, na.rm=T)
  indic_eff_surv$sup[which((indic_eff_surv$numb_sent == k) & (indic_eff_surv$indic == "n_f_detec") & (indic_eff_surv$score == "infpath"))] = quantile(n_f_d, 0.975, na.rm=T)
}

######################################################################
# Plot:
library(ggplot2)

which_indicator = "epi_detec" # or "t_detec" or "n_f_detec"

if(which_indicator == "epi_detec"){
  
  p = ggplot(data = indic_eff_surv[which(indic_eff_surv$indic =="epi_detec"),], aes(x = numb_sent, y = val, fill = score))
  p = p + ggtitle("Percentage of detected incursions")
  p = p + ylab("Percentage of detected incursions (%)")
  
}else if(which_indicator == "t_detec"){
  
  p = ggplot(data = indic_eff_surv[which(indic_eff_surv$indic =="t_detec"),], aes(x = numb_sent, y = val, fill = score))
  p = p + ggtitle("Time before detection")
  p = p + ylab("Time before detection (weeks)")
  
}else if(which_indicator == "n_f_detec"){

  p = ggplot(data = indic_eff_surv[which(indic_eff_surv$indic =="n_f_detec"),], aes(x = numb_sent, y = val, fill = score))
  p = p + ggtitle("Percentage of contaminated farms at detection")
  p = p + ylab("Percentage of contaminated farms at detection (%)")

}

p = p + xlab("Number of sentinel farms")
p = p + geom_bar(stat="identity", position="dodge", colour="black")
p = p + geom_errorbar(aes(ymin = inf, ymax = sup), position=position_dodge(0.9), width = 0.3, size = 0.8)
p = p + scale_fill_brewer(palette = "Paired",
                          name = "Method for sentinels selection:",
                          labels = c("Alternated method", "Highest betweeness centrality", "Highest closeness centrality", "Highest coreness centrality", "Highest eigenvector centrality", "Highest indegree centrality", "Highest influx centrality", "Invasion paths method", "Highest outdegree centrality", "Highest outflux centrality", "Random"))
p = p + geom_vline(xintercept = c(1.5,2.5), linetype = "dashed", size = 1)
p = p + theme(legend.position="bottom")
p

