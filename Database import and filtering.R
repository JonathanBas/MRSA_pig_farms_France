rm(list=ls(all=TRUE))
library(mc2d)

# Import
#####
load("~/siteDesc.rdata")
siteDesc$lieu=as.character(siteDesc$lieu)
table(factor(siteDesc$site)) # Categories

# Removed from "farms":
# 1) abattoir, TR, BS, SP, WB
# 2) categories SEL, MU, FA with null breeding herd size, and therefore null total herd size in the model
farms=siteDesc[which(!(siteDesc$site%in%c("abattoir","TR","BS","WB","SP"))),c(1,3,4,5,6,9)]
farms=farms[-which((farms$site%in%c("SEL","MU","FA"))&(farms$repro==0)),]
table(factor(farms$site))

load("~/nbAnimalMov.rdata")

movb=nbAnimalMov$Breeders[farms$lieu,farms$lieu]
movp=nbAnimalMov$piglets[farms$lieu,farms$lieu]+nbAnimalMov$weaners[farms$lieu,farms$lieu]+nbAnimalMov$finishers[farms$lieu,farms$lieu]
mov=movb+movp

# Removed from "farms": no outgoing nor ingoing movement
farms=farms[which(!(rowSums(mov)==0&colSums(mov)==0)),]
table(factor(farms$site))

movb=movb[farms$lieu,farms$lieu]
movp=movp[farms$lieu,farms$lieu]

rm(nbAnimalMov)
rm(siteDesc)
rm(mov)
#####

# farms: categories (FF, FA, FI, FPW, PW, PWF, SEL, MU), not isolated, with non null breeding herd for SEL, MU et FA
# movb: movements of breeding pigs
# movp: movements of fattening pigs

# Categories:
#####

# 289/3351 of FI import from FA.
# Among them, 27 export to FI, PWF, FF (and for 1, PW)
# ==> Renamed as PWF
farms_FA=farms[which(farms$site=="FA"),]
farms_FI=farms[which(farms$site=="FI"),]
FI_imp_FA=farms_FI$lieu[which(colSums(movp[as.character(farms_FA$lieu),as.character(farms_FI$lieu)])!=0)]
FI_to_PWF=names(which(rowSums(movp[as.character(FI_imp_FA),])!=0))
farms$site[which(farms$lieu%in%FI_to_PWF)]="PWF"
rm(FI_imp_FA)
rm(FI_to_PWF)
rm(farms_FA)
rm(farms_FI)

farms_A=farms[which(farms$site%in%c("FA","FPW","FF","MU","SEL")),]

farms_B=farms[which(farms$site%in%c("FPW","PW","PWF","FF","MU","SEL")),]

farms_MU=farms[which(farms$site=="MU"),]
farms_SEL=farms[which(farms$site=="SEL"),]
farms_PWF=farms[which(farms$site=="PWF"),]
farms_PW=farms[which(farms$site=="PW"),]
farms_PWF_FI=farms[which(farms$site%in%c("PWF","FI")),]
farms_SEL_MU_FF=farms[which(farms$site%in%c("SEL","MU","FF")),]
farms_SEL_MU_FF_FPW=farms[which(farms$site%in%c("SEL","MU","FF","FPW")),]

farms_repro=farms[which((farms$repro>=1)&(farms$site%in%c("MU","SEL","FA","FF","FPW"))),]
farms_nonrepro=farms[which(!(farms$lieu%in%farms_repro$lieu)),]

#####

# 2 movements periods for fattening pigs: A->B and B->C

# Integers matrices: movAB_non_prop and movBC_non_prop
# Proportion matrices (sum of rows = 1 or 0): movAB and movBC
# prop_export = Vector with proportion of a batch exported

#####

# 1) movAB
movAB_non_prop=movp
movAB_non_prop[,which(farms$site%in%c("FA","MU","SEL"))]=0 #FA (= strict A), MU and SEL receive no animal
movAB_non_prop[which(!(farms$site%in%c("FA"))),]=0 # Non-FA categories do not export at THIS step (4 weeks)

movAB=movAB_non_prop
for (i in 1:nrow(movAB_non_prop)){
  a=sum(movAB_non_prop[i,])
  if (a!=0){
    movAB[i,]=movAB_non_prop[i,]/a
  }else{
    if(any(movAB_non_prop[i,]!=0)){print("Erreur sur movAB_non_prop")}
  }
}

f_imp_AB=farms$lieu[which(colSums(movAB)!=0)]

table(factor(rowSums(movAB))) # Should be only 0 and 1
table(factor(farms$site[which(rowSums(movAB)==1)])) # Export
table(factor(farms$site[which(rowSums(movAB)==0)])) # Don't export

table(factor(farms$site[which(colSums(movAB)!=0)])) # Import
table(factor(farms$site[which(colSums(movAB)==0)])) # Don't import

# 2) movBC
movBC_non_prop=movp
movBC_non_prop[,which(farms$site%in%c("FA","FPW","PW","MU","SEL"))]=0 # FA, FPW, PW, MU, SEL receive no animal
movBC_non_prop[which(farms$site%in%c("FA","FI")),]=0 # FA, FI don't export at THIS step (4 weeks)

movBC=movBC_non_prop
for (i in 1:nrow(movBC_non_prop)){
  a=sum(movBC_non_prop[i,])
  if (a!=0){
    movBC[i,]=movBC_non_prop[i,]/a
  }else{
    if(any(movBC_non_prop[i,]!=0)){print("Erreur sur movBC_non_prop")}
  }
}

f_imp_BC=farms$lieu[which(colSums(movBC)!=0)]

table(factor(rowSums(movBC)))
table(factor(farms$site[which(rowSums(movBC)==1)]))
table(factor(farms$site[which(rowSums(movBC)==0)]))

table(factor(farms$site[which(colSums(movBC)!=0)]))
table(factor(farms$site[which(colSums(movBC)==0)]))

# 3) prop_export
prop_export=(rowSums(movAB_non_prop+movBC_non_prop)/13)/((farms$herdSize-farms$repro)/7)
prop_export[which(farms$site%in%c("FA","FPW","PW"))]=1 # FA, FPW, PW export all fattening pigs
prop_export[which(farms$site%in%c("FI"))]=0 # Value doesn't matter for FI because they never export
prop_export[which(farms$site%in%c("SEL","MU"))]=pmin(0.5,prop_export[which(farms$site%in%c("SEL","MU"))]) # Limited to 50%

f_parti=which(!(is.finite(prop_export)))
prop_export[f_parti]=1
prop_export[which(prop_export>1)]=1
print("Quantiles de prop_export:")
quantile(prop_export,probs=seq(0,1,0.05))
rm(f_parti)

table(factor(farms$site[which(prop_export==1)])) # Export 100%
table(factor(farms$site[which(!(prop_export%in%c(0,1)))])) # Export between 0% and 100% (excluded)
table(factor(farms$site[which(prop_export==0)])) # Export 0%
#####

# Breeding pigs:

# Integers matrix: mov_gilts_non_prop
# Proportion matrix (sum of columns = 1 ou 0): mov_gilts

#####
mov_gilts_non_prop=movb
mov_gilts_non_prop[,union(as.character(farms_nonrepro$lieu),farms$lieu[which(farms$site%in%c("PW","PWF","FI"))])]=0 # Receiving farms with no breeding pigs or not in ("MU","SEL","FA","FF","FPW") import nothing
mov_gilts_non_prop[farms$lieu[which(farms$site%in%c("FA","FPW","FF","PW","PWF","FI"))],]=0 # Non MU or SEL categories export no breeding pigs

mov_gilts=mov_gilts_non_prop
for (i in 1:ncol(mov_gilts_non_prop)){
  a=sum(mov_gilts_non_prop[,i])
  if (a!=0){
    mov_gilts[,i]=mov_gilts_non_prop[,i]/a
  }else{
    if(any(mov_gilts_non_prop[,i]!=0)){print("Erreur sur mov_gilts_non_prop")}
  }
}

# Self-renewal of sows in farms with breeding pigs but with no recorded ingoing gilts movements:
repro_without_imp=farms_repro$lieu[which(colSums(mov_gilts[,farms_repro$lieu])==0)]
for(i in repro_without_imp){mov_gilts[i,i]=1}
rm(repro_without_imp)

# Farms exporting gilts:
f_exp_gilts=farms$lieu[which(rowSums(mov_gilts[,farms_repro$lieu])!=0)]

table(factor(colSums(mov_gilts)))
table(factor(farms$site[which(colSums(mov_gilts)==1)]))
table(factor(farms$site[which(colSums(mov_gilts)==0)]))

table(factor(farms$site[which(rowSums(mov_gilts)!=0)])) # Include self-renewal
table(factor(farms$site[which(rowSums(mov_gilts)==0)]))
#####

rm(movb)
rm(movp)

######################################################################
# Functions and parameters:

epidemio=function(compart,contact_compart,beta){
  N=ncol(compart)
  
  # Number of PR animals in compartments in contact
  PR_contact=rep(0,N)
  P_tot=rep(0,N)
  for (co in contact_compart){
    PR_contact=PR_contact+co["PR",]
    P_tot=P_tot+colSums(co)
  }
  
  p1=beta*PR_contact/P_tot # Numeric approximation of probability
  p1[which(!(is.finite(p1)))]=0 # If null population
  p2=rep(VPR,N) # Numeric approximation of probability
  delta.1=rbinom(n=N,size=compart["PX",],prob=p1) # X to R
  delta.2=rbinom(n=N,size=compart["PR",],prob=p2) # R to X
  
  compart["PX",]=compart["PX",]-delta.1+delta.2
  compart["PR",]=compart["PR",]+delta.1-delta.2
  
  return(list(compart,delta.1))
}

epidemio_cont_meas=function(compart,contact_compart,beta,beta_in_farms_with_contmeas){
  N=ncol(compart)
  
  # Number of PR animals in compartments in contact
  PR_contact=rep(0,N)
  P_tot=rep(0,N)
  for (co in contact_compart){
    PR_contact=PR_contact+co["PR",]
    P_tot=P_tot+colSums(co)
  }
  
  p1=beta*PR_contact/P_tot # Numeric approximation of probability
  p1[which(!(is.finite(p1)))]=0 # If null population
  p2=rep(VPR,N) # Numeric approximation of probability
  delta.1=rbinom(n=N,size=compart["PX",],prob=p1) # X to R
  delta.2=rbinom(n=N,size=compart["PR",],prob=p2) # R to X
  
  compart["PX",]=compart["PX",]-delta.1+delta.2
  compart["PR",]=compart["PR",]+delta.1-delta.2
  
  # Specific to farms with control measures (different Beta: beta_in_farms_with_contmeas):
  
  Ncm=length(f_cont_meas)
  
  PR_contact_cm=rep(0,Ncm)
  P_tot_cm=rep(0,Ncm)
  for (co in contact_compart){
    PR_contact_cm=PR_contact_cm+co["PR",f_cont_meas]
    P_tot_cm=P_tot_cm+colSums(co[,f_cont_meas])
  }
  
  p1=beta_in_farms_with_contmeas*PR_contact_cm/P_tot_cm # Numeric approximation of probability
  p1[which(!(is.finite(p1)))]=0 # If null population
  p2=rep(VPR,Ncm) # Numeric approximation of probability
  delta.1cm=rbinom(n=Ncm,size=compart["PX",f_cont_meas],prob=p1) # X to R
  delta.2cm=rbinom(n=Ncm,size=compart["PR",f_cont_meas],prob=p2) # R to X
  
  compart["PX",f_cont_meas]=compart["PX",f_cont_meas]-delta.1cm+delta.2cm
  compart["PR",f_cont_meas]=compart["PR",f_cont_meas]+delta.1cm-delta.2cm
  
  return(list(compart,delta.1,delta.1cm))
}

# Function giving value of model's outputs at t = 52 weeks:
val_output=function(output,results,iter,simu){
  if(output=="n_f_cont"){
    m_output=results[[1]]
    t_fin=52
    val=sum(colSums(m_output[,"PR",,t_fin])!=0) *100/10542
  }
  if(output=="apop_prev"){
    # Prevalence in %:
    apop_prev = results[[1]]
    pop_tot = results[[3]]
    t_fin=52
    val = sum(apop_prev[,"PR",,t_fin]) *100 / sum(pop_tot[,t_fin])
  }
  if(output=="to_abat_prev"){
    # Prevalence in %:
    m_output=results[[4]]
    t_fin=52
    val = mean(colSums(results[[4]]["PR",,seq(t_fin-3,t_fin)],na.rm=T)) *100 / numb_pigs_slaughter[iter, simu]
  }
  return(val)
}

# Function giving value of model's outputs over time:
val_output_2=function(output,results){
  if(output=="n_f_cont"){
    m_output=results[[1]]
    t_fin=dim(m_output)[4]
    t1=1:t_fin
    val=rep(0,t_fin)
    for (t in t1){val[t]=sum(colSums(m_output[,"PR",,t])!=0)*100/10542}
  }
  if(output=="apop_prev"){
    # Prevalence in %:
    m_output=results[[1]]
    t_fin=dim(m_output)[4]
    pop_tot = results[[3]]
    val=colSums(m_output[,"PR",,],dims=2) *100 / colSums(pop_tot)
  }
  if(output=="to_abat_prev"){
    # Prevalence in %:
    m_output=results[[3]]
    t_fin=dim(m_output)[3]
    t1=seq(4,t_fin,4)
    val=rep(NA,t_fin)
    for (t in t1){val[t]=mean(colSums(results[[3]]["PR",,seq(t-3,t)],na.rm=T)) *100 / numb_pigs_slaughter[iter, simu]}
  }
  return(val)
}

load(file="~/Temps equilibre.rdata")
load(file = "~/numb_pigs_slaughter.rdata")

t_start=1
I=nrow(farms)

betaPPR=0.262 # Calibrated using scenario 1 (baseline)
VPR=1/4.3 # 1/[carriage duration] (baseline)

# Pig population on week t
apop=rep(list(matrix(0,nrow=2,ncol=I,dimnames=list(c("PX","PR"),as.character(farms$lieu)))),12)
names(apop)=c("G1","G2","G3","G4","LA","FA","PW1","PW2","F1","F2","F3","F4")

# Farms exporting pigs at step AB:
f_exportAB=Reduce(intersect,list(as.character(farms_A$lieu),names(which(prop_export!=0)),names(which(rowSums(movAB)!=0))))

# Farms exporting pigs at step BC:
f_exportBC=Reduce(intersect,list(as.character(farms_B$lieu),names(which(prop_export!=0)),names(which(rowSums(movBC)!=0))))

# Number of imported gilts in farms with breeding pigs:
n_imp_gilts=round(0.0319*farms$repro,0)
names(n_imp_gilts)=as.character(farms$lieu)

# Initialization sows:
#####
apop$G1["PX",as.character(farms_repro$lieu)]=pmax(1,round(0.2*farms_repro$repro,0))
apop$G2["PX",as.character(farms_repro$lieu)]=pmax(1,round(0.2*farms_repro$repro,0))
apop$G3["PX",as.character(farms_repro$lieu)]=pmax(1,round(0.2*farms_repro$repro,0))
apop$G4["PX",as.character(farms_repro$lieu)]=pmax(1,round(0.2*farms_repro$repro,0))
apop$LA["PX",as.character(farms_repro$lieu)]=pmax(1,round(0.2*farms_repro$repro,0))
#####

# Initialization fattening pigs:
# ==> "apop_initial"

# Farm age vector (randomly shifted):
agep_init=sample(x=1:4,size=I,replace=T)
names(agep_init)=as.character(farms$lieu)

# Farms without breeding pigs:
#####
apop$PW1["PX",as.character(farms_PWF$lieu)]=round((1/2)*farms_PWF$ps,0)
apop$PW2["PX",as.character(farms_PWF$lieu)]=round((1/2)*farms_PWF$ps,0)
apop$PW1["PX",as.character(farms_PW$lieu)]=round((1/2)*farms_PW$ps,0)
apop$PW2["PX",as.character(farms_PW$lieu)]=round((1/2)*farms_PW$ps,0)

apop$F1["PX",as.character(farms_PWF_FI$lieu)]=round((1/4)*farms_PWF_FI$eng,0)
apop$F2["PX",as.character(farms_PWF_FI$lieu)]=round((1/4)*farms_PWF_FI$eng,0)
apop$F3["PX",as.character(farms_PWF_FI$lieu)]=round((1/4)*farms_PWF_FI$eng,0)
apop$F4["PX",as.character(farms_PWF_FI$lieu)]=round((1/4)*farms_PWF_FI$eng,0)
#####

# Farms with breeding pigs:
#####
apop$FA["PX",as.character(farms_repro$lieu)]=round(11.8*0.2*farms_repro$repro,0)

apop$PW1["PX",as.character(farms_SEL_MU_FF_FPW$lieu)]=round((1/2)*farms_SEL_MU_FF_FPW$ps,0)
apop$PW2["PX",as.character(farms_SEL_MU_FF_FPW$lieu)]=round((1/2)*farms_SEL_MU_FF_FPW$ps,0)

apop$F1["PX",as.character(farms_SEL_MU_FF$lieu)]=round((1/4)*farms_SEL_MU_FF$eng,0)
apop$F2["PX",as.character(farms_SEL_MU_FF$lieu)]=round((1/4)*farms_SEL_MU_FF$eng,0)
apop$F3["PX",as.character(farms_SEL_MU_FF$lieu)]=round((1/4)*farms_SEL_MU_FF$eng,0)
apop$F4["PX",as.character(farms_SEL_MU_FF$lieu)]=round((1/4)*farms_SEL_MU_FF$eng,0)
#####

rm(list=c("farms_SEL","farms_MU","farms_PW","farms_PWF","farms_PWF_FI","farms_SEL_MU_FF","farms_SEL_MU_FF_FPW"))

# Model ran until the pig herd populations stabilize:
for (t in t_start:32){
  
  # Demographic dynamics:
  
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
    
    p_gilts_X=colSums(mov_gilts[f_exp_gilts_t,f_repro_chang]*prop_X_F4) # On multiplie mov_gilts colonne par colonne par le vecteur prop_X_F4 (prop_X_F4 ne change pas), puis on fait la somme de chaque colonne
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
    
    movBC_number=floor(movBC[f_chang_exportBC,f_imp_BC]*n_exportBC) # movBC(matrice proportion valable en tout t) mulitpliee colonne par colonne (pour les colonnes de fermes qui importent dans la matrice a BC) par n_exportBC (nombre de porcs exportes a ce tour ci par chaque ferme exportatrice)
    
    prop_X_BC=draw_leave_farm[,1]/(rowSums(draw_leave_farm))
    prop_X_BC[which(!(is.finite(prop_X_BC)))]=0 # Quand la population de l'elevage d'origine est nulle, cela cree un NaN ==> on considere que c'est zero
    prop_R_BC=draw_leave_farm[,2]/(rowSums(draw_leave_farm))
    prop_R_BC[which(!(is.finite(prop_R_BC)))]=0
    
    p_BC_X=colSums(movBC_number*prop_X_BC) # On multiplie movBC_number colonne par colonne par le vecteur prop_X_BC (prop_X_BC ne change pas), puis on fait la somme de chaque colonne
    p_BC_R=colSums(movBC_number*prop_R_BC)

    f_imp_BC_t=names(which((p_BC_X+p_BC_R)!=0))
    
    BC_draw=rmultinomial(n=length(f_imp_BC_t),size=colSums(movBC_number[,f_imp_BC_t]),prob=cbind(p_BC_X[f_imp_BC_t],p_BC_R[f_imp_BC_t]))

    apop$F1[,f_imp_BC_t]=apop$F1[,f_imp_BC_t]+t(BC_draw)
    
    #####
    
    for(i in 1:12){if(any(apop[[i]]<0)){print(paste("Population negative dans le compartiment:",names(apop)[i]))}}
    
  }
  t=t+1
  agep_init=agep_init+1
}

apop_initial=apop

# Farms for which the initial number of sows is not enough for this initial value of "n_imp_gilts"
n_imp_gilts[which((apop_initial$G1["PX",]-n_imp_gilts)<0)]=0

