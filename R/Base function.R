
CoDa <- function(dt, sex = c("M", "F", "T"), hstr=c(1,2), ih, years.fit){
  source("CoDa functions.R")
  
  # hstr=2
  # sex = "M"
  
  ages <- as.numeric(colnames(dt[,-c(1:3)]))
  m <- length(ages)
  years <- as.numeric(unique(dt$Year))
  n <- length(years)
  CoD <- unique(dt$Cause_Rev)
  k <- length(CoD)
  
  nomiCoD <- c("Infectious diseases",
               "Cancer Smoking",
               "Cancer No Smoking",
               "Diseases of the circulatory system",
               "Diseases of the respiratory system",
               "External diseases",
               "Diseases of the digestive system",
               "Other diseases")
  
  dati0 <- as.matrix(subset(dt, Sex== sex)[,1:m+3])
  
  dati1 <- array(c(t(dati0)), dim=c(m,n,k))
  dimnames(dati1) <- list(ages, years, nomiCoD)
  
  dati <- dati1[7:18,,]
  
  if(sex == "M"){LT.1 <-(read.table("LT/mltper_5x1.txt",skip=1,header=TRUE,stringsAsFactors=F))}
  if(sex == "F"){LT.1 <-(read.table("LT/fltper_5x1.txt",skip=1,header=TRUE,stringsAsFactors=F))}
  if(sex == "T"){LT.1 <-(read.table("LT/bltper_5x1.txt",skip=1,header=TRUE,stringsAsFactors=F))}
  
  dx.nat <- matrix(LT.1$dx,length(unique(LT.1$Age)),length(unique(LT.1$Year)))
  colnames(dx.nat) <- unique(LT.1$Year)
  rownames(dx.nat) <- unique(LT.1$Age)
  
  dx.nat1 <- dx.nat[,colnames(dx.nat) %in% years]
  
  rownames(dx.nat1) <- unique(LT.1$Age)
  
  age.temp<- unique(LT.1$Age)
  
  dx.nat2 <- dx.nat1
  
  dx.nat3 <- rbind(
    dx.nat2[7:17, , drop = FALSE],
    "80+" = colSums(dx.nat2[18:24, , drop = FALSE])
  )
  
  
  dati.total <- (dati[,,1]+dati[,,2]+dati[,,3]+dati[,,4]+
                    dati[,,5]+dati[,,6]+dati[,,7]+dati[,,8])
  
  ages2 <- ages[7:18]
  
  dati.rela1 <- array(NA,dim=c(ncol(dati[,,1]),nrow(dati[,,1]),8))
  dimnames(dati.rela1) <- list(years, ages2, nomiCoD)

  m.cause.i <- dati.temp <- dati *0 
  
  # calculate life table deaths 
  for(i in 1:8){
    dati.temp[,,i] <- ifelse(dati[,,i]==0,0.25,dati[,,i])

    dati.rela1[,,i] <- t(dx.nat3*(dati.temp[,,i]/dati.total))
  }
  
  #RISCALARE AFFINCHè SOMMINO A 10000
  for (j in 1:n) {
    somma <- sum(dati.rela1[j,,])
    fattore <- 100000 / somma
    dati.rela1[j,,] <- dati.rela1[j,,] * fattore
  }
  
  # life table deaths 
  if(hstr==1){
    dx.com <- cbind(dati.rela1[,,1],dati.rela1[,,2],dati.rela1[,,3],
                   dati.rela1[,,4],dati.rela1[,,5],dati.rela1[,,6],
                   dati.rela1[,,7],dati.rela1[,,8])}
  
  if(hstr==2){
    dx.com <- cbind(dati.rela1[,,1] + dati.rela1[,,2] + dati.rela1[,,3] +
                    dati.rela1[,,4] + dati.rela1[,,5] + dati.rela1[,,6] +
                    dati.rela1[,,7] + dati.rela1[,,8])}
  
  # Fit the models 
  years.fitfor <- c(years.fit,(max(years.fit)+1):(max(years.fit)+ih))
  years.for <- c((max(years.fit)+1):(max(years.fit)+ih))
  
  if(hstr == 2){
    nomiCoD = "All-causes"
  }
  
  model.fit.for.CT <- CoDa.CT(dx=dx.com[(1:length(years.fit)),],
                              ih=ih, 
                              k=length(nomiCoD),
                              years=years.fit,
                              ages=ages2,
                              ses=nomiCoD)
  
  fitfor <- as.data.frame(model.fit.for.CT$dx.forcast) %>% 
    mutate(Year = row.names(.), .before=1) %>% 
    pivot_longer(cols = c(2:ncol(.)), names_to = "AgeCause", values_to = "dx") %>% 
    separate(col = "AgeCause", into = c("Age", "Cause"), sep = "\\.") %>% 
    mutate("Type" = case_when(
      Year %in% years.fit ~ "Fitted",
      Year %in% years.for ~ "Forecasted",))
  
  if(hstr == 2){
    dati.rela1 <- dati.rela1[,,1]+dati.rela1[,,2]+dati.rela1[,,3]+
      dati.rela1[,,4]+dati.rela1[,,5]+dati.rela1[,,6]+
      dati.rela1[,,7]+dati.rela1[,,8]
    
    colnames(dati.rela1) <- paste(colnames(dati.rela1), "All-causes", sep=".")
  }
  
  inoutsample <- as.data.frame(dati.rela1) %>% 
    mutate(Year = row.names(.), .before=1) %>% 
    pivot_longer(cols = c(2:ncol(.)), names_to = "AgeCause", values_to = "dx") %>% 
    separate(col = "AgeCause", into = c("Age", "Cause"), sep = "\\.") %>% 
    mutate("Type" = case_when(
      Year %in% years.fit ~ "Observed_In_sample",
      Year %in% years.for ~ "Observed_Out_of_sample",))
  
  
  total <- rbind(inoutsample, fitfor)  %>% 
    group_by(Year, Type, Age) %>% 
    summarise(dx = sum(dx)) %>% 
    mutate(cumsum_dx = cumsum(dx))%>% 
    group_by(Year, Type) %>%
    mutate(lx = case_when(
      Age != 25 ~ 100000- lag(cumsum_dx),
      .default = 100000)) %>% 
    dplyr::select(-cumsum_dx) %>% 
    mutate(Lx = 2.5*(lx+lead(lx))) %>% 
    dplyr::select(-dx)
  
  CoDa_Output <- rbind(inoutsample, fitfor) %>% 
    left_join(total, by=c("Year", "Age", "Type")) %>% 
    mutate(mx = ifelse(Age!=80, dx/Lx, (2/5)*dx/lx), .before=Type) %>% 
    dplyr::select(-c(lx, Lx)) %>% 
    mutate(Cause_cod = as.numeric(as.factor(Cause)), .after=Cause,
           Type = factor(Type, levels =  c("Observed_In_sample",
                                           "Fitted",
                                           "Observed_Out_of_sample",
                                           "Forecasted"))) %>% 
    arrange(Type, Year, Cause_cod, Age)
  
  if(hstr == 2){
    CoDa_Output <- CoDa_Output %>% 
      mutate(Cause_cod = 999)
  }
    
  return(CoDa_Output)

}


