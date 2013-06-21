mantel.easy<-function(x,y,nperm=1000, method="pearson",na.action="complete.obs",hist=T){
  ##primeiro: verificar se ambos os objetos s„o matrizes quadradas
  if(ncol(x)!=nrow(x)|ncol(y)!=nrow(y)){
    cat("ERROR! You must provide two distance/similarity SQUARE matrices as input\n")
    stop()
    break
  }
  ##verificar se as matrizes tem mesmo tamanho
  if(dim(x)[1]!=dim(y)[1]){
    cat("ERROR! Distance matrices must have the same length\n")
    stop()
    break
  }
  
  n<-dim(x)[1]#n˙mero de colunas ou linhas de cada matriz
  
  if(nperm != 0){
    j<-seq(from=1, to=nperm,by=nperm/10)
    for(i in j){
      cat (i, "of", nperm, "\n")
    }
    cat(nperm, "of", nperm,"\n")
    m12<-x #coloca uma das matrizes(a primeira) no objeto m12
    m1<-x[row(x)!=col(x)]
    m2<-y[row(y)!=col(y)]
    na.m1<-sum(is.na(m1))
    na.m2<-sum(is.na(m2)) #coloca em m1 e m2 apenas os valores n„o diagonais de x e y (que ser„o sempre zero, por definiÁ„o)
    if(na.action=="overall"){
      mean(m1, na.rm=TRUE) -> mean         # mÈdia dos valores n„o-NA
      m1[is.na(m1)] <- mean #imputa media no lugar dos NAs
      mean(m2, na.rm=TRUE) -> mean2
      m2[is.na(m2)] <- mean2 
      mean(m12, na.rm=TRUE) -> mean3
      m12[is.na(m12)] <- mean3         
      cor(m1,m2, method=method, use="all.obs") ->cor
      resmat<-rep(NA, nperm) 
      for (i in 1:nperm) {   
        samp<-sample(1:n) #pra cada uma das 1000 perm, amostra de n (n˙mero de linhas da matriz quadrada) valores entre 1 e n
        a<-m12[samp,samp] 
        b<-m2 #a outra matriz vai ser a mesma definida anteriormente, porque so uma matriz eh permutada nesse teste
        a<-a[row(a)!=col(a)] #aqui elimina-se da matriz permutada os valores diagonais (ou seja, zero)
        resmat[i]=cor(a,b,use="all.obs", method=method) 
      }
    }
    else{
      cor(m1,m2, method=method, use=na.action) ->cor
      resmat<-rep(NA, nperm)
      for (i in 1:nperm) {
        samp<-sample(1:n) 
        a<-m12[samp,samp] 
        b<-m2 
        a<-a[row(a)!=col(a)]
        resmat[i]=cor(a,b,use=na.action, method=method) 
      }
    }
    if(cor>0){
      
      p<-sum(cor>resmat)/(nperm) 
      p <- min(c(p, 1 - p))
      #/(nperm)
    }
    else{
      p<-sum(cor<resmat)/(nperm)
      p <- min(c(p, 1 - p))
      #/(nperm)
    }
  }
  else{
    p<-NA
    cat ("I cannot perform a mantel test without a null distribution! 'nperm' should be a positive value!\n")
    stop()
    break()
  }
   
  if(hist==T){
    
    if(cor>max(resmat)){
      hist(resmat, xlim=c(min(resmat),cor),nclass=nperm/10,border="cornflowerblue", main="", xlab="Correlation", ylab="Frequency",col="cornflowerblue")
      abline(v=cor, col="orange", lty=2)
      
    }
    if(cor<min(resmat)){  
      hist(resmat, xlim=c(cor,max(resmat)),nclass=nperm/10,border="cornflowerblue", main="",xlab="Correlation", ylab="Frequency",col="cornflowerblue")
      abline(v=cor, col="orange", lty=2)
      
    }
    if(cor>=min(resmat) & cor<=max(resmat)){  
      hist(resmat, xlim=c(min(resmat),max(resmat)),nclass=nperm/10,border="cornflowerblue", main="", xlab="Correlation", ylab="Frequency",col="cornflowerblue")
      abline(v=cor, col="orange", lty=2)
      
    }
    
  }
  if(na.action=="overall"){
    cat("Overall mean imputation over",na.m1, "values from matrix x and", na.m2, "values from matrix y\n")
    final<-list(correlation=cor,p.value=p,NAs=c(na.m1,na.m2,"Overall Mean Imputation"))
  }
  if(na.action=="complete.obs"){
    
    cat("Casewise deletion of", na.m1+na.m2, "pairs of values from the two matrices\n")
    final<-list(correlation=cor,p.value=p,NAs=c(na.m1, na.m2, "Casewise Deletion"))
  }
  if(na.action!="overall" & na.action!="complete.obs"){
    cat("Nothing was done with the",na.m1+na.m2,"missing values\n")
    
    final<-list(correlation=cor,p.value=p, NAs=c(na.m1, na.m2, na.action))
  }
  class(final)="Mantel"
  return(final)
}


