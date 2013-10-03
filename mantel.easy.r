#############################################
################Mantel Easy##################
#############################################
#############################################
####Creation: 2010 ##########################
####Last modified: 02.10.2013################
#############################################


#mantel.easy                         package:NA                    R Documentation


#Mantel Test


#Description:

#A simple function which executes a Mantel Test. The input must consist of two distance/similarity square matrices. 
#The function allows visualization of a histogram os the correlations (null distribution) and the choice of the preferred method for #correlation calculation. Also, and most importantly, it provides alternatives for the treatment of missing data. 



#Uma função simples que executa um teste de Mantel baseado em permutações. Os dados de entrada devem ser duas matrizes de distância/#similaridade simétricas e quadradas. A função permite a visualização do histograma de correlações (distribuição nula) e a escolha do #método de cálculo de correlação, além de fornecer alternativas para o tratamento de dados faltantes.


#Usage:

#mantel.easy(x,y,nperm=1000, method="pearson", na.action="complete.obs",hist=T)


#Arguments:

#x,y: distance similarity matrix 1 and distance/similarity matrix 2.

#nperm: number of permutations to be performed. These will generate the null distribution for the test.

#method: correlation calculation method ('pearson', 'kendall', 'spearman').

#na.action: how to handle missing data (NAs).

#hist: logical. If 'TRUE', the function will plot an histogram with the correlations (null distribution) and show where the real correlation fits into that distribution.


################################################################################################
#the code.

mantel.easy<-function(x,y,nperm=1000, method="pearson",na.action="complete.obs",hist=T){
  ##first: check if both objects are square matrices.
  if(ncol(x)!=nrow(x)|ncol(y)!=nrow(y)){
    cat("ERROR! You must provide two distance/similarity SQUARE matrices as input\n")
    stop()
    break
  }
  ##check if the two matrices have the same dimensions.
  if(dim(x)[1]!=dim(y)[1]){
    cat("ERROR! Distance matrices must have the same length\n")
    stop()
    break
  }
  
  n<-dim(x)[1]#number of lines in each matrix (same as number of collumns)
  
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
    na.m2<-sum(is.na(m2)) #
    if(na.action=="overall"){
      mean(m1, na.rm=TRUE) -> mean         
      m1[is.na(m1)] <- mean #
      mean(m2, na.rm=TRUE) -> mean2
      m2[is.na(m2)] <- mean2 
      mean(m12, na.rm=TRUE) -> mean3
      m12[is.na(m12)] <- mean3         
      cor(m1,m2, method=method, use="all.obs") ->cor
      resmat<-rep(NA, nperm) 
      for (i in 1:nperm) {   
        samp<-sample(1:n) #
        a<-m12[samp,samp] 
        b<-m2 #
        a<-a[row(a)!=col(a)] #
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




##################################################
#example

#as.vector(rnorm(780,0.5, 7))->x1
#x2<-x1 + c(rep(0.05,390), rep(3,390))

#making.matrix(x1)->mat1
#making.matrix(x2)->mat2


#mantel.easy(mat1, mat2)
