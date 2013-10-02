################################################
#	### Making Matrix ######################
#	Creation: 2010
#	#Last modified: 02.10.2013
################################################	
################################################
################################################

#English:
#y:data vector. Could also be a matrix or dataframe (one collumn) with all values. Alternatively, a lower triangular matrix. 
#user must specify, thorugh the argument 'data', if the input is a vector ('vector') or a lower triangular matrix. In both cases, the #output will be a square simmetric matrix.


#Important: the way thorugh which the program converts the vector into a matrix (option 'vector') supposes that the vector is the output of paiwise analyses from PAML. In this type of output, the order of the pairs is well established, and, thanks to that, this function is useful.


#Remember CODEML's pairwise (runmode: -2) output looks like this:
#    [,1]
#[1,] v21
#[2,] v31
#[3,] v32
#[4,] v41
#[5,] v42
#[6,] v43

#and what we want, for mantel.easy.r, is:
#     1   2   3    4
#1    0  v21  v31  v41
#2   v21  0   v32  v42
#3   v31 v32   0   v43
#4   v41 v42  v43   0

#Portugues:

#y: vetor com os dados. 
#Pode ser tambem uma matriz ou dataframe de uma coluna com todos os valores. Ou, alternativamente, uma matriz triangular inferior. O usuario deve especificar, atraves do argumento "data" se o input eh um vetor ("vector") ou uma matriz triangular inferior("lt"). Em ambos os casos, o output sera uam mariz quadrada e simetrica

#importante: a forma como o programa converte o vetor em matriz, na opcao "vector" do argumento "data", supoe que o vetor eh um output de resultados de analises par a par do PAML. Nesse tipo de output, a ordem dos "pares" eh bem estabelecida e, gracas a isso, foi possivel escrever a funcao.


#lembrando que o ouput do paml vem assim
#    [,1]
#[1,] v21
#[2,] v31
#[3,] v32
#[4,] v41
#[5,] v42
#[6,] v43

#e o que queremos eh
#     1   2   3    4
#1    0  v21  v31  v41
#2   v21  0   v32  v42
#3   v31 v32   0   v43
#4   v41 v42  v43   0



#######################################################################################################################################

making.matrix<-function(y, data="vector"){
                                        
  

  if(data=="vector"){
    length(y)->obs
                                        #to solve a quadratic equation of type: ax^2 + bx + c=0, Bhaskara:
                                        #x=(-b +/- sqrt(b^2-4ac))/(2a)
                                        #ax^2 + bx +c
                                        #n*(n-1)/2=780
                                        #n^2-n=1560
                                        #n^2-n-1560
                                        #delta=b^2-4*a*c
    delta=(-1)^2-(4*1*(-(obs*2)))
    
    n1=-(-1+sqrt(delta))/(2*1)  #first solution
    n2=-(-1-sqrt(delta))/(2*1) 	#second solution
    
    if(n2>n1){
      n<-n2
    }

    else{
      n<-n1
    }  

    mat<-matrix(nrow=n,ncol=n) #empty matrix
    
    for(i in 1:n){			#diagonals are 0.
      mat[i,i]<-0
    }    	
					
    count<-rep(2:n,c(seq(from=1,to=n-1)))
                                        
    count2<-1
    for(i in 3:n-1){
      count2=c(count2,seq(1:i))
    }
                                      
    counter<-cbind(count, count2)  
                                      
    
    for(i in 1:obs){
      
      mat[counter[i,1], counter[i,2]]<-y[i]->mat[counter[i,2],counter[i,1]] 
	##put the value in its position and in the mirror position in the matrix.
    }
  }
  else{
    if(data!="lt"){   #if input is a lower triangular matrix...
      stop()
    }
    else{
      n <-dim(y)[1]
      obs<-(n*(n-1))/2
      for(i in 1:n){
        y[i,i]<-0
      }                                 #diagonals are 0.
      count<-rep(2:n,c(seq(from=1,to=n-1)))
                                        
      count2<-1
      for(i in 3:n-1){
        count2=c(count2,seq(1:i))
      }
                                        #bind the two  counters.
      counter<-cbind(count, count2)  
                                        
      obs=(n*(n-1))/2
      for(i in 1:obs){
        
        y[counter[i,1], counter[i,2]]->y[counter[i,2],counter[i,1]]#put the value in its position and in the mirror position in the matrix.
      }
      
    }
  }
  if(data=="vector"){
    return(mat)
  }
  else{
    
    return(y)
    
  }
}
#
##########################################
#example:

#as.vector(rnorm(780,0.5, 7))->x1
#x2<-x1 + c(rep(0.05,390), rep(3,390))

#making.matrix(x1)->mat1
#making.matrix(x2)->mat2
