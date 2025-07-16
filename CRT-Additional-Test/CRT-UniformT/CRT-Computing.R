#-----------------------------------------------------------
# Function: CSSCRT
# Purpose:  Conducts a conditional randomization test based on exchangeable sampling
# Arguments:
#   Y        - Response vector
#   X        - Data matrix
#   graph    - Adjacency matrix of the graph
#   setT     - indices of variables to be tested for conditional independence
#   CSS_samples - Optional, precomputed exchangeable samples (of setT) and x_hub
#   Fun_Stat - Function (X,Y, setT,...) to compute test statistic
#   Stat_Procedure - Choose among 1) directly apply Fun_Stat, 
#                                 2) Regress Y on X[-T] and then on X[T]
#                                 3) Regress Y on X[T] and X[-T] with feature selection
#   Fun_FitNull -  Regress Y on X[-T] and return fitted.value and imp.var
#   ...       - Arguments to pass to Fun_Stat and Fun_FitNull
# Returns:   Computed p-value
#-----------------------------------------------------------
CSSCRT <- function(Y, X, setT, CSS_samples = NULL, graph=NULL,
                   Fun_Stat, 
                   Stat_Procedure=c("Direct","Distill", "SelDist",
                                    "Residual","Select"),
                   Fun_FitNull=NULL, k_sel=NULL,
                   type='One-sided',
                   randomize=T,M=100,L=1, ...) {
  
  # Generate exchangeable samples if not provided
  if (is.null(CSS_samples)) {
    CSS_samples <- exchangeable_sampling(X, graph,I = setT, M=M, L=L,bReduce=T)
  }else{
    M=length(CSS_samples$X_copies)
  }
  
  X_copies <- CSS_samples$X_copies
  X_hub <- CSS_samples$X_hub
  
  p=ncol(X)
  size_T = length(setT)
  new_setT=1:size_T
  
  if(Stat_Procedure=="Direct"){
    T0 <- Fun_Stat(X,Y,setT,...)
    X_tilde=X
    Ttilde=list()
    for(j in 1:M){
      X_tilde[,setT]=X_copies[[j]];
      Ttilde[[j]]=Fun_Stat(X_tilde,Y,setT,...)
    }
  }else if(Stat_Procedure=="Residual"){
    Model_FitNull=Fun_FitNull(X[,-setT,drop=F],Y, k=-1,...)
    Residual_Null= as.numeric(Y) - as.numeric(Model_FitNull$fitted.value)
    
    T0 <- Fun_Stat(X[,setT,drop=F], Residual_Null,...)
    
    
    Ttilde <- lapply(1:M, function(j){
      Fun_Stat(X_copies[[j]],Residual_Null,...)
      })
    
  }else if(Stat_Procedure=="Select"){
    ## Specify the number of importance
    if (is.null(k_sel)){
      
      k_sel = min(as.integer(2 * log(p)), p - size_T ) 
    }else{
      k_sel = 0  # Let the function itself choose
    }
    
    setV=setdiff(1:p , setT)
    Model_FitNull=Fun_FitNull(X[,setV,drop=F],Y, k_sel=k_sel,...)
    setS=setV[Model_FitNull$imp.var]
    k_sel=length(setS)
    
    
    T0 <- Fun_Stat(X,Y, setT,sel_var=setS, ...)
    
    X_tilde=X
    Ttilde=list()
    for(j in 1:M){
      X_tilde[,setT]=X_copies[[j]];
      
      Ttilde[[j]]=Fun_Stat(X_tilde,Y, setT,sel_var=setS,...)
    }
    
  }else if(Stat_Procedure=="Distill"){
    
    Model_FitNull=Fun_FitNull(X[,-setT,drop=F],Y, k=-1,...)
    
    
    T0 <- Fun_Stat( cbind( X[,setT,drop=F], Model_FitNull$fitted.value),   
                   setT=new_setT, Y,...)
    
    
    Ttilde <- lapply(1:M, function(j){
      Fun_Stat(cbind(X_copies[[j]], Model_FitNull$fitted.value),   
               setT=new_setT, Y,...)
      
    })
    
    
  }
  
  
  return(ComputePValue(T0, Ttilde, type=type, randomize = randomize))
}
