
###################
# Source function #
###################
install.packages("devtools")
library(devtools)
path="~Downloads/eSMC2_2.0.1.tar.gz" # Path to the dowloaded eSMC package
devtools::install_local(path)
#to install without root permissions (in a local folder) use this command:
withr::with_libpaths(new = "~/R/your_local_dir", install_local("eSMC2_2.0.1.tar.gz"))
library(eSMC2) # If library already installed

########
#Script#
########


##############
# Parameters #
##############

M=10 # Haploid sequence number
NC=1 # Number of chrosome/scaffold simulated
nsim=10 # Number of simulation or iteration 
mu=10^-8# Mutation rate
r=5*10^-8 # Recombination rate
n=40 # Number of hidden state
rho=r/mu # ratio of recombination over mutation rate
for(x in 1:nsim){
  
################
# Getting data #
################



    Using_my_script_to_simulate=T # TRUE if using my python script to simulate data
    if(Using_my_script_to_simulate){
    chr=1
    O_total=list()
    setwd(" ") # Location of txt files
    path=paste("Stefan_4",mu,r,"x",(x-1),"chr",(chr-1),"self",self,"sc",scenario,"data",data,".txt",sep ="" ) # File name
    O_total= read_mspts_data(path,M)
    if(length(unique(O_total[dim(O_total)[1],]))<dim(O_total)[2]){
      O_total=O_total[,match(unique(O_total[dim(O_total)[1],]),O_total[dim(O_total)[1],])]
    }
    O_total=O_total[,which(as.numeric(O_total[dim(O_total)[1],])>0)]
    }
    
     
    Using_multihetsepfiles=T  # TRUE if input file is a multihetsep file
    if(Using_multihetsepfiles){
      path=" " # Lo
      filename=paste("  ",sep="")
      O_total=Get_real_data(path,M,filename,"\t")
    }

################
# Getting data #
################
    
    test_=teSMC(n=n,rho=rho,model=list("Free"),O_total,maxit =20,BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),Boxs=list(c(0,0.99)),Constant_Pop=F,estimate="SF",NC=1)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=test_$sigma
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc    
    
    test_=teSMC(n=n,rho=rho,model=list("Constant"),O_total,maxit =20,BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),Boxs=list(c(0,0.99)),Constant_Pop=F,estimate="SF",NC=1)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=rep(test_$sigma,n)
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc    
    
    #########
    # Given #
    #########
    current_selfing_value= # Give current selfing value 
    ancestral_selfing_value=  # Give ancestral selfing value
    given= c(current_selfing_value,ancestral_selfing_value)
    
    test_=teSMC(n=n,rho=rho,model=list("Given Transition"),O_total,maxit =20,BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),Boxs=list(given),Constant_Pop=F,estimate="SF",NC=1)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=test_$sigma
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc    
    
    test_=teSMC(n=n,rho=rho,model=list("One transition"),O_total,maxit =20,BoxB=list(c(0.05,1)),BoxP=c(3,3),Boxr=list(c(1,1)),Boxs=list(list(c(0.5,0.99),c(0,0.5))),Constant_Pop=F,estimate="SF",NC=1)
    total[[(1+length(total))]]=as.numeric(test_$Xi)
    s_ex[[(1+length(s_ex))]]=test_$sigma
    mu_ex[[(1+length(mu_ex))]]=mean(test_$mu)
    Tc[[(1+length(Tc))]]=test_$Tc    
    
}


########
# Plot #
########
if(T){
  gen <- 1
  setwd(" ") # Plot location
  name_sc=c("Constant","Bottleneck"," Expansion"," Decrease")
  col_u=c("red","orange","green","blue","purple")

  #eSMC
  gen <- 1 # number of year of one generation

  mat_save_p=matrix(NA,nrow=(4*4*nsim),ncol=40)
  mat_save_t=matrix(NA,nrow=(4*4*nsim),ncol=40)
  mat_save_s=matrix(NA,nrow=(4*4*nsim),ncol=40)
  count_p=0
  count_t=0
  count_s=0
  
  
  #eSMC
  #
  

  pdf_name=paste("Figure_Stefan_",sc,"_.pdf",sep = "")# Name of PDF
    pdf(pdf_name) 

        plot(c(10,10^6),c(1,1), log=c("x"), ylim =c(3,5) ,
             type="n", xlab= paste("Generations ago",sep=" "), ylab="population size (log10)",main = "")

          for(x in 1:nsim){
            
            test_<-total[[1+((x-1)*4)]]
            Pop_=mu_ex[[1+((x-1)*4)]]/(mu)
            lines((Tc[[1+((x-1)*4)]]*Pop_), log10((test_)*0.5*Pop_), type="s", col=col_u[1])
            count_p=count_p+1
            count_t=count_t+1
            mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
            mat_save_t[count_t,]=(Tc[[1+((x-1)*4)]]*Pop_)
            
            test_<-total[[2+((x-1)*4)]]
            Pop_=mu_ex[[2+((x-1)*4)]]/(mu)
            lines((Tc[[2+((x-1)*4)]]*Pop_), log10((test_)*0.5*Pop_), type="s", col=col_u[2])
            count_p=count_p+1
            count_t=count_t+1
            mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
            mat_save_t[count_t,]=(Tc[[2+((x-1)*4)]]*Pop_)
            
            test_<-total[[3+((x-1)*4)]]
            Pop_=mu_ex[[3+((x-1)*4)]]/(mu)
            lines((Tc[[3+((x-1)*4)]]*Pop_), log10((test_)*0.5*Pop_), type="s", col=col_u[3])
            count_p=count_p+1
            count_t=count_t+1
            mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
            mat_save_t[count_t,]=(Tc[[3+((x-1)*4)]]*Pop_)
            
            
            test_<-total[[(4*x)]]
            Pop_=mu_ex[[(4*x)]]/(mu)
            lines((Tc[[(4*x)]]*Pop_), log10((test_)*0.5*Pop_), type="s", col=col_u[4])
            count_p=count_p+1
            count_t=count_t+1
            mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
            mat_save_t[count_t,]=(Tc[[(4*x)]]*Pop_)
            
            
          }
        

        legend("topright",legend=c("Free","Constant","Given Transition","One transition"), col=col_u[1:4], lty=c(1,1),cex=0.75,x.intersp=0.5,y.intersp=0.8)
      
    
    dev.off()
  
  
  #eSMC
  # 
  

    pdf_namepaste("Figure_Stefan_selfing_",sc,"_.pdf",sep = "")# Name of PDF
    pdf(pdf_name) 

      plot(c(10,10^6),c(1,1), log=c("x"), ylim =c(0,1.1) ,
           type="n", xlab= paste("Generations ago",sep=" "), ylab="selfing rate",main = "")
          for(x in 1:nsim){
            Pop_=mu_ex[[1+((x-1)*4)]]/(mu)
            lines((Tc[[1+((x-1)*4)]]*Pop_),s_ex[[1+((x-1)*4)]], type="s", col=col_u[1])
            count_s=count_s+1
            mat_save_s[count_s,]=s_ex[[1+((x-1)*4)]]
            
            Pop_=mu_ex[[2+((x-1)*4)]]/(mu)
            lines((Tc[[2+((x-1)*4)]]*Pop_),s_ex[[2+((x-1)*4)]], type="s", col=col_u[2])
            count_s=count_s+1
            mat_save_s[count_s,]=s_ex[[2+((x-1)*4)]]
            
            Pop_=mu_ex[[3+((x-1)*4)]]/(mu)
            lines((Tc[[3+((x-1)*4)]]*Pop_),s_ex[[3+((x-1)*4)]], type="s", col=col_u[3])
            count_s=count_s+1
            mat_save_s[count_s,]=s_ex[[3+((x-1)*4)]]
            
            Pop_=mu_ex[[(4*x)]]/(mu)
            lines((Tc[[(4*x)]]*Pop_),s_ex[[(4*x)]], type="s", col=col_u[4])
            count_s=count_s+1
            mat_save_s[count_s,]=s_ex[[(4*x)]]
            
          }
      legend("topright",legend=c("Free","Constant","Given Transition","One transition"), col=col_u[1:4], lty=c(1,1),cex=0.75,x.intersp=0.5,y.intersp=0.8)

      
    
    
    
    
    dev.off()
    
  
  name_csv=paste("Stuett",collapse ="") # CSV name
  write.csv(mat_save_p,file = paste("Figure_Stefan_mat_save_p_real","sc",,".csv",sep="_")) # Saves population size
  write.csv(mat_save_t,file = paste("Figure_Stefan_mat_save_t_real","sc",name_csv,".csv",sep="_")) # Saves time (x axis)
  write.csv(mat_save_s,file = paste("Figure_Stefan_mat_save_s_real","sc",name_csv,".csv",sep="_")) # Saves selfing value
}
