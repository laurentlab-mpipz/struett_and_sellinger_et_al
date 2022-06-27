###################
# Source function #
###################
#library(utils)
library(BB)
library(parallel)
library(readr)
################################

###################
# Source function #
###################

################
# SNP finders ##
################

# Build signal for the model
# DNA : matrix of segregating sites for the whole sample
# n : sequence length
# output : sequence of  0 and 1, where 0 means that nucleotides are the same and 1 that they are different between the two selected sequences
seqSNP2_MH<-function(DNA,n){
  DNA=as.matrix(DNA)
  pos=which(DNA[1,]!=DNA[2,])
  seq=rep(0,n)
  seq[as.numeric(DNA[3,pos])]=1
  return(seq)
}



##############
#File reader #
##############

# Transforms a snp file  in a segregating sites matrix
# DNAfile : the snp file
# L : sequence length
# M: number of individual
# s : number of repetition/Different analysis
Seqlist2data<-function(DNAfile,L,M,s){
  Output=list()
  vect_t=vector()
  s_count=0
  for(i in 1:length(DNAfile)){
    DNA=matrix(0,nrow=M+1)
    if(!is.na(DNAfile[[i]][1])){
      if(nchar(DNAfile[[i]][1])>2){
        if(substr(DNAfile[[i]][1],1,3)=="seg"){
          if(as.numeric(DNAfile[[i]][2])==0){
            s_count=s_count+1
            Output[[s_count]]=DNA
          }
          if(as.numeric(DNAfile[[i]][2])==1){
            s_count=s_count+1
            position=round(L*as.numeric(DNAfile[[i+1]][2]))
            DNA=matrix(0,nrow=M+1,ncol=1)
            for(k in 1:M){
              if(DNAfile[[i+1+k]][1]=="1"){
                DNA[k,1]=1
              }
            }
            DNA[(M+1),1]=position
            Output[[s_count]]=DNA
          }
          if(as.numeric(DNAfile[[i]][2])>1){
            s_count=s_count+1
            if(substr(DNAfile[[i+1]][1],1,3)=="pos"){
              position=as.numeric(DNAfile[[i+1]][2:length(DNAfile[[i+1]])])
              r_position=round(L* position)
              r_position_u=as.numeric(as.matrix(unique(r_position)))
            }
            if(length(r_position)>length(unique(r_position))){
              print("Length problem!Mutation appeared more than once on same site!")
              print("Extra mutation will be removed!")
            }
            DNA=matrix(0,nrow=M+1,ncol=length(r_position_u))
            DNA[(M+1),]=r_position_u
            SEQ_data=vector()
            for(k in 1:M){
              SEQ_data=rbind(SEQ_data,as.numeric(as.vector(unlist(strsplit(DNAfile[[i+1+k]][1],"")))))
            }
            print(dim(SEQ_data))
            print(length(r_position_u))
            xy=match(r_position_u,r_position)
            DNA[1:M,]=SEQ_data[,xy]
            Output[[s_count]]=DNA
          }
        }
      }
    }
  }
  if(s_count==s){
    print("All sample were analysed")
  }
  
  return(Output)
  
  
}

# Find simulation output and transform it into a SNP file
# path : location of the simulation output and its name
old_Get_data<-function(path){
  DNAseqfile=list()
  count_DNA=0
  con = file(path, "r")
  while ( TRUE ) {
    count_DNA= count_DNA+1
    line = readLines(con, n = 1)
    if ( length(line) == 0 ){
      break
    }
    DNAseqfile[[count_DNA]]=as.vector(unlist(strsplit(as.vector(line)," ")))
  }
  close(con)
  return(DNAseqfile)
}

# Find simulation output and transform it into a SNP file
# path : location of the simulation output and its name
Get_data<-function(path,heavy=F){
  DNAseqfile=list()
  count_DNA=0
  con = file(path, "r")
  while ( TRUE ) {
    count_DNA= count_DNA+1
    line = readLines(con, n = 1)
    if ( length(line) == 0  ){
      break
    }
    DNAseqfile[[count_DNA]]=as.vector(unlist(strsplit(as.vector(line)," ")))
    if (heavy&substr(line,1,3)=="seg"){
      break
    }
  }
  close(con)
  return(DNAseqfile)
}

Get_vcf_data<-function(path){
  DNAseqfile=list()
  count_DNA=0
  con = file(path, "r")
  while ( TRUE ) {
    count_DNA= count_DNA+1
    line = readLines(con, n = 1)
    if ( length(line) == 0 ){
      break
    }
    DNAseqfile[[count_DNA]]=as.vector(unlist(strsplit(as.vector(line),"\t")))
  }
  close(con)
  return(DNAseqfile)
}

Process_vcf_data<-function(path){
  data=Get_vcf_data(path)
  M=length(data[[15]][c(10:length(data[[15]]))])*2
  SNP_mat=matrix(0,ncol = length(data),nrow = (M+1))
  i=0
  for(ii in 7:length(data)){
    i=i+1
    SNP_mat[(M+1),i]=as.numeric(data[[ii]][2])
    temp_seq=strsplit(paste(data[[ii]][c(10:length(data[[ii]]))],collapse = ""),"")[[1]]
    SNP_mat[1:M,i]=temp_seq[which(temp_seq!="|")]
  }
  return(SNP_mat)
}


# creates a snp file from a msmc input file
# path : location of the file
# M : number of  haplotypes
create_good_file<-function(path=NA,M,filename){
  old_wd=getwd()
  if(!is.na(path)){
    setwd(path)
  }
  con = file(filename, "r")
  a=system(paste("wc -l",filename,sep=" "),intern =T)
  nb_line=as.numeric(substr(a,1,(nchar(a)-1-nchar(filename))))
  pos=0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    pos=pos+1
    if(pos==1){
      pattern="\t"
      nb_col=length(as.vector(unlist(strsplit(as.vector(line),pattern))))
      if(nb_col==1){
        nb_col=length(as.vector(unlist(strsplit(as.vector(line)," "))))
        if(nb_col>1){
          pattern=" "
        }
      }
      DNAseqfile=matrix(NA,nrow=nb_line,ncol =nb_col )
    }
    DNAseqfile[pos,]=as.vector(unlist(strsplit(as.vector(line),pattern)))
  }
  close(con)
  L=as.numeric(DNAseqfile[dim(DNAseqfile)[1],2])
  DNAseqfile=DNAseqfile[,-1]
  output=matrix(NA,ncol=nb_line,nrow=(1+M))
  for(i in 1:dim(DNAseqfile)[1]){
    output[1:M,i]=as.vector(unlist(strsplit(DNAseqfile[i,3],"")))[1:M]
    output[(M+1),i]=as.numeric(DNAseqfile[i,1])
  }
  res=list()
  res$L=L
  res$output=output
  res$DNAseqfile=DNAseqfile
  setwd(old_wd)
  return(res)
}

create_good_file_new<-function(path,M,filename){
  old_wd=getwd()
  if(!is.na(path)){
    setwd(path)
  }
  DNAseqfile=readr::read_delim(filename," ",col_names = F)
  DNAseqfile=as.matrix(DNAseqfile)
  L=as.numeric(DNAseqfile[dim(DNAseqfile)[1],2])
  DNAseqfile=DNAseqfile[,-1]
  output=matrix(NA,ncol=dim(DNAseqfile)[1],nrow=(2+M))
  for(i in 1:dim(DNAseqfile)[1]){
    output[c(1:M),i]=as.vector(unlist(strsplit(DNAseqfile[i,3],"")))[1:M]
    output[(M+1),i]=as.numeric(DNAseqfile[i,2])
    output[(M+2),i]=as.numeric(DNAseqfile[i,1])
  }
  res=list()
  res$L=L
  res$output=output
  setwd(old_wd)
  return(res)
}


get_rec_smcpp<-function(file,path=NA){
  if(!is.na(path)){
    setwd(path)
  }
  test=read.delim(file, skipNul = T)
  test=as.vector(as.matrix(test))
  rho_vect=substr(test,5,7)
  rm=as.numeric(substr(test[which(rho_vect=="rho")],10,(nchar(test[which(rho_vect=="rho")])-1)))/as.numeric(substr(test[(1+which(rho_vect=="rho"))],12,nchar(test[(1+which(rho_vect=="rho"))])))
  return(rm)
}

#########################
# MSMC related function #
#########################

# Create a MSMC input file  from a segregating site matrix
create_realinput_msmc2<-function(O,name,num=F){
  options(scipen=999)
  O=as.matrix(O)
  M=dim(O)[1]-1
  n=dim(O)[2]
  vect_opti=as.numeric(O[M+1,])
  O=O[-(M+1),]
  if(num){
    for(seq in 1:dim(O)[1]){
      O1=as.numeric(O[seq,])
      Apos=which(O1==1)
      Tpos=which(O1==0)
      O1[Apos]="A"
      O1[Tpos]="T"
      O[seq,]=O1
    }
    
  }
  output=matrix(NA,nrow=length(vect_opti),ncol = 4)
  count=0
  for(p in 1:length(vect_opti)){
    diff=length(unique(O[,p]))
    if(diff>1){
      count=count+1
      output[count,1]=1
      output[count,2]=as.numeric(format(as.numeric(vect_opti[p]),scientific = F))
      if(count==1){
        output[count,3]=as.numeric((as.numeric(format(as.numeric(vect_opti[p]),scientific = F))-1),scientific = F)
      }
      if(count>1){
        output[count,3]=as.numeric(format(as.numeric(vect_opti[p])-as.numeric(output[(count-1),2]),scientific = F))
      }
      output[count,4]=paste(O[,p],collapse="")
    }
  }
  output=output[1:count,]
  if(length(vect_opti)>1){
    write.table(output, file=paste(name,".txt",sep=""), quote = FALSE, row.names=FALSE, col.names=FALSE)
  }
  if(length(vect_opti)==1){
    write.table(t(output), file=paste(name,".txt",sep=""), quote = FALSE, row.names=FALSE, col.names=FALSE)
  }
  options(scipen=0)
}

# Create a MSMC input file  from a segregating site matrix
create_realinput_msmc2_masked<-function(O,name,num=F,mask){
  options(scipen=999)
  O=as.matrix(O)
  M=dim(O)[1]-1
  n=dim(O)[2]
  vect_opti=as.numeric(O[M+1,])
  O=O[-(M+1),]
  if(num){
    for(seq in 1:dim(O)[1]){
      O1=as.numeric(O[seq,])
      Apos=which(O1==1)
      Tpos=which(O1==0)
      O1[Apos]="A"
      O1[Tpos]="T"
      O[seq,]=O1
    }
    
  }
  output=matrix(NA,nrow=length(vect_opti),ncol = 4)
  count=0
  for(p in 1:length(vect_opti)){
    diff=length(unique(O[,p]))
    if(diff>1){
      count=count+1
      output[count,1]=1
      output[count,2]=as.numeric(format(as.numeric(vect_opti[p]),scientific = F))
      if(count==1){
        output[count,3]=as.numeric((as.numeric(format(as.numeric(vect_opti[p]),scientific = F))-1),scientific = F)
      }
      if(count>1){
        output[count,3]=as.numeric(format(as.numeric(vect_opti[p])-as.numeric(output[(count-1),2]),scientific = F))
      }
      output[count,4]=paste(O[,p],collapse="")
    }
  }
  output=output[1:count,]
  
  for(mm in 1:length(mask)){
    if(mask[[mm]][1]<as.numeric(output[1,2])){
      output[(1),3]= as.numeric(output[1,3])-(mask[[mm]][2]-mask[[mm]][1])
    }else{
      pos_seg=max(which(as.numeric(output[,2])<mask[[mm]][1]))
      if(length(pos_seg)>0){
        if((pos_seg+1)<=dim(output)[1]){
          output[(pos_seg+1),3]= as.numeric(output[(pos_seg+1),3])-(mask[[mm]][2]-mask[[mm]][1])
        }
      }
    }
  }
  if(any(as.numeric(output[,3])<0)){
    stop("Problem in mask file")
  }
  if(length(vect_opti)>1){
    write.table(output, file=paste(name,".txt",sep=""), quote = FALSE, row.names=FALSE, col.names=FALSE)
  }
  if(length(vect_opti)==1){
    write.table(t(output), file=paste(name,".txt",sep=""), quote = FALSE, row.names=FALSE, col.names=FALSE)
  }
  options(scipen=0)
}


####################
# TIB PSMC FUNCTION#
####################

# Runs eSMC
# maxit : maximum number of iteration for the baum_welch algorithm
# symbol_number : number of symbol use for the  zipping  process ( must be  less than 30)
# BoxB : boundaries of the germination rate ( first value must be  bigger than 0)
# BoxP : logarithmic boundaries for the  demography e.g. c(3,3) means the  population size can grow up to a thousand time  and decrease up to a thousand time
# Boxr : logarithmic boundaries for the  recombination rate e.g. c(1,) means the recombination rate can be up to  ten times smaller of bigger than the initial given value
# Boxs : boundaries for the self-fertilization rate e.g. c(0.5,0.9) means the selfing rate is between 0.5 and 0.9
# MB : number of iteration for the model ( set to  2 if  SF or SB is TRUE set 1 otherwise)
# NC : number of different chromosome analysed
# pop_vect : vector of grouped demography estimator e.g. c(1,1,1,1,1,1) in every hidden state the population size will be estimated and c(2,2,2) population size is the same in the two hidden states ( sum of the vector should always be equal to the number of hidden states)
teSMC<-function(n=40,rho=1,O,mu_r,maxit =20,symbol_number=30,BoxB=c(0.05,1),BoxP=c(2,2),Boxr=c(1,1),Boxs=c(0,0.99),N_max=7,pop=F,SB=F,SF=F,Rho=T,Check=F,BW=F,FS=F,NC=1,mu_b=1,pop_vect=NA,rec_vect=NA,sigma=0,beta=1,SCALED=F,Big_Window=F,window_scaling=c(1,0)){
  gamma=rho
  sigma=max(Boxs[1],sigma)
  sigma=min(Boxs[2],sigma)
  beta=min(BoxB[2],beta)
  beta=max(BoxB[1],beta)
  
  if(SF|SB){
    BaWe=2
  }else{
    BaWe=1
  }
  if(is.na(pop_vect)){
    pop_vect=rep(2,(n*0.5))
  }
  
  if(NC==1){
    M=dim(O)[1]-1
    L=as.numeric(O[(M+1),dim(O)[2]])
    Os=list()
    count=0
    theta_W=0
    s=dim(O)[1]
    
    for(k in 1:(M-1)){
      for(l in (k+1):M){
        Os_=seqSNP2_MH(O[c(k,l,s),],L)
        
        print(sum(Os_))
        if(sum(Os_)>=(L/(10^5))){#
          print(sum(Os_))
          count=count+1
          print(count)
          theta_W=theta_W+(sum(Os_))
          if(count==1){
            if(length(Os_)>10^6){
              Os_temporary=Build_zip_ID_2seq(Os_[1:1000000],max_count=symbol_number)
              Mat_symbol=Os_temporary[[2]]
              rm(Os_temporary)
              Os[[count]]=Zip_seq(Os_,Mat_symbol)
              
              
            }else{
              Os[[count]]=Build_zip_ID_2seq(Os_,max_count=symbol_number)
              Mat_symbol=Os[[count]][[2]]
            }
            
          }
          if(count>1){
            Os[[count]]=Zip_seq(Os_,Mat_symbol)
          }
          Os[[count]]=symbol2Num(Os[[count]][[1]],Os[[count]][[2]])
        }
      }
    }
    if(length(Os)>=1){
      theta_W=theta_W/length(Os)
      Ne=theta_W/(4*L*mu_r)
      if(log10(Ne)>N_max){
        beta=sqrt((10^(N_max))/Ne)
        if(beta<=BoxB[1]){
          print("Problem in Prior. Extra diversity cannot come from seed bank.")
          beta=BoxB[1]
        }
        if(beta>=BoxB[2]){
          beta=BoxB[2]
        }
      }
      theta=theta_W*(beta*beta)*2/((2-sigma)*(beta+((1-beta)*mu_b)))
      mu=theta/(2*L)
      rho=rho*theta
      
      O=Os
      rm(Os_)
    }
    if(length(Os)==0){
      stop("data too poor")
    }
    rm(Os)
    
  }
  if(NC>1){
    if(length(O)!=NC){
      stop("Not good number of chromosome given")
    }
    Os=list()
    theta_W_V=vector(length=NC)
    L_total=vector()
    for(chr in 1:NC){
      theta_W=0
      count=0
      OST_=list()
      M=(dim(O[[chr]])[1]-1)
      L=as.numeric(O[[chr]][(M+1),dim(O[[chr]])[2]])
      L_total=c(L_total,L)
      for(k in 1:(M-1)){
        for(l in (k+1):M){
          s=dim(O[[chr]])[1]
          Os_=seqSNP2_MH(O[[chr]][c(k,l,s),],L)
          if(sum(Os_)>=(L/(10^5))){
            count=count+1
            print(count)
            theta_W=theta_W+(sum(Os_))
            if(count==1){
              if(length(Os_)>10^6){
                Os_temporary=Build_zip_ID_2seq(Os_[1:1000000],max_count=symbol_number)
                Mat_symbol=Os_temporary[[2]]
                rm(Os_temporary)
                OST_[[count]]=Zip_seq(Os_,Mat_symbol)
                
                
              }else{
                OST_[[count]]=Build_zip_ID_2seq(Os_,max_count=symbol_number)
                Mat_symbol=OST_[[count]][[2]]
              }
              
            }
            if(count>1){
              OST_[[count]]=Zip_seq(Os_,Mat_symbol)
            }
            OST_[[count]]=symbol2Num(OST_[[count]][[1]],OST_[[count]][[2]])
          }
        }
      }
      
      rm(Os_)
      Os[[chr]]=OST_
      theta_W=theta_W/length(OST_)
      theta_W_V[chr]=theta_W
      rm(OST_)
    }
    
    Ne=mean(theta_W_V/(2*L_total*mu_r))
    if(log10(Ne)>N_max){
      beta=sqrt((10^(N_max))/Ne)
      if(beta<=BoxB[1]){
        print("Problem in Prior. Extra diversity cannot come from seed bank.")
        beta=BoxB[1]
      }
      if(beta>=BoxB[2]){
        beta=BoxB[2]
      }
    }
    theta=theta_W_V*(beta*beta)*2/((2-sigma)*(beta+((1-beta)*mu_b)))
    mu=theta/(2*L_total)
    rho=rho*theta
    O=Os
    L=L_total
    rm(Os)
    theta_W=theta_W_V
  }
  if(length(O)>=1){
    
    
    print(theta_W)
    
    
    if(BaWe==1){
      print('Searching for average recombination rate')
      
      results=Baum_Welch_algo(Os=O, maxIt =maxit,L=L,mu=mu,theta_W=theta_W,Rho=rho,beta=1,sigma=0,Popfix=pop,SB=F,SF=F,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=c(1,1),Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=T,NC=NC,BW=F,mu_b=mu_b,FS=FS,SCALED=T,Big_Window=Big_Window,window_scaling=window_scaling,redo_R=T)
      r=results$rho[1:NC]
      mu_=results$mu
      gamma_=r/mu_
      rho=gamma_*theta
      
      print('Searching for change of recombination rate in time')
      
      
      results=Baum_Welch_algo_t(Os=O, maxIt =maxit,L=L,mu=mu,theta_W =theta_W,Rho=rho,beta=beta,sigma=sigma,Popfix=pop,SB=SB,SF=SF,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=Boxr,Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=Rho,NC=NC,BW=BW,redo_R=T,mu_b=mu_b,SCALED=SCALED,Big_Window=Big_Window,window_scaling=window_scaling)
    }
    if(BaWe==2){
      if(SB){
      beta=1
    }
          if(SF){
      sigma=0
      }
         results=Baum_Welch_algo_t(Os=O, maxIt =maxit,L=L,mu=mu,theta_W =theta_W,Rho=rho,beta=beta,sigma=sigma,Popfix=pop,SB=SB,SF=SF,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=Boxr,Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=Rho,NC=NC,BW=BW,mu_b=mu_b,SCALED=SCALED,Big_Window=Big_Window,window_scaling=window_scaling)
    }
  }
  
  return(results)
}



teSMC_old<-function(n=40,rho=1,O,mu_r,maxit =20,symbol_number=30,BoxB=c(0.05,1),BoxP=c(2,2),Boxr=c(1,1),Boxs=c(0,0.99),N_max=7,pop=F,SB=F,SF=F,Rho=T,Check=F,BW=F,FS=F,NC=1,mu_b=1,pop_vect=NA,rec_vect=NA,sigma=0,beta=1,SCALED=F,Big_Window=F,window_scaling=c(1,0)){
  gamma=rho
  sigma=max(Boxs[1],sigma)
  sigma=min(Boxs[2],sigma)
  beta=min(BoxB[2],beta)
  beta=max(BoxB[1],beta)
  
  if(SF|SB){
    BaWe=2
  }else{
    BaWe=1
  }
  if(is.na(pop_vect)){
    pop_vect=rep(2,(n*0.5))
  }
  
  if(NC==1){
    M=dim(O)[1]-1
    L=as.numeric(O[(M+1),dim(O)[2]])
    Os=list()
    count=0
    theta_W=0
    s=dim(O)[1]
    
    for(k in 1:(M-1)){
      for(l in (k+1):M){
        Os_=seqSNP2_MH(O[c(k,l,s),],L)
        
        print(sum(Os_))
        if(sum(Os_)>=(L/(10^5))){#
          print(sum(Os_))
          count=count+1
          print(count)
          theta_W=theta_W+(sum(Os_))
          if(count==1){
            if(length(Os_)>10^6){
              Os_temporary=Build_zip_ID_2seq(Os_[1:1000000],max_count=symbol_number)
              Mat_symbol=Os_temporary[[2]]
              rm(Os_temporary)
              Os[[count]]=Zip_seq(Os_,Mat_symbol)
              
              
            }else{
              Os[[count]]=Build_zip_ID_2seq(Os_,max_count=symbol_number)
              Mat_symbol=Os[[count]][[2]]
            }
            
          }
          if(count>1){
            Os[[count]]=Zip_seq(Os_,Mat_symbol)
          }
          Os[[count]]=symbol2Num(Os[[count]][[1]],Os[[count]][[2]])
        }
      }
    }
    if(length(Os)>=1){
      theta_W=theta_W/length(Os)
      Ne=theta_W/(4*L*mu_r)
      if(log10(Ne)>N_max){
        beta=sqrt((10^(N_max))/Ne)
        if(beta<=BoxB[1]){
          print("Problem in Prior. Extra diversity cannot come from seed bank.")
          beta=BoxB[1]
        }
        if(beta>=BoxB[2]){
          beta=BoxB[2]
        }
      }
      theta=theta_W*(beta*beta)*2/((2-sigma)*(beta+((1-beta)*mu_b)))
      mu=theta/(2*L)
      rho=rho*theta
      
      O=Os
      rm(Os_)
    }
    if(length(Os)==0){
      stop("data too poor")
    }
    rm(Os)
    
  }
  if(NC>1){
    if(length(O)!=NC){
      stop("Not good number of chromosome given")
    }
    Os=list()
    theta_W_V=vector(length=NC)
    L_total=vector()
    for(chr in 1:NC){
      theta_W=0
      count=0
      OST_=list()
      M=(dim(O[[chr]])[1]-1)
      L=as.numeric(O[[chr]][(M+1),dim(O[[chr]])[2]])
      L_total=c(L_total,L)
      for(k in 1:(M-1)){
        for(l in (k+1):M){
          s=dim(O[[chr]])[1]
          Os_=seqSNP2_MH(O[[chr]][c(k,l,s),],L)
          if(sum(Os_)>=(L/(10^5))){
            count=count+1
            print(count)
            theta_W=theta_W+(sum(Os_))
            if(count==1){
              if(length(Os_)>10^6){
                Os_temporary=Build_zip_ID_2seq(Os_[1:1000000],max_count=symbol_number)
                Mat_symbol=Os_temporary[[2]]
                rm(Os_temporary)
                OST_[[count]]=Zip_seq(Os_,Mat_symbol)
                
                
              }else{
                OST_[[count]]=Build_zip_ID_2seq(Os_,max_count=symbol_number)
                Mat_symbol=OST_[[count]][[2]]
              }
              
            }
            if(count>1){
              OST_[[count]]=Zip_seq(Os_,Mat_symbol)
            }
            OST_[[count]]=symbol2Num(OST_[[count]][[1]],OST_[[count]][[2]])
          }
        }
      }
      
      rm(Os_)
      Os[[chr]]=OST_
      theta_W=theta_W/length(OST_)
      theta_W_V[chr]=theta_W
      rm(OST_)
    }
    
    Ne=mean(theta_W_V/(2*L_total*mu_r))
    if(log10(Ne)>N_max){
      beta=sqrt((10^(N_max))/Ne)
      if(beta<=BoxB[1]){
        print("Problem in Prior. Extra diversity cannot come from seed bank.")
        beta=BoxB[1]
      }
      if(beta>=BoxB[2]){
        beta=BoxB[2]
      }
    }
    theta=theta_W_V*(beta*beta)*2/((2-sigma)*(beta+((1-beta)*mu_b)))
    mu=theta/(2*L_total)
    rho=rho*theta
    O=Os
    L=L_total
    rm(Os)
    theta_W=theta_W_V
  }
  if(length(O)>=1){
    
    
    print(theta_W)
    
    
    if(BaWe==1){
      print('Searching for average recombination rate')
      
      results=Baum_Welch_algo(Os=O, maxIt =maxit,L=L,mu=mu,theta_W=theta_W,Rho=rho,beta=1,sigma=0,Popfix=pop,SB=F,SF=F,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=c(1,1),Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=T,NC=NC,BW=F,mu_b=mu_b,FS=FS,SCALED=T,Big_Window=Big_Window,window_scaling=window_scaling,redo_R=T)
      r=results$rho[1:NC]
      mu_=results$mu
      gamma_=r/mu_
      rho=gamma_*theta
      
      print('Searching for change of recombination rate in time')
      
      
      results=Baum_Welch_algo_t(Os=O, maxIt =maxit,L=L,mu=mu,theta_W =theta_W,Rho=rho,beta=beta,sigma=sigma,Popfix=pop,SB=SB,SF=SF,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=Boxr,Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=Rho,NC=NC,BW=BW,redo_R=T,mu_b=mu_b,SCALED=SCALED,Big_Window=Big_Window,window_scaling=window_scaling)
    }
    if(BaWe==2){
      print('Searching for average recombination rate')
      
      results=Baum_Welch_algo(Os=O, maxIt =maxit,L=L,mu=mu,theta_W=theta_W,Rho=theta,beta=1,sigma=0,Popfix=pop,SB=F,SF=F,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=c(1,1),Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=T,NC=NC,BW=F,mu_b=mu_b,FS=FS,SCALED=T,Big_Window=Big_Window,window_scaling=window_scaling,redo_R=T)
      r=results$rho[1:NC]
      mu_=results$mu
      gamma_=r/mu_
      rho=gamma_*theta
      
      print(r)
      print(mu_)
      print(gamma_)
      if(mean(gamma_)>1){
        print("Results might not be reliable")
      }
      effect=mean(gamma_/gamma)
      if(SF&!SB){
        sigma=(1-effect)/(1-(effect/2))
        sigma=max(Boxs[1],sigma)
        sigma=min(Boxs[2],sigma)
        if(sigma>=Boxs[1]&sigma<=Boxs[2]){
          Boxs[1]=sigma
        }
      }
      if(SB&!SF){
        beta=effect
        beta=min(BoxB[2],beta)
        beta=max(BoxB[1],beta)
        if(beta>=BoxB[1]&beta<=BoxB[2]){
          BoxB[2]=beta
        }
      }
      if(SF&SB){
        if(min(BoxB)>(1-max(Boxs))){
          sigma=(1-(effect/beta))/(1-(effect/(2*beta)))
          sigma=max(Boxs[1],sigma)
          sigma=min(Boxs[2],sigma)
          beta=effect/(1-sigma)
          beta=max(BoxB[1],beta)
          beta=min(BoxB[2],beta)
        }
        if(min(BoxB)<=(1-max(Boxs))){
          sigma=(1-(effect/beta))/(1-(effect/(2*beta)))
          sigma=max(Boxs[1],sigma)
          sigma=min(Boxs[2],sigma)
          beta=effect/(1-sigma)
          beta=max(BoxB[1],beta)
          beta=min(BoxB[2],beta)
        }
        if(gamma_!=(gamma*beta*2*(1-sigma)/(2-sigma))){
          print("Prior might disagree with results.")
        }
      }
      theta=theta_W*(beta*beta)*2/((2-sigma)*(beta+((1-beta)*mu_b)))
      mu=theta/(2*L)
      rho=gamma*theta
      results=Baum_Welch_algo_t(Os=O, maxIt =maxit,L=L,mu=mu,theta_W =theta_W,Rho=rho,beta=beta,sigma=sigma,Popfix=pop,SB=SB,SF=SF,k=n,BoxB=BoxB,BoxP=BoxP,Boxr=Boxr,Boxs = Boxs,maxBit = 1,pop_vect=pop_vect,ER=Rho,NC=NC,BW=BW,mu_b=mu_b,SCALED=SCALED,Big_Window=Big_Window,window_scaling=window_scaling)
    }
  }
  
  return(results)
}

####################
# Ziping functions #
####################

# Zip the  signal to reduce signal length
# Os : original signal
# max_count : maximum number of symbol
Build_zip_ID_2seq<-function(Os,max_count){
  count=0
  ini=Os[1]
  Os=Os[-1]
  output=list()
  criteria=T
  mat_symbol=vector()
  l=vector()
  while(criteria){
    count=count+1
    if(count<9){
      new_symbol=1+count
    }
    if(count>8& (count-8)<27){
      new_symbol=letters[count-8]
    }
    if((count-8)>=27){
      criteria=F
    }
    if((count-8)<27){
      L=length(Os)
      impair=cbind(seq(1,(L-1),2),seq(2,L,2))
      pair=cbind(seq(2,(L-1),2),seq(3,L,2))
      tot_pair=rbind(impair,pair)
      rm(impair)
      rm(pair)
      pattern=rbind(as.vector(as.matrix(paste(Os[tot_pair[,1]],Os[tot_pair[,2]],sep=""))),tot_pair[,1],tot_pair[,2])
      rm(tot_pair)
      symbol=unique(pattern[1,])
      real_symb=vector()
      nb_symbol=length(symbol)
      for(k in 1:nb_symbol){
        sym_dec=strsplit(symbol[k],"")[[1]]
        if(all(as.character(sym_dec)!="1")){
          real_symb=c(real_symb,symbol[k])
        }
      }
      symbol=real_symb
      if(length(symbol)==0){
        criteria=F
      }
      if(length(symbol)>0){
        nb_symbol=length(symbol)
        count_symbol=rep(0,nb_symbol)
        pos_s=vector()
        for(k in 1:nb_symbol){
          test=as.numeric(as.matrix(pattern[2,which(pattern[1,]==symbol[k])]))
          x=length(test)
          test=c(test,c(test+1))
          coef=length(unique(test))/length(test)
          count_symbol[k]=x*coef
        }
        rm(pattern)
        rm(test)
        rep_symbol=as.vector(symbol[which(count_symbol==max(count_symbol))])[1]
        pos_s=gregexpr(rep_symbol,paste(as.character(Os),collapse = ""))
        pos_s=as.vector(as.matrix(pos_s[[1]]))
        mat_symbol=rbind(mat_symbol,c(new_symbol,rep_symbol),l[length(l)])
        Os[pos_s]=new_symbol
        pos_s=pos_s+1
        Os=Os[-pos_s]
        if(count>=max_count)
        {
          criteria=F
        }
      }
    }
  }
  Os=c(ini,Os)
  output[[1]]=Os
  output[[2]]=mat_symbol
  return(output)
}

# Builds a numerical zip signal from the zip signal
# Os : signal
# Mat_symbol : matrix with the zipping pattern
symbol2Num<-function(Os,Mat_symbol,count=0){
  for(i in 1:dim(Mat_symbol)[1]){
    sym=strsplit(Mat_symbol[i,2],"")
    sym=sym[[1]]
    sym=paste(sym,collapse=" ")
    Mat_symbol[i,2]=sym
  }
  
  for( i in letters){
    if(any(Mat_symbol[,1]%in%i)){
      pos=which(Os%in%i)
      x=as.numeric(which(letters%in%i))
      Os[pos]=(9+x)
      pos_i=which(Mat_symbol[,1]==i)
      Mat_symbol[pos_i,1]=9+x
      sym=strsplit(Mat_symbol[pos_i,2]," ")
      sym=sym[[1]]
      if(sym[1]%in%letters){
        sym[1]=9+as.numeric(which(letters%in%sym[1]))
      }
      if(sym[2]%in%letters){
        sym[2]=9+as.numeric(which(letters%in%sym[2]))
      }
      sym=paste(sym,collapse=" ")
      Mat_symbol[pos_i,2]=sym
    }
  }
  output=list()
  output[[1]]=Os
  output[[2]]=Mat_symbol
  mat_length=c(as.character(0:(count+1)),Mat_symbol[,1])
  mat_length=cbind(mat_length,rep(0,length(mat_length)))
  l=rep(1,(count+2))
  for(i in 1:dim(Mat_symbol)[1]){
    sx=strsplit(Mat_symbol[i,2]," ")
    sx=as.numeric(as.matrix(sx[[1]]))
    a=sx[1]
    b=sx[2]
    new_l=l[a+1]+l[b+1]
    l=c(l,new_l)
  }
  mat_length[,2]=l
  output[[3]]=mat_length
  return(output)
}

# Builds a zip signal from the signal
# Os : signal
# Mat_symbol : matrix with the zipping pattern
Zip_seq<-function(Os,Mat_symbol){
  ini=Os[1]
  Os=Os[-1]
  output=list()
  criteria=T
  l=vector()
  for(jj in 1:dim(Mat_symbol)[1]){
    new_symbol= Mat_symbol[jj,1]
    L=length(Os)
    rep_symbol= Mat_symbol[jj,2]
    pos_s=gregexpr(rep_symbol,paste(as.character(Os),collapse = ""))
    pos_s=as.vector(as.matrix(pos_s[[1]]))
    if(length(pos_s)>1){
      Os[pos_s]=new_symbol
      pos_s=pos_s+1
      Os=Os[-pos_s]
    }
  }
  Os=c(ini,Os)
  output[[1]]=Os
  output[[2]]=Mat_symbol
  return(output)
}

# Builds the HMM matrices to analysed the zip signal
# A : Transition matrix
# E : emission matrix
# Mat_symbol : matrix with the zipping pattern
# q : equilibrium probability
Build_zip_Matrix_mailund<-function(A,E,Mat_symbol,q){
  C=list()
  TO=list()
  l=vector(length =(2+(dim(Mat_symbol)[1])))
  Q=list()
  Q_=list()
  k=length(q)
  
  for(i in 1:dim(E)[2]){
    B=diag(x=E[,i])
    C[[i]]=B%*%A
    l[i]=1
    TO[[i]]=t(A)%*%B
    Q[[i]]=diag(1,k,k)
    Q_[[i]]=diag(1,k,k)
  }
  x=length(C)
  mat_trick=eigen(TO[[1]])
  D=diag(mat_trick$values)
  for(i in 1:dim(Mat_symbol)[1]){
    sx=strsplit(Mat_symbol[i,2]," ")
    sx=as.numeric(as.matrix(sx[[1]]))
    a=sx[1]
    b=sx[2]
    C[[(x+i)]]=as.matrix(C[[(b+1)]])%*%as.matrix(C[[(a+1)]])
    l[(i+2)]=(l[(a+1)]+l[(b+1)])
    l_temp=(l[(a+1)]+l[(b+1)])
    TO[[x+i]]=TO[[(a+1)]]%*%TO[[(b+1)]]
    Q_temp=matrix(0,dim(D)[1],dim(D)[1])
    Q_temp_=matrix(0,dim(D)[1],dim(D)[1])
    for(x1 in 1:dim(D)[1]){
      for(x2 in 1:dim(D)[1]){
        if(D[x1,x1]!=D[x2,x2]){
          Q_temp[x1,x2]=((D[x1,x1]^(l_temp))-(D[x2,x2]^(l_temp)))/(D[x1,x1]-D[x2,x2] )
          Q_temp_[x1,x2]=sum((D[x1,x1]^seq(1,l_temp,1))*(D[x2,x2]^(l_temp-seq(1,l_temp,1))))
        }
        if(D[x1,x1]==D[x2,x2]){
          Q_temp[x1,x2]=(l_temp)*(D[x2,x2]^(l_temp-1))
          Q_temp_[x1,x2]=(l_temp)*(D[x2,x2]^(l_temp))
        }
      }
    }
    Q[[i+2]]=Q_temp
    Q_[[i+2]]=Q_temp_
  }
  output=list()
  output[[1]]=C
  output[[2]]=Q_
  output[[3]]=TO
  output[[4]]=Q
  
  return(output)
}


# Zip the  signal to reduce signal length
# Os : original signal
# max_count : maximum number of symbol
Build_zip_2seq<-function(Os,max_count){
  ini=Os[1]
  Os=Os[-1]
  count=0
  output=list()
  criteria=T
  mat_symbol=vector()
  l=vector()
  while(criteria){
    count=count+1
    #print(paste("count : ",count))
    if(count<9){
      new_symbol=1+count
    }
    if(count>8& (count-8)<27){
      new_symbol=letters[count-8]
    }
    L=length(Os)
    impair=cbind(seq(1,(L-1),2),seq(2,L,2))
    pair=cbind(seq(2,(L-1),2),seq(3,L,2))
    tot_pair=rbind(impair,pair)
    rm(impair)
    rm(pair)
    pattern=rbind(as.vector(as.matrix(paste(Os[tot_pair[,1]],Os[tot_pair[,2]],sep=""))),tot_pair[,1],tot_pair[,2])
    rm(tot_pair)
    symbol=unique(pattern[1,])
    nb_symbol=length(symbol)
    count_symbol=rep(0,nb_symbol)
    pos_s=vector()
    for(k in 1:nb_symbol){
      test=as.numeric(as.matrix(pattern[2,which(pattern[1,]==symbol[k])]))
      x=length(test)
      test=c(test,c(test+1))
      coef=length(unique(test))/length(test)
      count_symbol[k]=x*coef
    }
    
    rm(pattern)
    rm(test)
    rep_symbol=as.vector(symbol[which(count_symbol==max(count_symbol))])[1]
    pos_s=gregexpr(rep_symbol,paste(as.character(Os),collapse = ""))
    pos_s=as.vector(as.matrix(pos_s[[1]]))
    #  print(paste("replaced symbol: ",rep_symbol))
    #  print(paste("news symbol: ",new_symbol))
    
    mat_symbol=rbind(mat_symbol,c(new_symbol,rep_symbol),l[length(l)])
    Os[pos_s]=new_symbol
    pos_s=pos_s+1
    Os=Os[-pos_s]
    
    if(count>=max_count)#Check_time(M,i,ti))
    {
      criteria=F
    }
    
  }
  Os=c(ini,Os)
  output[[1]]=Os
  output[[2]]=mat_symbol
  return(output)
}



Unzip_2seq<-function(Os){
  mat_symbol=Os[[2]]
  Os=Os[[1]]
  for(i in seq(dim(mat_symbol)[1],1,-1)){
    symbol=mat_symbol[i,1]
    pos=which(Os%in%symbol)
    Os[pos]=mat_symbol[i,2]
    Os=paste(Os,collapse = "")
    Os=strsplit(Os,"")
    Os=as.matrix(Os[[1]])
  }
  return(Os)
}



#####################
# Forward algorithm #
#####################

# Does the forward alagorithm :
# Os: zip signal
# E : emission matrix
# q : equilibrium probability
# C : matrix for the zip signal
forward_zip_mailund <- function(Os,E,q,C){
  Os=as.numeric(as.vector(Os))
  T_prime=length(Os)
  output=list()
  alpha=matrix(0,ncol=T_prime,nrow=length(q))
  d=numeric(length = T_prime)
  alpha[,1]=diag(x=E[,(1+as.numeric(Os[1]))])%*%q
  d_n=sum(alpha[,1])
  alpha[,1]=alpha[,1]/sum(alpha[,1])
  d[1]=log(d_n)
  explore=2:T_prime
  for(i in explore){
    
    truc=C[[(1+as.numeric(Os[i]))]]%*%alpha[,(i-1)]
    
    dn=sum(truc)
    alpha[,i]=truc/dn
    d[i]=log(dn)
  }
  output[[1]]=alpha
  output[[2]]=d
  output[[3]]=(sum(d))
  
  return(output)
}





######################
#Backward algorithm ##
######################

# does the backward algorithm
# y : zip signal
# TO : matrix for the zip signal
# n : number of  hidden states
# c : log-lh outputed by the forward algorihtm
Backward_zip_mailund<-function(y,TO,n,c){
  dims <- n;
  n<- length(y);
  beta <- matrix(1, nrow=dims[1], ncol=n);
  d=numeric(n)
  for (t in seq((n-1), 1, -1)){
    beta[,t] = (TO[[(1+as.numeric(y[(t+1)]))]]%*%beta[, (t+1)])/c[t+1];
  }
  return(beta);
  
}


# does the backward algorithm
# y : zip signal
# TO : matrix for the zip signal
# n : number of  hidden states
# c : log-lh outputed by the forward algorihtm
Backward_zip_mailund_phased<-function(y,TO,n,c,phasing,A,E){
  dims <- n;
  n<- length(y);
  beta <- matrix(1, nrow=dims[1], ncol=n);
  d=numeric(n)
  count_2=length(phasing)
  B1=diag(E[,1])
  B2=diag(E[,2])
  for (t in seq((n-1), 1, -1)){
    if(as.numeric(y[(t+1)])!=2){
      beta[,t] = (TO[[(1+as.numeric(y[(t+1)]))]]%*%beta[, (t+1)])/c[t+1];
    }else{
      beta[,t] = ((phasing[count_2]*((t(A)%*%B2)%*%beta[,(t+1)]))+((1-phasing[count_2])*((t(A)%*%B1)%*%beta[,(t+1)])))/c[t+1];
      count_2=count_2-1
    }
  }
  return(beta);
  
}


####################
# Build HMM matrix #
####################

# build all the required matrices for the  HMM
build_HMM_matrix<-function(n=20,rho,beta=1,L,Pop=T,Xi=NA,Beta=1,scale=c(1,0),sigma=0,Sigma=0,FS=F,Big_Window=F,tmax=15,alpha_t=0.1,npair=1,cut_edge=F){
  FS=F
  if(cut_edge){
    n=n+1
  }
  if(as.numeric(Big_Window)==0){
    Vect=0:(n-1)
    Tc= scale[2] -(0.5*(2-Sigma)*log(1-(Vect/n))/((Beta^2)*scale[1]))
  }
  if(as.numeric(Big_Window)==1){
    Vect=1:(n-1)
    alpha_t=alpha_t/(npair)
    tmax=tmax*npair
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])))
    
  }
  if(as.numeric(Big_Window)==2){
    tmax=50
    Vect=1:(n-1)
    alpha_t=0.001
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])) )
  }
  if(as.numeric(Big_Window)==3){
    tmax=20
    Vect=1:(n-1)
    alpha_t=0.1
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])) )
  }
  
  if(as.numeric(Big_Window)==4){
    tmax=20
    alpha_t=0.005
    Vect=1:(n-1)
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])) )
  }
  
  t=vector()
  q=vector()
  r=rho/(2*(L-1))
  D=Tc[2:n]-Tc[1:(n-1)]
  if(Pop){
    Xi=rep(1,n)
  }
  t[1:(n-1)]=(Xi[1:(n-1)]*(2-sigma)/(beta*beta*2))+ (Tc[1:(n-1)]-(Tc[2:n]*(exp(-beta*beta*(D[1:(n-1)]/Xi[1:(n-1)])*2/(2-sigma)))))/(1-exp(-beta*beta*((D[1:(n-1)])/Xi[1:(n-1)])*2/(2-sigma)))
  t[n]=(Xi[n]*(2-sigma)/(beta*beta*2)) +Tc[n]
  q[1]=(1-exp(-D[1]*2*beta*beta/(Xi[1]*(2-sigma))))
  for(alpha in 2:(n-1)){q[alpha]=exp(-(2*beta*beta/(2-sigma))*sum(D[1:(alpha-1)]/Xi[1:(alpha-1)]))*(1-exp(-D[alpha]*2*beta*beta/(Xi[alpha]*(2-sigma))))}
  q[n]=exp(-(2*beta*beta/(2-sigma))*sum(D[1:(n-1)]/Xi[1:(n-1)]))
  Q=matrix(0,ncol = n, nrow = n)
  if(!FS){
    for(i in 1:(n-1)){
      if(i==1){
        Q[1,(2:n)]=(1-exp(-2*r*t[(2:n)]*beta*2*(1-sigma)/(2-sigma)))*(1/(2*t[(2:n)]))*( D[1]-((1-exp(-D[1]*beta*beta*4/(Xi[1]*(2-sigma))))/(4*beta*beta/(Xi[1]*(2-sigma)))))
      }
      if(i==2){
        eta=1
        truc= ((1-exp(-D[eta]*beta*beta*4/(Xi[eta]*(2-sigma))))/((4*beta*beta)/(Xi[eta]*(2-sigma))))
        truc=truc*(1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))
        Q[i,((i+1):n)]=(1-exp(-2*r*t[((i+1):n)]*beta*2*(1-sigma)/(2-sigma)))*(1/(2*t[((i+1):n)]))*( truc + (D[i]-((1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))/(4*beta*beta/(Xi[i]*(2-sigma))))))
        Q[i,1]=(1-exp(-2*r*t[1]*beta*2*(1-sigma)/(2-sigma)))*(1/t[1])*exp((t[1]-Tc[2])*beta*beta*2/(Xi[1]*(2-sigma)))*(1-exp(-D[2]*beta*beta*2/(Xi[2]*(2-sigma))))*(((1-exp(-t[1]*4*beta*beta/(Xi[1]*(2-sigma))))/(4*beta*beta/(Xi[1]*(2-sigma)))))
      }
      if( i>2){
        truc=0
        for(eta in 1:(i-2)){
          truc= truc+ ((( exp(-4*beta*beta*sum(D[(eta+1):(i-1)]/(Xi[(eta+1):(i-1)]*(2-sigma))))*(1-exp(-D[eta]*beta*beta*4/(Xi[eta]*(2-sigma))))/((4*beta*beta)/(Xi[eta]*(2-sigma))))))
        }
        eta=i-1
        truc= truc+ ((((1-exp(-D[eta]*beta*beta*4/(Xi[eta]*(2-sigma))))/((4*beta*beta)/(Xi[eta]*(2-sigma))))))
        truc=truc*(1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))
        Q[i,((i+1):n)]=(1-exp(-2*r*t[((i+1):n)]*beta*2*(1-sigma)/(2-sigma)))*(1/(2*t[((i+1):n)]))*( truc + ((D[i]-((1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))/(4*beta*beta/(Xi[i]*(2-sigma)))))))
        truc=rep(0,length(1:(i-1)))
        exp_truc=rep(0,length(1:(i-1)))
        for(gamma in 1:(i-1)){
          if(gamma<(i-1)){
            exp_truc[gamma]=exp(-beta*beta*2*(sum(D[(gamma+1):(i-1)]/(Xi[(gamma+1):(i-1)]*(2-sigma)))+((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
          }
          if(gamma==(i-1)){
            exp_truc[gamma]=exp(-beta*beta*2*(((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
          }
          if(gamma>2){
            sub_truc=0
            for(eta in 1:(gamma-2)){
              sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(sum(D[(eta+1):(gamma-1)]/(Xi[(eta+1):(gamma-1)]*(2-sigma)))+((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
            }
            eta=gamma-1
            sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
            truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
          }
          if(gamma==2){
            sub_truc=0
            eta=1
            sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma)))))
            truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
          }
          if(gamma==1){
            truc[gamma]=((1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
          }
        }
        Q[i,1:(i-1)]=(1-exp(-2*r*t[(1:(i-1))]*beta*2*(1-sigma)/(2-sigma)))*(1/t[(1:(i-1))])*(1-exp(-D[i]*beta*beta*2/(Xi[i]*(2-sigma))))*exp_truc*truc
      }
    }
    truc=rep(0,length(1:(n-1)))
    exp_truc=rep(0,length(1:(n-1)))
    for(gamma in 1:(n-1)){
      if(gamma<(n-1)){
        exp_truc[gamma]=exp(-beta*beta*2*( sum(D[(gamma+1):(n-1)]/(Xi[(gamma+1):(n-1)]*(2-sigma)) ) + ((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
      }
      if(gamma==(n-1)){
        exp_truc[gamma]=exp(-beta*beta*2*(((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
      }
      if(gamma>2){
        sub_truc=0
        for(eta in 1:(gamma-2)){
          sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(sum(D[(eta+1):(gamma-1)]/(Xi[(eta+1):(gamma-1)]*(2-sigma))) + ((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
        }
        eta=gamma-1
        sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
        truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
      }
      if(gamma==2){
        sub_truc=0
        eta=1
        sub_truc=sub_truc+(( (1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
        truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
      }
      if(gamma==1){
        truc[gamma]=((1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
      }
    }
    Q[n,1:(n-1)]=(1-exp(-2*r*t[(1:(n-1))]*beta*2*(1-sigma)/(2-sigma)))*(1/t[(1:(n-1))])*exp_truc*truc
    diag(Q)=rep(1,n)-apply(Q,2,sum)
  }
  if(FS){
    if(F){
      e_t=0
      for(i in 1:n){
        if(i>1){
          truc_exp=exp(-(2*beta*beta/(2-sigma))*sum(D[1:(i-1)]/Xi[1:(i-1)]))
        }else{
          truc_exp=1
        }
        if(i<n){
          Lambda_gamma=beta*beta*2/(Xi[i]*(2-sigma))
          e_t=e_t+((truc_exp*( ( Tc[i]*Lambda_gamma+1)-((Tc[i+1]*Lambda_gamma+1)*exp(-Lambda_gamma*D[i])  )  ))/Lambda_gamma)
        }else{
          truc_exp=exp(-(2*beta*beta/(2-sigma))*sum(D[1:(i-1)]/Xi[1:(i-1)]))
          e_t=e_t+((truc_exp*( ( Tc[i]*Lambda_gamma+1)  ))/Lambda_gamma)
        }
        
        
      }
    }
    
    
    k=ceiling(1/(2*r*t[n]*beta*2*(1-sigma)/(2-sigma)))
    tot_prec=numeric(n)
    p_rec=list()
    for(k_t in 0:k){
      p_rec_t=((2*r*t*beta*2*(1-sigma)/(2-sigma))^(k_t))*exp(-2*r*t*beta*2*(1-sigma)/(2-sigma))/factorial(k_t)
      tot_prec=p_rec_t+tot_prec
      p_rec[[(k_t+1)]]=p_rec_t
    }
    for(k_t in 1:k+1){
      p_rec[[k_t]]=p_rec[[k_t]]/tot_prec
    }
    
    #diag(Q)=p_rec[[1]]
    Q_base=Transition_Mat(n,rho,beta,L,Pop,Xi,Beta,scale,sigma,Sigma)
    for(k_t in 1:k){
      Q=Q+( diag(p_rec[[(k_t+1)]])%*%(Q_base%^%k_t))
    }
    #corrector_Q=rowSums(Q)
    #Q=diag((1-p_rec[[1]])/corrector_Q)%*%Q
    diag(Q)=p_rec[[1]]+diag(Q)
    #stop()
    #print(sum(Q))
  }
  if(cut_edge){
    Tc=Tc[-n]
    q=q[-n]
    q=q/sum(q)
    t=t[-n]
    Q=Q[-n,-n]
    n=n-1
  }
  output=list()
  output[[1]]=Q
  output[[2]]=q
  output[[3]]=t
  output[[4]]=Tc
  return(output)
}

Transition_Mat<-function(n,rho,beta,L,Pop,Xi,Beta,scale,sigma,Sigma,Big_Window=F,tmax=15,alpha_t=0.1,npair=2,cut_edge=F){
  if(cut_edge){
    n=n+1
  }
  if(!Big_Window){
    Vect=0:(n-1)
    alpha_t=alpha_t/npair
    Tc= scale[2] -(0.5*(2-Sigma)*log(1-(Vect/n))/((Beta^2)*scale[1]))
  }
  if(Big_Window){
    Vect=1:(n-1)
    alpha_t=alpha_t/npair
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])) )
    
  }
  t=vector()
  q=vector()
  r=rho/(2*(L-1))
  D=Tc[2:n]-Tc[1:(n-1)]
  if(Pop){
    Xi=rep(1,n)
  }
  t[1:(n-1)]=(Xi[1:(n-1)]*(2-sigma)/(beta*beta*2))+ (Tc[1:(n-1)]-(Tc[2:n]*(exp(-beta*beta*(D[1:(n-1)]/Xi[1:(n-1)])*2/(2-sigma)))))/(1-exp(-beta*beta*((D[1:(n-1)])/Xi[1:(n-1)])*2/(2-sigma)))
  t[n]=(Xi[n]*(2-sigma)/(beta*beta*2)) +Tc[n]
  q[1]=(1-exp(-D[1]*2*beta*beta/(Xi[1]*(2-sigma))))
  for(alpha in 2:(n-1)){q[alpha]=exp(-(2*beta*beta/(2-sigma))*sum(D[1:(alpha-1)]/Xi[1:(alpha-1)]))*(1-exp(-D[alpha]*2*beta*beta/(Xi[alpha]*(2-sigma))))}
  q[n]=exp(-(2*beta*beta/(2-sigma))*sum(D[1:(n-1)]/Xi[1:(n-1)]))
  Q=matrix(0,ncol = n, nrow = n)
  for(i in 1:(n-1)){
    if(i==1){
      Q[1,(2:n)]=(1/(2*t[(2:n)]))*( D[1]-((1-exp(-D[1]*beta*beta*4/(Xi[1]*(2-sigma))))/(4*beta*beta/(Xi[1]*(2-sigma)))))
    }
    if(i==2){
      eta=1
      truc= ((1-exp(-D[eta]*beta*beta*4/(Xi[eta]*(2-sigma))))/((4*beta*beta)/(Xi[eta]*(2-sigma))))
      truc=truc*(1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))
      Q[i,((i+1):n)]=(1/(2*t[((i+1):n)]))*( truc + (D[i]-((1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))/(4*beta*beta/(Xi[i]*(2-sigma))))))
      Q[i,1]=(1/t[1])*exp((t[1]-Tc[2])*beta*beta*2/(Xi[1]*(2-sigma)))*(1-exp(-D[2]*beta*beta*2/(Xi[2]*(2-sigma))))*(((1-exp(-t[1]*4*beta*beta/(Xi[1]*(2-sigma))))/(4*beta*beta/(Xi[1]*(2-sigma)))))
    }
    if( i>2){
      truc=0
      for(eta in 1:(i-2)){
        truc= truc+ ((( exp(-4*beta*beta*sum(D[(eta+1):(i-1)]/(Xi[(eta+1):(i-1)]*(2-sigma))))*(1-exp(-D[eta]*beta*beta*4/(Xi[eta]*(2-sigma))))/((4*beta*beta)/(Xi[eta]*(2-sigma))))))
      }
      eta=i-1
      truc= truc+ ((((1-exp(-D[eta]*beta*beta*4/(Xi[eta]*(2-sigma))))/((4*beta*beta)/(Xi[eta]*(2-sigma))))))
      truc=truc*(1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))
      Q[i,((i+1):n)]=(1/(2*t[((i+1):n)]))*( truc + ((D[i]-((1-exp(-D[i]*beta*beta*4/(Xi[i]*(2-sigma))))/(4*beta*beta/(Xi[i]*(2-sigma)))))))
      truc=rep(0,length(1:(i-1)))
      exp_truc=rep(0,length(1:(i-1)))
      for(gamma in 1:(i-1)){
        if(gamma<(i-1)){
          exp_truc[gamma]=exp(-beta*beta*2*(sum(D[(gamma+1):(i-1)]/(Xi[(gamma+1):(i-1)]*(2-sigma)))+((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
        }
        if(gamma==(i-1)){
          exp_truc[gamma]=exp(-beta*beta*2*(((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
        }
        if(gamma>2){
          sub_truc=0
          for(eta in 1:(gamma-2)){
            sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(sum(D[(eta+1):(gamma-1)]/(Xi[(eta+1):(gamma-1)]*(2-sigma)))+((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
          }
          eta=gamma-1
          sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
          truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
        }
        if(gamma==2){
          sub_truc=0
          eta=1
          sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma)))))
          truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
        }
        if(gamma==1){
          truc[gamma]=((1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
        }
      }
      Q[i,1:(i-1)]=(1/t[(1:(i-1))])*(1-exp(-D[i]*beta*beta*2/(Xi[i]*(2-sigma))))*exp_truc*truc
    }
  }
  truc=rep(0,length(1:(n-1)))
  exp_truc=rep(0,length(1:(n-1)))
  for(gamma in 1:(n-1)){
    if(gamma<(n-1)){
      exp_truc[gamma]=exp(-beta*beta*2*( sum(D[(gamma+1):(n-1)]/(Xi[(gamma+1):(n-1)]*(2-sigma)) ) + ((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
    }
    if(gamma==(n-1)){
      exp_truc[gamma]=exp(-beta*beta*2*(((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma)))))
    }
    if(gamma>2){
      sub_truc=0
      for(eta in 1:(gamma-2)){
        sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(sum(D[(eta+1):(gamma-1)]/(Xi[(eta+1):(gamma-1)]*(2-sigma))) + ((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
      }
      eta=gamma-1
      sub_truc=sub_truc+(((1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
      truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
    }
    if(gamma==2){
      sub_truc=0
      eta=1
      sub_truc=sub_truc+(( (1-exp(-D[eta]*4*beta*beta/(Xi[eta]*(2-sigma))))/(4*beta*beta/(Xi[eta]*(2-sigma))))*exp(-4*beta*beta*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma))))))
      truc[gamma]=sub_truc+( (1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
    }
    if(gamma==1){
      truc[gamma]=((1-exp((Tc[gamma]-t[gamma])*4*beta*beta/(Xi[gamma]*(2-sigma))))/(4*beta*beta/(Xi[gamma]*(2-sigma))))
    }
  }
  Q[n,1:(n-1)]=(1/t[(1:(n-1))])*exp_truc*truc
  diag(Q)=rep(1,n)-apply(Q,2,sum)
  
  if(cut_edge){
    Q=Q[-n,-n]
    corrector_Q=colSums(Q)
    Q=diag(1/corrector_Q)%*%Q
  }
  
  return(Q)
}

# build all the required matrices for the  HMM
build_HMM_matrix_t<-function(n=40,rho,L,Pop=T,Xi=NA,scale=c(1,0),beta=1,sigma=0,Beta=1,Sigma=0,FS=F,Big_Window=F,tmax=15,alpha_t=0.1,npair=2){
  if(!Big_Window){
    Vect=0:(n-1)
    Tc= scale[2] -(0.5*(2-Sigma)*log(1-(Vect/n))/((Beta^2)*scale[1]))
  }
  if(as.numeric(Big_Window)==1){
    Vect=1:(n-1)
    alpha_t=alpha_t/npair
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])) )
    
  }
  if(as.numeric(Big_Window)==2){
    tmax=100
    Vect=1:(n-1)
    npair=10
    alpha_t=alpha_t/npair
    Tc= c(0,scale[2] + (0.5*(2-Sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((Beta^2)*scale[1])) )
  }
  t=vector()
  q=vector()
  p_gamma=vector()
  Pi_gamma=vector()
  pi=vector()
  if(length(beta)==1){
    beta=rep(beta,n)
  }
  if(length(sigma)==1){
    sigma=rep(sigma,n)
  }
  
  r=rho/(2*(L-1))
  if(length(r)==1){
    r=rep(r,n)
  }
  D=Tc[2:n]-Tc[1:(n-1)]
  if(Pop){
    Xi=rep(1,n)
  }
  t[1:(n-1)]=(Xi[1:(n-1)]*(2-sigma[1:(n-1)])/(beta[1:(n-1)]*beta[1:(n-1)]*2))+ (Tc[1:(n-1)]-(Tc[2:n]*(exp(-beta[1:(n-1)]*beta[1:(n-1)]*(D[1:(n-1)]/Xi[1:(n-1)])*2/(2-sigma[1:(n-1)])))))/(1-exp(-beta[1:(n-1)]*beta[1:(n-1)]*((D[1:(n-1)])/Xi[1:(n-1)])*2/(2-sigma[1:(n-1)])))
  t[n]=(Xi[n]*(2-sigma[n])/(beta[n]*beta[n]*2)) +Tc[n]
  q[1]=(1-exp(-D[1]*2*beta[1]*beta[1]/(Xi[1]*(2-sigma[1]))))
  for(alpha in 2:(n-1)){q[alpha]=exp(-sum((2*beta[1:(alpha-1)]*beta[1:(alpha-1)]/(2-sigma[1:(alpha-1)]))*D[1:(alpha-1)]/Xi[1:(alpha-1)]))*(1-exp(-D[alpha]*2*beta[alpha]*beta[alpha]/(Xi[alpha]*(2-sigma[alpha]))))}
  q[n]=exp(-sum((2*beta[1:(n-1)]*beta[1:(n-1)]/(2-sigma[1:(n-1)]))*D[1:(n-1)]/Xi[1:(n-1)]))
  
  for(ii in 1:n){
    if(ii>1){
      p_gamma[ii]=(1-exp(-2*(sum(D[1:(ii-1)]*r[1:(ii-1)]*beta[1:(ii-1)]*2*(1-sigma[1:(ii-1)])/(2-sigma[1:(ii-1)]))+ ( (t[ii]-Tc[ii])*r[ii]*beta[ii]*2*(1-sigma[ii])/(2-sigma[ii]) ) )))
      Pi_gamma[ii]=(sum(D[1:(ii-1)]*r[1:(ii-1)]*beta[1:(ii-1)]*2*(1-sigma[1:(ii-1)])/(2-sigma[1:(ii-1)]))+( (t[ii]-Tc[ii])*r[ii]*beta[ii]*2*(1-sigma[ii])/(2-sigma[ii]) ))
    }
    if(ii==1){
      p_gamma[ii]=(1-exp(-2*(( (t[ii]-Tc[ii])*r[ii]*beta[ii]*2*(1-sigma[ii])/(2-sigma[ii]) ) )))
      Pi_gamma[ii]=( (t[ii]-Tc[ii])*r[ii]*beta[ii]*2*(1-sigma[ii])/(2-sigma[ii]) )
    }
    pi[ii]=r[ii]*beta[ii]*2*(1-sigma[ii])/(2-sigma[ii])
    
  }
  Q=matrix(0,ncol = n, nrow = n)
  
  
  
  
  for(i in 1:(n-1)){
    if(i==1){
      Q[1,(2:n)]=p_gamma[(2:n)]*(pi[1]/(2*Pi_gamma[(2:n)]))*( D[1]-((1-exp(-D[1]*beta[1]*beta[1]*4/(Xi[1]*(2-sigma[1]))))/(4*beta[1]*beta[1]/(Xi[1]*(2-sigma[1])))))
    }
    if(i==2){
      eta=1
      truc= (pi[eta]*(1-exp(-D[eta]*beta[eta]*beta[eta]*4/(Xi[eta]*(2-sigma[eta]))))/((4*beta[eta]*beta[eta])/(Xi[eta]*(2-sigma[eta]))))
      truc=truc*(1-exp(-D[i]*beta[i]*beta[i]*4/(Xi[i]*(2-sigma[i]))))
      Q[i,((i+1):n)]=p_gamma[((i+1):n)]*(1/(2*Pi_gamma[((i+1):n)]))*( truc + pi[i]*(D[i]-((1-exp(-D[i]*beta[i]*beta[i]*4/(Xi[i]*(2-sigma[i]))))/(4*beta[i]*beta[i]/(Xi[i]*(2-sigma[i]))))))
      Q[i,1]=p_gamma[1]*(1/Pi_gamma[1])*exp((t[1]-Tc[2])*beta[1]*beta[1]*2/(Xi[1]*(2-sigma[1])))*(1-exp(-D[2]*beta[2]*beta[2]*2/(Xi[2]*(2-sigma[2]))))*(( pi[1]*(1-exp(-t[1]*4*beta[1]*beta[1]/(Xi[1]*(2-sigma[1]))))/(4*beta[1]*beta[1]/(Xi[1]*(2-sigma[1])))))
    }
    if( i>2){
      truc=0
      for(eta in 1:(i-2)){
        truc= truc+ ((( exp(-4*sum(beta[(eta+1):(i-1)]*beta[(eta+1):(i-1)]*D[(eta+1):(i-1)]/(Xi[(eta+1):(i-1)]*(2-sigma[(eta+1):(i-1)]))))*pi[eta]*(1-exp(-D[eta]*beta[eta]*beta[eta]*4/(Xi[eta]*(2-sigma[eta]))))/((4*beta[eta]*beta[eta])/(Xi[eta]*(2-sigma[eta]))))))
      }
      eta=i-1
      truc= truc+ (((pi[eta]*(1-exp(-D[eta]*beta[eta]*beta[eta]*4/(Xi[eta]*(2-sigma[eta]))))/((4*beta[eta]*beta[eta])/(Xi[eta]*(2-sigma[eta]))))))
      truc=truc*(1-exp(-D[i]*beta[i]*beta[i]*4/(Xi[i]*(2-sigma[i]))))
      Q[i,((i+1):n)]=p_gamma[((i+1):n)]*(1/(2*Pi_gamma[((i+1):n)]))*( truc + (pi[i]*(D[i]-((1-exp(-D[i]*beta[i]*beta[i]*4/(Xi[i]*(2-sigma[i]))))/(4*beta[i]*beta[i]/(Xi[i]*(2-sigma[i])))))))
      truc=rep(0,length(1:(i-1)))
      exp_truc=rep(0,length(1:(i-1)))
      for(gamma in 1:(i-1)){
        if(gamma<(i-1)){
          exp_truc[gamma]=exp(-2*(sum(beta[(gamma+1):(i-1)]*beta[(gamma+1):(i-1)]*D[(gamma+1):(i-1)]/(Xi[(gamma+1):(i-1)]*(2-sigma[(gamma+1):(i-1)])))+((Tc[gamma+1]-t[gamma])*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma])))))
        }
        if(gamma==(i-1)){
          exp_truc[gamma]=exp(-beta[gamma]*beta[gamma]*2*(((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma[gamma])))))
        }
        if(gamma>2){
          sub_truc=0
          for(eta in 1:(gamma-2)){
            sub_truc=sub_truc+((pi[eta]*(1-exp(-D[eta]*4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))/(4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))*exp(-4*(sum(beta[(eta+1):(gamma-1)]*beta[(eta+1):(gamma-1)]*D[(eta+1):(gamma-1)]/(Xi[(eta+1):(gamma-1)]*(2-sigma[(eta+1):(gamma-1)])))+((t[gamma]-Tc[gamma])*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))))
          }
          eta=gamma-1
          sub_truc=sub_truc+((pi[eta]*(1-exp(-D[eta]*4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))/(4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))*exp(-4*beta[gamma]*beta[gamma]*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma[gamma]))))))
          truc[gamma]=sub_truc+( pi[gamma]*(1-exp((Tc[gamma]-t[gamma])*4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))/(4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))
        }
        if(gamma==2){
          sub_truc=0
          eta=1
          sub_truc=sub_truc+((pi[eta]*(1-exp(-D[eta]*4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))/(4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))*exp(-4*beta[gamma]*beta[gamma]*((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma[gamma])))))
          truc[gamma]=sub_truc+( pi[gamma]*(1-exp((Tc[gamma]-t[gamma])*4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))/(4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))
        }
        if(gamma==1){
          truc[gamma]=(pi[gamma]*(1-exp((Tc[gamma]-t[gamma])*4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))/(4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))
        }
      }
      Q[i,1:(i-1)]=p_gamma[(1:(i-1))]*(1/Pi_gamma[(1:(i-1))])*(1-exp(-D[i]*beta[i]*beta[i]*2/(Xi[i]*(2-sigma[i]))))*exp_truc*truc
    }
  }
  truc=rep(0,length(1:(n-1)))
  exp_truc=rep(0,length(1:(n-1)))
  for(gamma in 1:(n-1)){
    if(gamma<(n-1)){
      exp_truc[gamma]=exp(-2*( sum(beta[(gamma+1):(n-1)]*beta[(gamma+1):(n-1)]*D[(gamma+1):(n-1)]/(Xi[(gamma+1):(n-1)]*(2-sigma[(gamma+1):(n-1)])) ) + (beta[gamma]*beta[gamma]*(Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma[gamma])))))
    }
    if(gamma==(n-1)){
      exp_truc[gamma]=exp(-beta[gamma]*beta[gamma]*2*(((Tc[gamma+1]-t[gamma])/(Xi[gamma]*(2-sigma[gamma])))))
    }
    if(gamma>2){
      sub_truc=0
      for(eta in 1:(gamma-2)){
        sub_truc=sub_truc+((pi[eta]*(1-exp(-D[eta]*4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))/(4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))*exp(-4*(sum(beta[(eta+1):(gamma-1)]*beta[(eta+1):(gamma-1)]*D[(eta+1):(gamma-1)]/(Xi[(eta+1):(gamma-1)]*(2-sigma[(eta+1):(gamma-1)]))) + (beta[gamma]*beta[gamma]*(t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma[gamma]))))))
      }
      eta=gamma-1
      sub_truc=sub_truc+((pi[eta]*(1-exp(-D[eta]*4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))/(4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))*exp(-4*beta[gamma]*beta[gamma]*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma[gamma]))))))
      truc[gamma]=sub_truc+(pi[gamma]*(1-exp((Tc[gamma]-t[gamma])*4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))/(4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))
    }
    if(gamma==2){
      sub_truc=0
      eta=1
      sub_truc=sub_truc+(( pi[eta]*(1-exp(-D[eta]*4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))/(4*beta[eta]*beta[eta]/(Xi[eta]*(2-sigma[eta]))))*exp(-4*beta[eta]*beta[eta]*(((t[gamma]-Tc[gamma])/(Xi[gamma]*(2-sigma[gamma]))))))
      truc[gamma]=sub_truc+( pi[gamma]*(1-exp((Tc[gamma]-t[gamma])*4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))/(4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))
    }
    if(gamma==1){
      truc[gamma]=(pi[gamma]*(1-exp((Tc[gamma]-t[gamma])*4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))/(4*beta[gamma]*beta[gamma]/(Xi[gamma]*(2-sigma[gamma]))))
    }
  }
  Q[n,1:(n-1)]=p_gamma[(1:(n-1))]*(1/Pi_gamma[(1:(n-1))])*exp_truc*truc
  diag(Q)=rep(1,n)-apply(Q,2,sum)
  output=list()
  output[[1]]=Q
  output[[2]]=q
  output[[3]]=t
  output[[4]]=Tc
  return(output)
}

#######################
# Baum-Welch functions#
#######################

sum_chi_corrected<-function(fo,ba,Q,ob,E,W,c){
  if(as.numeric(ob)<=1){
    truc=(fo%*%(t(ba*E[,(as.numeric(ob)+1)]/c)))
  }
  if(as.numeric(ob)>1){
    U=t(W$P)%*%fo%*%t(ba/c)%*%t(W$P_)
    Q=Q[[(as.numeric(ob)+1)]]
    truc=(t(W$P_)%*%(U*Q)%*%(t(W$P))%*%diag(E[,1]))
  }
  return(truc)
}

sum_M_cor<-function(fo,ba,ob,Q,W,c,TO){
  M=matrix(0,nrow=length(fo),ncol=2)
  if(as.numeric(ob)<=1){
    truc=fo*(TO[[(as.numeric(ob)+1)]]%*%(ba/c))
    M[,(as.numeric(ob)+1)]=truc
  }
  if(as.numeric(ob)>1){
    U=(t(W$P)%*%fo%*%t(ba/c)%*%t(W$P_))
    Q=Q[[(as.numeric(ob)+1)]]
    A=(U*Q)
    truc=diag(t(W$P_)%*%A%*%t(W$P))
    M[,1]=truc
  }
  output=M
  return(output)
}

Baum_Welch_algo<-function(Os, maxIt =20,L,mu,theta_W,Rho,beta=1,Popfix=T,SB=F,k=20,BoxB=c(0.1,1),BoxP=c(3,3),Boxr=c(1,1),maxBit=1,pop_vect=NA,window_scaling=c(1,0),sigma=0.00,SF=F,Boxs=c(0,0.97),ER=F,BW=F,NC=1,redo_R=F,mu_b=1,FS=F,SCALED=F,Big_Window=F,cut_edge=F,Share_r=T,mu_r=NA){
  maxIt_o=maxIt
  Xi=NA
  print(paste("sequence length :",L,sep=" "))
  if(length(Rho)>1&length(Rho)!=NC){
    stop("Problem in recombination definition")
  }
  if(FS){
    Ne=mean( (log((1+((log((1-(theta_W/(L*0.75))))/2))))/(log(1-(mu_r*4/3)))))
  }else{
    Ne=NA
  }
  old_list=list()
  n <- length(Os[[1]][[1]]);
  Pop=Popfix
  theta=mu*2*L
  if(FS){
    theta=(0.75-0.75*exp(-mu*2))*L
  }
  gamma=Rho/theta
  if(Share_r){
    print("former Rho:")
    print(Rho)
    print("old mu :")
    print(mu)
    mu=mean(mu)
    print("new mu :")
    print(mu)
    theta=mu*2*L
    Rho=theta*gamma
    print("new Rho:")
    print(Rho)
  }
  gamma_o=gamma
  print(gamma)
  test.env <- new.env()
  test.env$L <- L
  test.env$k <- k
  test.env$mu <- mu
  test.env$Rho <- Rho
  test.env$window_scaling <- window_scaling
  test.env$BW<-BW
  test.env$Pop<-Popfix
  test.env$NC<-NC
  test.env$FS<-FS
  test.env$mu_b <- mu_b
  test.env$Big_Window <- Big_Window
  test.env$cut_edge <- cut_edge
  test.env$Share_r<-Share_r
  if(NC>1){
    npair=NC
    test.env$npair <- NC
  }
  if(NC==1){
    npair=length(Os)
    test.env$npair <- length(Os)
  }
  mb=0
  if(SB){
    BoxB[1]=max(sqrt(0.01),sqrt(BoxB[1]))
    BoxB[2]=min(sqrt(1),sqrt(BoxB[2]))
    beta=max((BoxB[1]^2),beta)
    beta=min(beta,(BoxB[2]^2))
    Beta=beta
  }
  if(SF){
    sigma=min(Boxs[2],sigma)
    sigma=max(sigma,Boxs[1])
    Self=sigma
  }
  if(any(!is.na(pop_vect))){
    Klink=length(pop_vect)
  }
  if((all(is.na(pop_vect))|sum(pop_vect)!=k)){
    Klink=0.5*k
    pop_vect=rep(2, Klink)
    print("Default pop vector")
  }
  test.env$pop_vect <- pop_vect
  if(!SB&!SF){
    maxBit=1
  }
  if(!SB){
    Beta=beta
    test.env$beta <- beta
  }
  if(!SF){
    Self=sigma
    oldSelf=Self
    test.env$sigma <- sigma
  }
  if(!ER){
    Boxr=c(0,0)
  }
  while(mb<maxBit){
    maxIt=maxIt_o
    print(paste("Beta:",Beta))
    print(paste("Self:",Self))
    test.env$Beta <- Beta
    test.env$Self <- Self
    test.env$BoxB <- BoxB
    test.env$Boxs <- Boxs
    mb=mb+1
    diff=1
    it <- 0
    if(mb==1){
      if(SB){
        oldbeta=(sqrt(beta)-BoxB[1])/(BoxB[2]-BoxB[1])
      }
      if(SF){
        oldsigma=(sigma-Boxs[1])/(Boxs[2]-Boxs[1])
      }
      
    }
    if(mb>1){
      if(SB){
        beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1])^2
      }
      if(SF){
        sigma=oldsigma*(Boxs[2]-Boxs[1])
        sigma=sigma+Boxs[1]
      }
      if(!Popfix){
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          oldXi[x:xx]=oldXi_[ix]
        }
        Xi_=oldXi*sum(BoxP)
        Xi_=Xi_-(BoxP[1])
        Xi_=10^Xi_
      }
      
      if(NC==1){
        theta=(((theta_W*(beta^2)))*2/((2-sigma)*(beta+((1-beta)*mu_b))))
        mu=theta/(2*L)
        Rho=gamma*theta
      }
      if(NC>1){
        
        if(Share_r){
          theta=((theta_W*(beta^2)))*2/(2-sigma)
          mu=mean(theta/(2*L))
          theta=mu*2*L
          Rho=gamma*theta
        }else{
          theta=((theta_W*(beta^2)))*2/(2-sigma)
          mu=theta/(2*L)
          Rho=gamma*theta
        }
        
        
      }
      
      test.env$mu <- mu
      test.env$Rho <- Rho
    }
    if(!ER){
      if(NC==1|Share_r){
        oldrho=0
      }else{
        oldrho=rep(0,NC)
      }
      
    }
    if(ER){
      if(NC==1|Share_r){
        oldrho=(Boxr[1]/sum(Boxr))
      }else{
        oldrho=rep((Boxr[1]/sum(Boxr)),NC)
      }
    }
    print(length(oldrho))
    oldXi_=rep((BoxP[1]/sum(BoxP)),Klink)
    oldXi=vector()
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      oldXi[x:xx]=oldXi_[ix]
    }
    Do_BW=T
    restart=F
    diff_conv=vector()
    while (it<maxIt&!restart){
      start_time <- Sys.time()
      restart=F
      if(!Do_BW){
        it=0
        oldXi_=rep((BoxP[1]/sum(BoxP)),Klink)
        oldXi=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          oldXi[x:xx]=oldXi_[ix]
        }
      }
      it <- it+1;
      print(paste("It:",it))
      if(Popfix){
        rho_=oldrho*sum(Boxr)
        rho_=rho_-(Boxr[1])
        rho_=10^(rho_)
        rho_=rho_*Rho
        if(SB){
          beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1])^2
        }
        if(SF){
          sigma=oldsigma*(Boxs[2]-Boxs[1])
          sigma=sigma+Boxs[1]
        }
        print(c("sigma:",sigma,"beta :",beta))
        if(Share_r){
          
          print(c("rho/theta:",mean(rho_/theta)))
        }else{
          print(c("rho/theta:",rho_/theta))
        }
        Keep_going=F
        if(it==1){
          diff_o=0
        }
        if(it>1){
          
          
          gamma_temp=mean(rho_/theta)
          diff_o=max(abs(sigma_o-sigma),abs(beta_o-beta),abs(gamma_temp_o-gamma_temp))
          count_diff_o=0
          if(diff_o>=0.005){
            if(it==maxIt){
              count_diff_o=count_diff_o+1
              maxIt=maxIt+1
            }
          }
        }
        sigma_o=sigma
        beta_o=beta
        gamma_temp_o=mean(rho_/theta)
        if(NC==1){
          builder=build_HMM_matrix(k,(rho_),beta,L=L,Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
        }
        if(NC>1){
          if(Share_r){
            builder=build_HMM_matrix(k,(rho_[1]),beta,L=L[1],Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
          }else{
            builder=list()
            for(chr in 1:NC){
              builder[[chr]]=build_HMM_matrix(k,(rho_[chr]),beta,L=L[chr],Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
            }
          }
        }
        
        
      }
      if(!Popfix){
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          oldXi[x:xx]=oldXi_[ix]
        }
        Xi_=oldXi*sum(BoxP)
        Xi_=Xi_-(BoxP[1])
        Xi_=10^Xi_
        rho_=oldrho*sum(Boxr)
        rho_=rho_-(Boxr[1])
        rho_=10^(rho_)
        rho_=rho_*Rho
        if(SB){
          beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1])^2
        }
        if(SF){
          sigma=oldsigma*(Boxs[2]-Boxs[1])
          sigma=sigma+Boxs[1]
        }
        print(c("sigma:",sigma,"beta :",beta))
        if(Share_r){
          
          print(c("rho/theta:",mean(rho_/theta)))
        }else{
          print(c("rho/theta:",rho_/theta))
        }
        Keep_going=F
        if(it==1){
          diff_o=0
        }
        if(it>1){
          gamma_temp=mean(rho_/theta)
          diff_o=max(abs(sigma_o-sigma),abs(beta_o-beta),abs(gamma_temp_o-gamma_temp))
          count_diff_o=0
          if(diff_o>=0.005){
            if(it==maxIt){
              count_diff_o=count_diff_o+1
              maxIt=maxIt+1
            }
          }
        }
        sigma_o=sigma
        beta_o=beta
        gamma_temp_o=mean(rho_/theta)
        if(NC==1){
          builder=build_HMM_matrix(k,(rho_),beta,L=L,Pop=Pop,Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
        }
        if(NC>1){
          if(Share_r){
            builder=build_HMM_matrix(k,(rho_[1]),beta,L=L[1],Pop=Pop,Xi=Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
          }else{
            builder=list()
            for(chr in 1:NC){
              builder[[chr]]=build_HMM_matrix(k,(rho_[chr]),beta,L=L[chr],Pop=Pop,Xi=Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
            }
          }
        }
        
        
      }
      if(NC==1){
        Q = builder[[1]]
        nu= builder[[2]]
        Tc=builder[[3]]
        g=matrix(0,nrow=length(Tc),ncol=2)
        if(!FS){
          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
        }
        if(FS){
          a=0.8125
          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
        }
        M=matrix(0,nrow=length(Tc),ncol=2)
        N=matrix(0,length(Tc),length(Tc))
        M=matrix(0,nrow=length(Tc),ncol=2)
        N=matrix(0,length(Tc),length(Tc))
        MLH=0
        q_=rep(0,length(Tc))
        s_t=Sys.time()
        test=Build_zip_Matrix_mailund(Q,g,Os[[1]][[2]],nu)
        e_t=Sys.time()
        print(e_t-s_t)
        
        for(i in 1:length(Os)){
          fo=forward_zip_mailund(Os[[i]][[1]],g,nu,test[[1]])
          MLH=MLH+fo[[3]]
          c=exp(fo[[2]])
          ba=Backward_zip_mailund(Os[[i]][[1]],test[[3]],length(Tc),c)
          
          W=list()
          int=t(Q)%*%diag(g[,1])
          int=eigen(int)
          W$P=int$vectors
          W$P_=solve(W$P)
          
          symbol= Os[[i]][[3]][,1]
          for(ob in 1){
            truc_M=matrix(0,nrow=length(Tc),ncol=2)
            if(Os[[i]][[1]][(ob+1)]<=1){
              truc_N=(fo[[1]][,ob]%*%(t(ba[,(ob+1)]*g[,(as.numeric(Os[[i]][[1]][(ob+1)])+1)]/c[(ob+1)])))
              truc_M[,(as.numeric(Os[[i]][[1]][(ob+1)])+1)]=fo[[1]][,ob]*(test[[2]][[(as.numeric(Os[[i]][[1]][(ob+1)])+1)]]%*%(ba[,(ob+1)]/c[(ob+1)]))
            }
            if(Os[[i]][[1]][(ob+1)]>1){
              truc_N=(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_)*test[[4]][[(as.numeric(Os[[i]][[1]][(ob+1)])+1)]])%*%(t(W$P))%*%diag(g[,1]))
              truc_M[,1]=diag(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_))*test[[2]][[(as.numeric(Os[[i]][[1]][(ob+1)])+1)]]%*%t(W$P))
            }
            #truc_N=sum_chi_corrected(fo[[1]][,ob],ba[,(ob+1)],test[[4]],Os[[i]][[1]][(ob+1)],g,W,c[(ob+1)])
            N=N+truc_N
            #truc_M=sum_M_cor(fo[[1]][,ob],ba[,(ob+1)],Os[[i]][[1]][(ob+1)],test[[2]],W,(c[(ob+1)]),test[[1]])
            M=M+truc_M
          }
          for(sym in symbol){
            ob=as.numeric(sym)
            pos=which(as.numeric(Os[[i]][[1]][-c(1,length(Os[[i]][[1]]))])==ob)
            pos=pos+1
            if(ob<2){
              ba_t=t(t(ba[,pos])/c[(pos)])
              truc=c(rowSums(fo[[1]][,(pos-1)]*(test[[1]][[(ob+1)]]%*%ba_t)))
              truc=(truc/(sum(truc)))*length(pos)
              M[,(ob+1)]=M[,(ob+1)]+ truc
              N=N+(fo[[1]][,(pos-1)]%*%(t(diag(g[,(ob+1)])%*%ba_t)))
            }
            if(ob>=2){
              A=(t(W$P)%*%fo[[1]][,(pos-1)]%*%t(t(t(ba[,pos])/c[(pos)]))%*%t(W$P_))
              A_=A*test[[2]][[(ob+1)]]
              A_=(t(W$P_)%*%A_%*%t(W$P))
              M[,1]=M[,1]+(diag(A_))
              A_=A*test[[4]][[(ob+1)]]
              A_=(t(W$P_)%*%A_%*%t(W$P))
              N=N+((A_%*%diag(g[,1])))
            }
            
          }
          if(as.numeric(Os[[i]][[1]][length(Os[[i]][[1]])])==1){
            M[,2]=M[,2]+(fo[[1]][,length(Os[[i]][[1]])]*ba[,length(Os[[i]][[1]])])
          }
          if(as.numeric(Os[[i]][[1]][length(Os[[i]][[1]])])!=1){
            M[,1]=M[,1]+((fo[[1]][,length(Os[[i]][[1]])]*ba[,length(Os[[i]][[1]])])/sum(fo[[1]][,length(Os[[i]][[1]])]*ba[,length(Os[[i]][[1]])]))
          }
          q_=q_+((fo[[1]][,1]*ba[,1])/sum(fo[[1]][,1]*ba[,1]))
        }
        
        N=N*t(Q)
        if(is.complex(N)){
          N=Re(N)
        }
        
        if(is.complex(M)){
          M=Re(M)
        }
        Scale_N=(L-1)/sum(N)
        
        N=N*Scale_N
        q_=q_/sum(q_)
        Scale_M=L/sum(M)
        M=M*Scale_M
        
        if(SCALED){
          corrector_N=rowSums(N)
          N=diag(1/corrector_N)%*%N
          corrector_M=rowSums(M)
          M=diag(1/corrector_M)%*%M
        }
        
        
      }
      if(NC>1){
        
        if(Share_r){
          Q = builder[[1]]
          nu= builder[[2]]
          Tc=builder[[3]]
          g=matrix(0,nrow=length(Tc),ncol=2)
          if(!FS){
            g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
            g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
          }
          if(FS){
            a=0.8125
            g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
            g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
          }
          
          M=matrix(0,nrow=length(Tc),ncol=2)
          N=matrix(0,length(Tc),length(Tc))
          MLH=0
          q_=rep(0,length(Tc))
          s_t=Sys.time()
          test=Build_zip_Matrix_mailund(Q,g,Os[[1]][[1]][[2]],nu)
          e_t=Sys.time()
          print(e_t-s_t)
          for(chr in 1:NC ){
            for(i in 1:length(Os[[chr]])){
              fo=forward_zip_mailund(Os[[chr]][[i]][[1]],g,nu,test[[1]])
              MLH=MLH+fo[[3]]
              c=exp(fo[[2]])
              ba=Backward_zip_mailund(Os[[chr]][[i]][[1]],test[[3]],length(Tc),c)
              W=list()
              int=t(Q)%*%diag(g[,1])
              int=eigen(int)
              W$P=int$vectors
              W$P_=solve(W$P)
              
              symbol= Os[[chr]][[i]][[3]][,1]
              for(ob in 1){
                truc_M=matrix(0,nrow=length(Tc),ncol=2)
                if(Os[[chr]][[i]][[1]][(ob+1)]<=1){
                  truc_N=(fo[[1]][,ob]%*%(t(ba[,(ob+1)]*g[,(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]/c[(ob+1)])))
                  truc_M[,(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]=fo[[1]][,ob]*(test[[2]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]]%*%(ba[,(ob+1)]/c[(ob+1)]))
                }
                if(Os[[chr]][[i]][[1]][(ob+1)]>1){
                  truc_N=(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_)*test[[4]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]])%*%(t(W$P))%*%diag(g[,1]))
                  truc_M[,1]=diag(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_))*test[[2]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]]%*%t(W$P))
                }
                N=N+truc_N
                M=M+truc_M
              }
              for(sym in symbol){
                ob=as.numeric(sym)
                pos=which(as.numeric(Os[[chr]][[i]][[1]][-c(1,length(Os[[chr]][[i]][[1]]))])==ob)
                pos=pos+1
                if(ob<2){
                  ba_t=t(t(ba[,(pos)])/c[(pos)])
                  truc=c(rowSums(fo[[1]][,(pos-1)]*(test[[1]][[(ob+1)]]%*%ba_t)))
                  truc=(truc/(sum(truc)))*length(pos)
                  M[,(ob+1)]=M[,(ob+1)]+ truc
                  N=N+(fo[[1]][,(pos-1)]%*%(t(diag(g[,(ob+1)])%*%ba_t)))
                }
                if(ob>=2){
                  A=(t(W$P)%*%fo[[1]][,(pos-1)]%*%t(t(t(ba[,(pos)])/c[(pos)]))%*%t(W$P_))
                  A_=A*test[[2]][[(ob+1)]]
                  A_=(t(W$P_)%*%A_%*%t(W$P))
                  M[,1]=M[,1]+(diag(A_))
                  A_=A*test[[4]][[(ob+1)]]
                  A_=(t(W$P_)%*%A_%*%t(W$P))
                  N=N+((A_%*%diag(g[,1])))
                }
                
              }
              if(as.numeric(Os[[chr]][[i]][[1]][length(Os[[chr]][[i]][[1]])])==1){
                M[,2]=M[,2]+(fo[[1]][,length(Os[[chr]][[i]][[1]])]*ba[,length(Os[[chr]][[i]][[1]])])
              }
              if(as.numeric(Os[[chr]][[i]][[1]][length(Os[[chr]][[i]][[1]])])!=1){
                M[,1]=M[,1]+((fo[[1]][,length(Os[[chr]][[i]][[1]])]*ba[,length(Os[[chr]][[i]][[1]])])/sum(fo[[1]][,length(Os[[chr]][[i]][[1]])]*ba[,length(Os[[chr]][[i]][[1]])]))
              }
              q_=q_+((fo[[1]][,1]*ba[,(1)])/sum(fo[[1]][,1]*ba[,(1)]))
            }
            
          }
          
          N=N*t(Q)
          Scale_N=(sum(L)-length(L))/sum(N)
          N=N*Scale_N
          Scale_M=sum(L)/sum(M)
          M=M*Scale_M
          q_=q_/sum(q_)
          
          if(SCALED){
            corrector_N=rowSums(N)
            N=diag(1/corrector_N)%*%N
            corrector_M=rowSums(M)
            M=diag(1/corrector_M)%*%M
          }
          
        }else{
          
          Q=list()
          nu=list()
          Tc=list()
          g=list()
          M=list()
          N=list()
          MLH=list()
          q_=list()
          for(chr in 1:NC){
            Q[[chr]] = builder[[chr]][[1]]
            nu[[chr]]= builder[[chr]][[2]]
            Tc[[chr]]=builder[[chr]][[3]]
            g[[chr]]=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
            if(!FS){
              g[[chr]][,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]])
              g[[chr]][,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]])
            }
            if(FS){
              
              a=0.8125
              g[[chr]][,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]])))
              g[[chr]][,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]])))
            }
            
            M[[chr]]=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
            N[[chr]]=matrix(0,length(Tc[[chr]]),length(Tc[[chr]]))
            MLH[[chr]]=0
            q_[[chr]]=rep(0,length(Tc[[chr]]))
            test=Build_zip_Matrix_mailund(Q[[chr]],g[[chr]],Os[[chr]][[1]][[2]],nu[[chr]])
            
            for(i in 1:length(Os[[chr]])){
              fo=forward_zip_mailund(Os[[chr]][[i]][[1]],g[[chr]],nu[[chr]],test[[1]])
              MLH[[chr]]=MLH[[chr]]+fo[[3]]
              c=exp(fo[[2]])
              ba=Backward_zip_mailund(Os[[chr]][[i]][[1]],test[[3]],length(Tc[[chr]]),c)
              W=list()
              int=t(Q[[chr]])%*%diag(g[[chr]][,1])
              int=eigen(int)
              W$P=int$vectors
              W$P_=solve(W$P)
              for(ob in 1){
                truc_M=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
                if(Os[[chr]][[i]][[1]][(ob+1)]<=1){
                  truc_N=(fo[[1]][,ob]%*%(t(ba[,(ob+1)]*g[[chr]][,(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]/c[(ob+1)])))
                  truc_M[,(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]=fo[[1]][,ob]*(test[[2]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]]%*%(ba[,(ob+1)]/c[(ob+1)]))
                }
                if(Os[[chr]][[i]][[1]][(ob+1)]>1){
                  truc_N=(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_)*test[[4]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]])%*%(t(W$P))%*%diag(g[[chr]][,1]))
                  truc_M[,1]=diag(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_))*test[[2]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]]%*%t(W$P))
                }
                N[[chr]]=N[[chr]]+truc_N
                M[[chr]]=M[[chr]]+truc_M
              }
              symbol= Os[[chr]][[i]][[3]][,1]
              for(sym in symbol){
                ob=as.numeric(sym)
                pos=which(as.numeric(Os[[chr]][[i]][[1]][-c(1,length(Os[[chr]][[i]][[1]]))])==ob)
                pos=pos+1
                if(ob<2){
                  ba_t=t(t(ba[,pos])/c[(pos)])
                  truc=c(rowSums(fo[[1]][,(pos-1)]*(test[[1]][[(ob+1)]]%*%ba_t)))
                  truc=(truc/(sum(truc)))*length(pos)
                  M[[chr]][,(ob+1)]=M[[chr]][,(ob+1)]+ truc
                  N[[chr]]=N[[chr]]+(fo[[1]][,(pos-1)]%*%(t(diag(g[[chr]][,(ob+1)])%*%ba_t)))
                }
                if(ob>=2){
                  A=(t(W$P)%*%fo[[1]][,(pos-1)]%*%t(t(t(ba[,pos])/c[(pos)]))%*%t(W$P_))
                  A_=A*test[[2]][[(ob+1)]]
                  A_=(t(W$P_)%*%A_%*%t(W$P))
                  M[[chr]][,1]=M[[chr]][,1]+(diag(A_))
                  A_=A*test[[4]][[(ob+1)]]
                  A_=(t(W$P_)%*%A_%*%t(W$P))
                  N[[chr]]=N[[chr]]+((A_%*%diag(g[[chr]][,1])))
                }
              }
              if(as.numeric(Os[[chr]][[i]][[1]][length(Os[[chr]][[i]][[1]])])==1){
                M[[chr]][,2]=M[[chr]][,2]+(fo[[1]][,length(Os[[chr]][[i]][[1]])]*(ba[,length(Os[[chr]][[i]][[1]])]/(c[length(Os[[chr]][[i]][[1]])])))
              }
              if(as.numeric(Os[[chr]][[i]][[1]][length(Os[[chr]][[i]][[1]])])!=1){
                M[[chr]][,1]=M[[chr]][,1]+(fo[[1]][,length(Os[[chr]][[i]][[1]])]*(ba[,length(Os[[chr]][[i]][[1]])]/(c[length(Os[[chr]][[i]][[1]])])))
              }
              q_[[chr]]=q_[[chr]]+((fo[[1]][,1]*ba[,1])/sum(fo[[1]][,1]*ba[,1]))
              
            }
            
            N[[chr]]=N[[chr]]*t(Q[[chr]])
            if(is.complex(N[[chr]])){
              N[[chr]]=Re(N[[chr]])
            }
            
            if(is.complex(M[[chr]])){
              M[[chr]]=Re(M[[chr]])
            }
            print("N:")
            print(sum(N[[chr]]))
            Scale_N=(L[chr]-1)/sum(N[[chr]])
            N[[chr]]=N[[chr]]*Scale_N
            print("N*:")
            print(sum(N[[chr]]))
            Scale_M=L[chr]/sum(M[[chr]])
            M[[chr]]=M[[chr]]*Scale_M
            q_[[chr]]=q_[[chr]]/sum(q_[[chr]])
            
            if(SCALED){
              corrector_M=rowSums(t(M[[chr]]))
              M[[chr]]=diag(1/corrector_M)%*%M[[chr]]
              corrector_N=rowSums(N[[chr]])
              N[[chr]]=diag(1/corrector_N)%*%N[[chr]]
            }
          }
        }
      }
      if(it>1){
        print(paste(" old Likelihood: ",oldMLH))
      }
      if(NC==1){
        print(paste(" New Likelihood: ",MLH))
        
        if(it>1){
          
          if(oldMLH > MLH|MLH=="NaN"){
            
            if(it>2){
              restart=T
              if(!Popfix){
                oldXi_=oldXi_s
              }
              if(ER){
                oldrho=oldrho_s
              }
              if(SB){
                oldbeta=oldbeta_s
              }
              if(SF){
                oldsigma=oldsigma_s
              }
            }
            if(it==2){
              print("Algortihm has converge but may need more data for better results ")
              restart=T
              MB=maxBit
            }
          }
          
        }
        oldMLH=MLH
      }
      if(NC>1){
        if(Share_r){
          print(paste("New Likelihood: ",MLH))
          if(it>1){
            if(oldMLH > MLH|MLH=="NaN"){
              if(it>1){
                restart=T
                if(!Popfix){
                  oldXi_=oldXi_s
                }
                if(ER){
                  oldrho=oldrho_s
                }
                if(SB){
                  oldbeta=oldbeta_s
                }
                if(SF){
                  oldsigma=oldsigma_s
                }
              }
              if(it==2){
                print("Algortihm has converge but may need more data for better results ")
                restart=T
                MB=maxBit
              }
            }
            
          }
          oldMLH=MLH
        }else{
          MLH1=0
          for(chr in 1:length(MLH)){
            MLH1=MLH[[chr]]+MLH1
          }
          print(paste("New Likelihood: ",MLH1))
          if(it>1){
            if(oldMLH > MLH1|MLH1=="NaN"){
              if(it>1){
                restart=T
                if(!Popfix){
                  oldXi_=oldXi_s
                }
                if(ER){
                  oldrho=oldrho_s
                }
                if(SB){
                  oldbeta=oldbeta_s
                }
                if(SF){
                  oldsigma=oldsigma_s
                }
              }
              if(it==2){
                print("Algortihm has converge but may need more data for better results ")
                restart=T
                MB=maxBit
              }
            }
            
          }
          oldMLH=MLH1
        }
      }
      if(!restart){
        
        if(NC==1){
          x=as.vector(g)
          keep=which(x>0)
          x=x[keep]
          m=as.vector(M)
          m=m[keep]
          A=as.vector(t(Q))
          keep=which(A>0&as.vector(N)>0)
          A=A[keep]
          Big_Xi=as.vector(N)
          Big_Xi=Big_Xi[keep]
          if(BW){
            LH=Re(sum(log(A)*Big_Xi)+sum(log(x)*m)+sum(log(nu)*q_))
          }
          if(!BW){
            LH=Re(sum(log(A)*Big_Xi))
          }
          print(paste(" old Complete likelihood : ", LH ))
          oldLH=-LH
          test.env$Big_Xi <- N
          test.env$Big_M <-M
          test.env$q_ <-q_
        }
        if(NC>1){
          test.env$Big_Xi <- N
          test.env$Big_M <-M
          test.env$q_ <-q_
        }
        Do_BW=T
        if(!Popfix){
          oldXi_s=oldXi_
        }
        if(ER){
          oldrho_s=oldrho
        }
        if(SB){
          oldbeta_s=oldbeta
        }
        if(SF){
          oldsigma_s=oldsigma
        }
        lr=length(oldrho)
        test.env$lr<-lr
        if(Do_BW){
          if(NC==1){
            if(ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[3]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      
                      builder=build_HMM_matrix(n,rho_,beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge) #
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:3,1])
                    rho=sol[1]
                    beta_=sol[2]
                    sigma_=sol[3]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldsigma=sigma_
                    
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[3]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[4:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,rho_,beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(3+(Klink)),1])
                    rho=sol[1]
                    beta_=sol[2]
                    sigma_=sol[3]
                    Xi_=sol[4:length(sol)]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma',envir = test.env)
                      Self=get('Self',envir = test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,rho_,beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:2,1])
                    rho=sol[1]
                    beta_=sol[2]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_)))
                    oldrho=rho
                    oldbeta=beta_
                    
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma',envir = test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BETA=get('BETA',envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[3:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2+(Klink)),1])
                    rho=sol[1]
                    beta_=sol[2]
                    Xi_=sol[3:length(sol)]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,Xi_-oldXi_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldXi_=Xi_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldsigma), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:2,1])
                    rho=sol[1]
                    sigma_=sol[2]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho,oldsigma-sigma_)))
                    oldrho=rho
                    oldsigma=sigma_
                    
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[3:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2+(Klink)),1])
                    rho=sol[1]
                    sigma_=sol[2]
                    Xi_=sol[3:length(sol)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_,oldsigma-sigma_)))
                    oldrho=rho
                    oldXi_=Xi_
                    oldsigma=sigma_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1,1])
                    rho=sol[1]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho)))
                    oldrho=rho
                    
                  }
                  
                  if(!Popfix){
                    
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      rho_=param[1]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[2:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+(Klink)),1])
                    rho=sol[1]
                    Xi_=sol[2:length(sol)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
                    oldrho=rho
                    oldXi_=Xi_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
            }
            if(!ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:2,1])
                    beta_=sol[1]
                    sigma_=sol[2]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_)))
                    oldbeta=beta_
                    oldsigma=sigma_
                    
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[3:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2+(Klink)),1])
                    beta_=sol[1]
                    sigma_=sol[2]
                    Xi_=sol[3:length(sol)]
                    diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1,1])
                    beta_=sol[1]
                    print(paste(" new Complete likelihood : ", LH ))
                    
                    diff=max(abs(c(oldbeta-beta_)))
                    oldbeta=beta_
                    
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[2:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge) #
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+(Klink)),1])
                    beta_=sol[1]
                    Xi_=sol[2:length(sol)]
                    
                    diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1,1])
                    sigma_=sol[1]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(oldsigma-sigma_)))
                    oldsigma=sigma_
                    
                  }
                  
                  if(!Popfix){
                    
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[2:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+(Klink)),1])
                    sigma_=sol[1]
                    Xi_=sol[2:length(sol)]
                    diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
                    oldXi_=Xi_
                    oldsigma=sigma_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[1:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        a=0.8125
                        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((Klink)),1])
                    Xi_=sol[1:length(sol)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
                    oldXi_=Xi_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                
                
              }
            }
          }
          if(NC>1){
            if(ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[lr+2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      LH=0
                      
                      Share_r=get('Share_r', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(length(oldrho)+2),1])
                    rho=sol[1:length(oldrho)]
                    beta_=sol[length(oldrho)+1]
                    sigma_=sol[length(oldrho)+2]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[lr+2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[3+lr:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      cut_edge=get('cut_edge', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2+length(oldrho)+(Klink)),1])
                    rho=sol[1:length(oldrho)]
                    beta_=sol[length(oldrho)+1]
                    sigma_=sol[length(oldrho)+2]
                    Xi_=sol[(length(oldrho)+3):length(sol)]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      LH=0
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)}
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                        
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])}
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(length(oldrho)+1),1])
                    rho=sol[1:length(oldrho)]
                    beta_=sol[length(oldrho)+1]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_)))
                    oldrho=rho
                    oldbeta=beta_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(lr+2):length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for( chr in 1:chr){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+length(oldrho)+(Klink)),1])
                    rho=sol[1:length(oldrho)]
                    beta_=sol[1+length(oldrho)]
                    Xi_=sol[(2+length(oldrho)):length(sol)]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,Xi_-oldXi_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      sigma=param[1+lr]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      LH=0
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){  LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(length(oldrho)+1),1])
                    rho=sol[1:length(oldrho)]
                    sigma_=sol[length(oldrho)+1]
                    diff=max(abs(c(rho- oldrho,oldsigma-sigma_)))
                    oldrho=rho
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(2+lr):length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1+lr]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      LH=0
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){  LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){  LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for( chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+length(oldrho)+(Klink)),1])
                    rho=sol[1:length(oldrho)]
                    sigma_=sol[1+length(oldrho)]
                    Xi_=sol[(2+length(oldrho)):length(sol)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_,oldsigma-sigma_)))
                    oldrho=rho
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      LH=0
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for( chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(length(oldrho)),1])
                    rho=sol[1:length(oldrho)]
                    diff=max(abs(c(rho- oldrho)))
                    oldrho=rho
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      lr=get('lr', envir=test.env)
                      rho_=param[1:lr]
                      rho_=rho_*sum(Boxr)
                      rho_=rho_-(Boxr[1])
                      rho_=10^(rho_)
                      Rho=get('Rho', envir=test.env)
                      rho_=rho_*Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(1+lr):length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      LH=0
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0)#&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        
                        for( chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(length(oldrho)+(Klink)),1])
                    rho=sol[1:length(oldrho)]
                    Xi_=sol[(1+length(oldrho)):length(sol)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
                    oldrho=rho
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
            }
            if(!ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      LH=0
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:2,1])
                    beta_=sol[1]
                    sigma_=sol[2]
                    diff=max(abs(c(oldbeta- beta,oldsigma-sigma)))
                    oldbeta=sol[1]
                    oldsigma=sol[2]
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC',envir = test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[2]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[3:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Big_M=get('Big_M', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2+(Klink)),1])
                    beta_=sol[1]
                    sigma_=sol[2]
                    Xi_=sol[3:length(sol)]
                    diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir = test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1,1])
                    beta_=sol[1]
                    diff=max(abs(c(oldbeta-beta_)))
                    oldbeta=beta_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC',envir = test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[2:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Big_M=get('Big_M', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                        }
                      }
                      
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+(Klink)),1])
                    beta_=sol[1]
                    Xi_=sol[2:length(sol)]
                    diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1,1])
                    sigma_=sol[1]
                    diff=max(abs(c(oldsigma-sigma_)))
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir = test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[2:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      cut_edge=get('cut_edge', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(1+(Klink)),1])
                    sigma_=sol[1]
                    Xi_=sol[2:length(sol)]
                    diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(sigma_)
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC',envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[1:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Share_r=get('Share_r', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Big_M=get('Big_M', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      if(Share_r){
                        builder=build_HMM_matrix(n,(rho_[1]),beta,Pop = Pop,Xi=Xi,L=L[1],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi)>0)
                        A=A[keep]
                        Big_Xi=as.vector(Big_Xi)
                        Big_Xi=Big_Xi[keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          a=0.8125
                          
                          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
                          
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M)
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi)
                        }
                      }else{
                        for(chr in 1:NC){
                          builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            a=0.8125
                            
                            g[,2]= (a - (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            g[,1]= a*((1-a) + (a*exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc)))
                            
                            
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                      }
                      
                      return(LH)
                    }
                    sol= BBoptim(c(oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((Klink)),1])
                    Xi_=sol[1:length(sol)]
                    diff=max(abs(Xi_-oldXi_))
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                    
                  }
                }
              }
            }
          }
        }
      }
      if(diff_o>=0.003){
        diff=(max(diff_o,diff))
      }
      diff_conv=c(diff_conv,diff_o)
      end_time <- Sys.time()
      print(end_time-start_time)
    }
    Vect=0:(k-1)
    
    res_t <- list();
    if(Popfix){
      rho_=oldrho*sum(Boxr)
      rho_=rho_-(Boxr[1])
      rho_=10^(rho_)
      rho_=rho_*Rho
      if(SB){
        beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
      }
      if(SF){
        sigma=oldsigma*(Boxs[2]-Boxs[1])
        sigma=sigma+Boxs[1]
      }
      if(NC==1){
        builder=build_HMM_matrix(k,(rho_),beta,L=L,Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge )
      }
      if(NC>1){
        builder=list()
        for(chr in 1:NC){
          builder[[chr]]=build_HMM_matrix(k,(rho_[chr]),beta,L=L[chr],Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge )
        }
      }
    }
    if(!Popfix){
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        oldXi[x:xx]=oldXi_[ix]
      }
      Xi_=oldXi*sum(BoxP)
      Xi_=Xi_-(BoxP[1])
      Xi_=10^Xi_
      rho_=oldrho*sum(Boxr)
      rho_=rho_-(Boxr[1])
      rho_=10^(rho_)
      rho_=rho_*Rho
      if(SB){
        beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
      }
      if(SF){
        sigma=oldsigma*(Boxs[2]-Boxs[1])
        sigma=sigma+Boxs[1]
      }
      if(NC==1){
        builder=build_HMM_matrix(k,(rho_),beta,L=L,Pop=Pop,Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
      }
      if(NC>1){
        if(Share_r){
          builder=build_HMM_matrix(k,(rho_[1]),beta,L=L[1],Pop=Pop,Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
        }else{
          builder=list()
          for(chr in 1:NC){
            builder[[chr]]=build_HMM_matrix(k,(rho_[chr]),beta,L=L[chr],Pop=Pop,Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
          }
        }
        
      }
    }
    if(NC==1){
      Q = builder[[1]]
      nu= builder[[2]]
      Tc=builder[[3]]
      g=matrix(0,nrow=length(Tc),ncol=2)
      if(!FS){
        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
      }
      if(FS){
        a=0.8125
        
        g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
        g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
        
      }
      
      Tc_r=builder[[4]]
    }
    if(NC>1){
      
      if(Share_r){
        Q = builder[[1]]
        nu= builder[[2]]
        Tc=builder[[3]]
        g=matrix(0,nrow=length(Tc),ncol=2)
        if(!FS){
          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
        }
        if(FS){
          a=0.8125
          
          g[,2]= (a - (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
          g[,1]= a*((1-a) + (a*exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)))
          
        }
        
        Tc_r=builder[[4]]
      }else{
        Q=list()
        nu=list()
        Tc=list()
        Tc_r=list()
        g=list()
        for(chr in 1:NC){
          Q[[chr]] = builder[[chr]][[1]]
          nu[[chr]]= builder[[chr]][[2]]
          Tc[[chr]]=builder[[chr]][[3]]
          Tc_r[[chr]]=builder[[chr]][[4]]
          g[[chr]]=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
          if(!FS){
            g[[chr]][,2]=1-exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]])
            g[[chr]][,1]=exp(-mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]])
          }else{
            g[[chr]][,2]= 0.75 - (0.75*exp(-4*mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]]/3))
            g[[chr]][,1]= 0.25 + (0.75*exp(-4*mu[chr]*(beta+((1-beta)*mu_b))*2*Tc[[chr]]/3))
          }
          
        }
      }
      
    }
    res_t$LH<-LH
    res_t$Tc<-Tc_r
    Tc_or=Tc_r
    Tc_o=Tc
    if(!Popfix){
      Xi_t=oldXi_
      Xit=vector()
      pop_vect=get('pop_vect', envir=test.env)
      xx=0
      for(ix in 1:length(Xi_t)){
        x=xx+1
        xx = xx + pop_vect[ix]
        Xit[x:xx]=Xi_t[ix]
      }
      Xit=Xit*sum(BoxP)
      Xit=Xit-(BoxP[1])
      Xit=10^Xit
    }
    rho_t=oldrho*sum(Boxr)
    rho_t=rho_t-(Boxr[1])
    rho_t=10^(rho_t)
    rho_t=rho_t*Rho
    if(SB){
      beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1])^2
    }
    if(SF){
      sigma=oldsigma*(Boxs[2]-Boxs[1])
      sigma=sigma+Boxs[1]
    }
    if(!Popfix){
      res_t$Xi=Xit
    }
    res_t$mu<-mu
    res_t$beta=beta
    
    res_t$N=test.env$Big_Xi
    res_t$M=test.env$Big_M
    res_t$q_=test.env$q_
    
    res_t$sigma=sigma
    rho_t=rho_t/(2*L)
    res_t$rho=rho_t
    res_t$Beta<- Beta
    res_t$Self<- Self
    old_list[[length(old_list)+1]]=res_t
    if(!SB){
      oldBeta=Beta
    }
    if(SB){
      oldBeta=Beta
      Beta=beta
    }
    if(SF){
      oldSelf=Self
      Self=sigma
    }
    if(!SF){
      oldSelf=Self
    }
    if(Beta==oldBeta&Self==oldSelf){
      mb=maxBit
      print("no need  for more iteration")
    }
    gamma_t=rho_t/mu
    print(gamma_t)
    print(gamma)
    if(any(gamma_t<=0.5*gamma|gamma_t>=2*gamma)&redo_R){
      maxBit=maxBit+1
      gamma=gamma_t
      
    }
    
  }
  if(NC==1){
    
    if(!restart){
      LH=0
      test=Build_zip_Matrix_mailund(Q,g,Os[[1]][[2]],nu)
      for(i in 1:length(Os)){
        fo=forward_zip_mailund(Os[[i]][[1]],g,nu,test[[1]])
        LH=LH+(sum(fo[[2]]))
      }
      res<-list()
      res$LH=LH
      res$Tc=Tc_or
      Q=t(Q)
      #    res$Q <- Q;
      #    res$g <- g;
      rho_=oldrho*sum(Boxr)
      rho_=rho_-(Boxr[1])
      rho_=10^(rho_)
      rho_=rho_*Rho
      rho_=rho_/(2*L)
      if(SB){
        beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
      }
      if(SF){
        sigma=oldsigma*(Boxs[2]-Boxs[1])
        sigma=sigma+Boxs[1]
      }
      if(!Popfix){
        res$Xi=Xi_
      }
      #     res$Self<-oldSelf
      
      res$mu=mu
      res$L<-L
      res$beta=beta
      res$sigma=sigma
      res$rho=rho_
      #      res$Beta<- oldBeta
      
      if(FS){
        res$mu=(1-exp(log(1-mu)/Ne))*3/4
        res$Ne=Ne
      }
      
      res$Os<-Os
      res$N=test.env$Big_Xi
      res$M=test.env$Big_M
      res$q_=test.env$q_
      #     res$old<-old_list
    }
    
    if(restart){
      res<-old_list[[length(old_list)]]
      res$Os<-Os
      res$old<-old_list
      res$L<-L
    }
  }
  if(NC>1){
    LH=0
    if(Share_r){
      test=Build_zip_Matrix_mailund(Q,g,Os[[1]][[1]][[2]],nu)
      for(chr in 1:NC){
        for(i in 1:length(Os[[chr]])){
          fo=forward_zip_mailund(Os[[chr]][[i]][[1]],g,nu,test[[1]])
          LH=LH+(sum(fo[[2]]))
        }
      }
      Q=t(Q)
    }else{
      for(chr in 1:NC){
        test=Build_zip_Matrix_mailund(Q[[chr]],g[[chr]],Os[[chr]][[1]][[2]],nu[[chr]])
        for(i in 1:length(Os[[chr]])){
          fo=forward_zip_mailund(Os[[chr]][[i]][[1]],g[[chr]],nu[[chr]],test[[1]])
          LH=LH+(sum(fo[[2]]))
        }
        Q[[chr]]=t(Q[[chr]])
      }
      
    }
    res <- list();
    res$LH=LH
    if(length(Tc_or[[1]])>1){
      res$Tc=Tc_or[[1]]
    }else{
      res$Tc=Tc_or
    }
    
    #  res$Q <- Q;
    #  res$g <- g;
    res$L<-L
    rho_=oldrho*sum(Boxr)
    rho_=rho_-(Boxr[1])
    rho_=10^(rho_)
    rho_=rho_*Rho
    rho_=rho_/(2*L)
    if(SB){
      beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
    }
    if(SF){
      sigma=oldsigma*(Boxs[2]-Boxs[1])
      sigma=sigma+Boxs[1]
    }
    if(!Popfix){
      res$Xi=Xi_
    }
    res$beta=beta
    res$sigma=sigma
    res$rho=rho_
    res$N=N
    res$M=test.env$Big_M
    res$q_=test.env$q_
    # res$Self<-oldSelf
    
    res$mu=mu
    if(FS){
      res$mu=(1-exp(log(1-mu)/Ne))*3/4
      res$Ne=Ne
    }
    res$Os<-Os
    #  res$old<-old_list
  }
  return(res)
}


Baum_Welch_algo_t<-function(Os, maxIt =20,L,mu,theta_W,Rho,beta=1,Popfix=T,SB=F,k=20,BoxB=c(0.1,1),BoxP=c(3,3),Boxr=c(1,1),maxBit=1,pop_vect=NA,window_scaling=c(1,0),sigma=0.00,SF=F,Boxs=c(0,0.97),ER=F,BW=F,NC=1,redo_R=F,mu_b=1,SCALED=F,Big_Window=F,EM=T,FS=F){
  Xi=NA
  print(paste("sequence length :",L,sep=" "))
  if(length(Rho)>1&length(Rho)!=NC){
    stop("Problem in recombination definition")
  }
  old_list=list()
  n <- length(Os[[1]][[1]]);
  Pop=Popfix
  theta=mu*2*L
  gamma=Rho/theta
  if(NC>1&length(Rho)==1){
    Rho=rep(Rho,NC)
  }
  gamma_o=gamma
  print(gamma)
  test.env <- new.env()
  test.env$L <- L
  test.env$k <- k
  test.env$mu <- mu
  test.env$Rho <- Rho
  test.env$window_scaling <- window_scaling
  test.env$BW<-BW
  test.env$Pop<-Popfix
  test.env$NC<-NC
  test.env$FS<-FS
  test.env$mu_b <- mu_b
  test.env$Big_Window <- Big_Window
  if(NC>1){
    npair=NC
    test.env$npair <- NC
  }
  if(NC==1){
    npair=length(Os)
    test.env$npair <- length(Os)
  }
  mb=0
  if(any(!is.na(pop_vect))){
    Klink=length(pop_vect)
  }
  if((all(is.na(pop_vect))|sum(pop_vect)!=k)){
    Klink=0.5*k
    pop_vect=rep(2, Klink)
    print("Default pop discretizaion vector")
  }
  test.env$Klink <- Klink
  test.env$pop_vect <- pop_vect
  if(SB){
    BoxB[1]=max(sqrt(0.01),sqrt(BoxB[1]))
    BoxB[2]=min(sqrt(1),sqrt(BoxB[2]))
    beta=max((BoxB[1]^2),beta)
    beta=min(beta,(BoxB[2]^2))
    Beta=beta
    beta=rep(beta,length(pop_vect))
  }
  if(SF){
    sigma=min(Boxs[2],sigma)
    sigma=max(sigma,Boxs[1])
    Self=sigma
    sigma=rep(sigma,length(pop_vect))
  }
  if(ER){
    oldrho=rep((Boxr[1]/sum(Boxr)),length(pop_vect))
  }
  if(!SB&!SF){
    maxBit=1
  }
  if(!SB){
    Beta=beta
    test.env$beta <- beta
  }
  if(!SF){
    Self=sigma
    oldSelf=Self
    test.env$sigma <- sigma
  }
  if(!ER){
    Boxr=c(0,0)
  }
  while(mb<maxBit){
    print(paste("Beta:",Beta))
    print(paste("Self:",Self))
    test.env$Beta <- Beta
    test.env$Self <- Self
    test.env$BoxB <- BoxB
    test.env$Boxs <- Boxs
    mb=mb+1
    diff=1
    it <- 0
    if(mb==1){
      if(SB){
        oldbeta=(sqrt(beta)-BoxB[1])/(BoxB[2]-BoxB[1])
      }
      if(SF){
        oldsigma=(sigma-Boxs[1])/(Boxs[2]-Boxs[1])
      }
      
    }
    if(mb>1){
      if(SB){
        beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1])^2
        xx=0
        beta_=vector()
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          beta_[x:xx]=beta[ix]
        }
      }
      if(SF){
        sigma=oldsigma*(Boxs[2]-Boxs[1])
        sigma=sigma+Boxs[1]
        xx=0
        sigma_=vector()
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          sigma_[x:xx]=sigma[ix]
        }
      }
      if(!Popfix){
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          oldXi[x:xx]=oldXi_[ix]
        }
        Xi_=oldXi*sum(BoxP)
        Xi_=Xi_-(BoxP[1])
        Xi_=10^Xi_
      }
      
      if(NC==1){
        theta=(((theta_W*(beta^2)))*2/((2-sigma)*(beta+((1-beta)*mu_b))))
        mu=theta/(2*L)
        Rho=gamma*theta
      }
      if(NC>1){
        theta=mean(theta_W/(2*L))
        mu=(((theta*(beta^2)))*2/((2-sigma)*(beta+((1-beta)*mu_b))))
        Rho=gamma*theta*2*L
      }
      
      test.env$mu <- mu
      test.env$Rho <- Rho
    }
    if(!ER){
      oldrho=rep(0,Klink)
      if(NC>1){
        oldrho=list()
        for(chr in 1:chr){
          oldrho[[chr]]=rep(0,Klink)
        }
        
      }
    }
    if(ER){
      oldrho=rep((Boxr[1]/sum(Boxr)),length(pop_vect))
      if(NC>1){
        oldrho=list()
        for(rr in 1:NC){
          oldrho[[rr]]=rep((Boxr[1]/sum(Boxr)),length(pop_vect))
        }
      }
      
    }
    oldXi_=rep((BoxP[1]/sum(BoxP)),Klink)
    oldXi=vector()
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      oldXi[x:xx]=oldXi_[ix]
    }
    Do_BW=T
    restart=F
    diff_conv=vector()
    while (it<maxIt&!restart){
      start_time <- Sys.time()
      restart=F
      if(!Do_BW){
        it=0
        oldXi_=rep((BoxP[1]/sum(BoxP)),Klink)
        oldXi=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          oldXi[x:xx]=oldXi_[ix]
        }
      }
      it <- it+1;
      print(paste("It:",it))
      if(Popfix){
        if(!is.list(oldrho)){
          rho_=oldrho*sum(Boxr)
          rho_=rho_-(Boxr[1])
          rho_=10^(rho_)
          rho_=rho_*Rho
          rho=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            rho[x:xx]= rho_[ix]
          }
        }else{
          rho_=list()
          rho=list()
          for(rr in 1:NC){
            rho_[[rr]]=oldrho[[rr]]*sum(Boxr)
            rho_[[rr]]=rho_[[rr]]-(Boxr[1])
            rho_[[rr]]=10^(rho_[[rr]])
            rho_[[rr]]=rho_[[rr]]*Rho[rr]
            
            xx=0
            for(ix in 1:Klink){
              x=xx+1
              xx = xx + pop_vect[ix]
              rho[[rr]][x:xx]= rho_[[rr]][ix]
            }
          }
          
        }
        
        if(SB){
          beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
          beta=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            beta[x:xx]=beta_[ix]
          }
        }
        if(SF){
          sigma_=oldsigma*(Boxs[2]-Boxs[1])
          sigma_=sigma_+Boxs[1]
          sigma=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            sigma[x:xx]= sigma_[ix]
          }
        }
        print(c("sigma:",sigma,"beta :",beta))
        print(c("rho/theta:",rho_/theta))
        Keep_going=F
        if(it==1){
          diff_o=0
        }
        if(it>1){
          diff_o=max(abs(sigma_o-sigma),abs(beta_o-beta))
          count_diff_o=0
          if(diff_o>=0.005){
            if(it==maxIt){
              count_diff_o=count_diff_o+1
              maxIt=maxIt+1
              if(count_diff_o>10){
                maxBit=maxBit-1
              }
            }
          }
        }
        sigma_o=sigma
        beta_o=beta
        if(NC==1){
          builder=build_HMM_matrix_t(k,(rho),beta=beta,L=L,Pop=Pop,Xi=NA,Beta=Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
        }
        if(NC>1){
          builder=list()
          for(chr in 1:NC){
            builder[[chr]]=build_HMM_matrix_t(k,(rho[[chr]]),beta=beta,L=L[chr],Pop=Pop,Xi=NA,Beta=Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
          }
        }
      }
      if(!Popfix){
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          oldXi[x:xx]=oldXi_[ix]
        }
        Xi_=oldXi*sum(BoxP)
        Xi_=Xi_-(BoxP[1])
        Xi_=10^Xi_
        if(!is.list(oldrho)){
          rho_=oldrho*sum(Boxr)
          rho_=rho_-(Boxr[1])
          rho_=10^(rho_)
          rho_=rho_*Rho
          rho=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            rho[x:xx]= rho_[ix]
          }
        }else{
          rho_=list()
          rho=list()
          for(rr in 1:NC){
            rho_[[rr]]=oldrho[[rr]]*sum(Boxr)
            rho_[[rr]]=rho_[[rr]]-(Boxr[1])
            rho_[[rr]]=10^(rho_[[rr]])
            rho_[[rr]]=rho_[[rr]]*Rho[rr]
            
            xx=0
            for(ix in 1:Klink){
              x=xx+1
              xx = xx + pop_vect[ix]
              rho[[rr]][x:xx]= rho_[[rr]][ix]
            }
          }
          
        }
        if(SB){
          beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
          beta=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            beta[x:xx]=beta_[ix]
          }
        }
        if(SF){
          sigma_=oldsigma*(Boxs[2]-Boxs[1])
          sigma_=sigma_+Boxs[1]
          sigma=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            sigma[x:xx]= sigma_[ix]
          }
        }
        print(c("sigma:",sigma,"beta :",beta))
        print(c("rho/theta:",rho/theta))
        Keep_going=F
        if(it==1){
          diff_o=0
        }
        if(it>1){
          diff_o=max(abs(sigma_o-sigma),abs(beta_o-beta))
          count_diff_o=0
          if(diff_o>=0.005){
            if(it==maxIt){
              count_diff_o=count_diff_o+1
              maxIt=maxIt+1
              if(count_diff_o>10){
                maxBit=maxBit-1
              }
            }
          }
        }
        sigma_o=sigma
        beta_o=beta
        
        if(NC==1){
          builder=build_HMM_matrix_t(k,(rho),beta=beta,L=L,Pop=Pop,Xi_,Beta=Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
        }
        if(NC>1){
          builder=list()
          for(chr in 1:NC){
            builder[[chr]]=build_HMM_matrix_t(k,(rho[[chr]]),beta=beta,L=L[chr],Pop=Pop,Xi_,Beta=Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
          }
        }
      }
      if(NC==1){
        Q = builder[[1]]
        nu= builder[[2]]
        Tc=builder[[3]]
        g=matrix(0,nrow=length(Tc),ncol=2)
        if(!FS){
          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
        }
        if(FS){
          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
        }
        M=matrix(0,nrow=length(Tc),ncol=2)
        N=matrix(0,length(Tc),length(Tc))
        M=matrix(0,nrow=length(Tc),ncol=2)
        N=matrix(0,length(Tc),length(Tc))
        MLH=0
        q_=rep(0,length(Tc))
        s_t=Sys.time()
        test=Build_zip_Matrix_mailund(Q,g,Os[[1]][[2]],nu)
        e_t=Sys.time()
        print(e_t-s_t)
        
        for(i in 1:length(Os)){
          fo=forward_zip_mailund(Os[[i]][[1]],g,nu,test[[1]])
          MLH=MLH+fo[[3]]
          c=exp(fo[[2]])
          ba=Backward_zip_mailund(Os[[i]][[1]],test[[3]],length(Tc),c)
          W=list()
          int=t(Q)%*%diag(g[,1])
          int=eigen(int)
          W$P=int$vectors
          W$P_=solve(W$P)
          
          symbol= Os[[i]][[3]][,1]
          for(ob in 1){
            truc_M=matrix(0,nrow=length(Tc),ncol=2)
            if(Os[[i]][[1]][(ob+1)]<=1){
              truc_N=(fo[[1]][,ob]%*%(t(ba[,(ob+1)]*g[,(as.numeric(Os[[i]][[1]][(ob+1)])+1)]/c[(ob+1)])))
              truc_M[,(as.numeric(Os[[i]][[1]][(ob+1)])+1)]=fo[[1]][,ob]*(test[[2]][[(as.numeric(Os[[i]][[1]][(ob+1)])+1)]]%*%(ba[,(ob+1)]/c[(ob+1)]))
            }
            if(Os[[i]][[1]][(ob+1)]>1){
              truc_N=(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_)*test[[4]][[(as.numeric(Os[[i]][[1]][(ob+1)])+1)]])%*%(t(W$P))%*%diag(g[,1]))
              truc_M[,1]=diag(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_))*test[[2]][[(as.numeric(Os[[i]][[1]][(ob+1)])+1)]]%*%t(W$P))
            }
            #truc_N=sum_chi_corrected(fo[[1]][,ob],ba[,(ob+1)],test[[4]],Os[[i]][[1]][(ob+1)],g,W,c[(ob+1)])
            N=N+truc_N
            #truc_M=sum_M_cor(fo[[1]][,ob],ba[,(ob+1)],Os[[i]][[1]][(ob+1)],test[[2]],W,(c[(ob+1)]),test[[1]])
            M=M+truc_M
          }
          for(sym in symbol){
            ob=as.numeric(sym)
            pos=which(as.numeric(Os[[i]][[1]][-c(1,length(Os[[i]][[1]]))])==ob)
            pos=pos+1
            if(ob<2){
              ba_t=t(t(ba[,(pos)])/c[(pos)])
              truc=c(rowSums(fo[[1]][,(pos-1)]*(test[[1]][[(ob+1)]]%*%ba_t)))
              truc=(truc/(sum(truc)))*length(pos)
              M[,(ob+1)]=M[,(ob+1)]+ truc
              N=N+(fo[[1]][,(pos-1)]%*%(t(diag(g[,(ob+1)])%*%ba_t)))
            }
            if(ob>=2){
              A=(t(W$P)%*%fo[[1]][,(pos-1)]%*%t(t(t(ba[,(pos)])/c[(pos)]))%*%t(W$P_))
              A_=A*test[[2]][[(ob+1)]]
              A_=(t(W$P_)%*%A_%*%t(W$P))
              M[,1]=M[,1]+(diag(A_))
              A_=A*test[[4]][[(ob+1)]]
              A_=(t(W$P_)%*%A_%*%t(W$P))
              N=N+((A_%*%diag(g[,1])))
            }
          }
          if(as.numeric(Os[[i]][[1]][length(Os[[i]][[1]])])==1){
            M[,2]=M[,2]+(fo[[1]][,length(Os[[i]][[1]])]*ba[,length(Os[[i]][[1]])])
          }
          if(as.numeric(Os[[i]][[1]][length(Os[[i]][[1]])])!=1){
            M[,1]=M[,1]+((fo[[1]][,length(Os[[i]][[1]])]*ba[,length(Os[[i]][[1]])])/sum(fo[[1]][,length(Os[[i]][[1]])]*ba[,length(Os[[i]][[1]])]))
          }
          q_=q_+((fo[[1]][,1]*ba[,(1)])/sum(fo[[1]][,1]*ba[,(1)]))
        }
        N=N*t(Q)
        if(is.complex(N)){
          N=Re(N)
        }
        
        if(is.complex(M)){
          M=Re(M)
        }
        Scale_N=(L-1)/sum(N)
        N=N*Scale_N
        Scale_M=L/sum(M)
        M=M*Scale_M
        q_=q_/sum(q_)
        if(SCALED){
          corrector_M=rowSums(M)
          M=diag(1/corrector_M)%*%M
          corrector_N=rowSums(N)
          N=diag(1/corrector_N)%*%N
        }
        
        
      }
      if(NC>1){
        Q=list()
        nu=list()
        Tc=list()
        g=list()
        M=list()
        N=list()
        MLH=list()
        q_=list()
        for(chr in 1:NC){
          Q[[chr]] = builder[[chr]][[1]]
          nu[[chr]]= builder[[chr]][[2]]
          Tc[[chr]]=builder[[chr]][[3]]
          g[[chr]]=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
          if(!FS){
            g[[chr]][,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]])
            g[[chr]][,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]])          }
          if(FS){
            g[[chr]][,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]]/3))
            g[[chr]][,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]]/3))
          }
          
          
          M[[chr]]=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
          N[[chr]]=matrix(0,length(Tc[[chr]]),length(Tc[[chr]]))
          MLH[[chr]]=0
          q_[[chr]]=rep(0,length(Tc[[chr]]))
          test=Build_zip_Matrix_mailund(Q[[chr]],g[[chr]],Os[[chr]][[1]][[2]],nu[[chr]])
          
          for(i in 1:length(Os[[chr]])){
            fo=forward_zip_mailund(Os[[chr]][[i]][[1]],g[[chr]],nu[[chr]],test[[1]])
            MLH[[chr]]=MLH[[chr]]+fo[[3]]
            c=exp(fo[[2]])
            ba=Backward_zip_mailund(Os[[chr]][[i]][[1]],test[[2]],length(Tc[[chr]]),c)
            W=list()
            int=t(Q[[chr]])%*%diag(g[[chr]][,1])
            int=eigen(int)
            W$P=int$vectors
            W$P_=solve(W$P)
            for(ob in 1){
              truc_M=matrix(0,nrow=length(Tc),ncol=2)
              if(Os[[chr]][[i]][[1]][(ob+1)]<=1){
                truc_N=(fo[[1]][,ob]%*%(t(ba[,(ob+1)]*g[[chr]][,(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]/c[(ob+1)])))
                truc_M[,(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]=fo[[1]][,ob]*(test[[2]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]]%*%(ba[,(ob+1)]/c[(ob+1)]))
              }
              if(Os[[chr]][[i]][[1]][(ob+1)]>1){
                truc_N=(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_)*test[[4]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]])%*%(t(W$P))%*%diag(g[[chr]][,1]))
                truc_M[,1]=diag(t(W$P_)%*%(t(W$P)%*%fo[[1]][,ob]%*%t(ba[,(ob+1)]/c[(ob+1)])%*%t(W$P_))*test[[2]][[(as.numeric(Os[[chr]][[i]][[1]][(ob+1)])+1)]]%*%t(W$P))
              }
              N=N+truc_N
              M=M+truc_M
            }
            symbol= Os[[chr]][[i]][[3]][,1]
            for(sym in symbol){
              ob=as.numeric(sym)
              pos=which(as.numeric(Os[[chr]][[i]][[1]][-c(1,length(Os[[chr]][[i]][[1]]))])==ob)
              pos=pos+1
              if(ob<2){
                ba_t=t(t(ba[,(pos)])/c[(pos)])
                truc=c(rowSums(fo[[1]][,(pos-1)]*(test[[1]][[(ob+1)]]%*%ba_t)))
                truc=(truc/(sum(truc)))*length(pos)
                M[[chr]][,(ob+1)]=M[[chr]][,(ob+1)]+ truc
                N[[chr]]=N[[chr]]+(fo[[1]][,(pos-1)]%*%(t(diag(g[[chr]][,(ob+1)])%*%ba_t)))
              }
              if(ob>=2){
                A=(t(W$P)%*%fo[[1]][,(pos-1)]%*%t(t(t(ba[,(pos)])/c[(pos)]))%*%t(W$P_))
                A_=A*test[[2]][[(ob+1)]]
                A_=(t(W$P_)%*%A_%*%t(W$P))
                M[[chr]][,1]=M[[chr]][,1]+(diag(A_))
                A_=A*test[[4]][[(ob+1)]]
                A_=(t(W$P_)%*%A_%*%t(W$P))
                N[[chr]]=N[[chr]]+((A_%*%diag(g[[chr]][,1])))
              }
            }
            if(as.numeric(Os[[chr]][[i]][[1]][length(Os[[chr]][[i]][[1]])])==1){
              M[[chr]][,2]=M[[chr]][,2]+(fo[[1]][,length(Os[[chr]][[i]][[1]])]*(ba[,length(Os[[chr]][[i]][[1]])]/(c[length(Os[[chr]][[i]][[1]])])))
            }
            if(as.numeric(Os[[chr]][[i]][[1]][length(Os[[chr]][[i]][[1]])])!=1){
              M[[chr]][,1]=M[[chr]][,1]+(fo[[1]][,length(Os[[chr]][[i]][[1]])]*(ba[,length(Os[[chr]][[i]][[1]])]/(c[length(Os[[chr]][[i]][[1]])])))
            }
            q_[[chr]]=q_[[chr]]+((fo[[1]][,1]*ba[,(1)])/sum(fo[[1]][,1]*ba[,(1)]))
            
          }
          
          N[[chr]]=N[[chr]]*t(Q[[chr]])
          if(is.complex(N[[chr]])){
            N[[chr]]=Re(N[[chr]])
          }
          
          if(is.complex(M[[chr]])){
            M[[chr]]=Re(M[[chr]])
          }
          Scale_N=(L[chr]-1)/sum(N[[chr]])
          N[[chr]]=N[[chr]]*Scale_N
          
          Scale_M=L[chr]/sum(M[[chr]])
          M[[chr]]=M[[chr]]*Scale_M
          q_[[chr]]=q_[[chr]]/sum(q_[[chr]])
          if(SCALED){
            corrector_M=rowSums(M[[chr]])
            M[[chr]]=diag(1/corrector_M)%*%M[[chr]]
            corrector_N=rowSums(N[[chr]])
            N[[chr]]=diag(1/corrector_N)%*%N[[chr]]
          }
          
          
        }
      }
      if(it>1){
        print(paste(" old Likelihood: ",oldMLH))
      }
      if(NC==1){
        print(paste(" New Likelihood: ",MLH))
        
        if(it>1){
          
          if(oldMLH > MLH|MLH=="NaN"){
            
            if(it>2){
              restart=T
              if(!Popfix){
                oldXi_=oldXi_s
              }
              if(ER){
                oldrho=oldrho_s
              }
              if(SB){
                oldbeta=oldbeta_s
              }
              if(SF){
                oldsigma=oldsigma_s
              }
            }
            if(it==2){
              print("Algortihm has converge but may need more data for better results ")
              restart=T
              MB=maxBit
            }
          }
          
        }
        oldMLH=MLH
      }
      if(NC>1){
        MLH1=0
        for(chr in 1:length(MLH)){
          MLH1=MLH[[chr]]+MLH1
        }
        print(paste("New Likelihood: ",MLH1))
        if(it>1){
          if(oldMLH > MLH1|MLH1=="NaN"){
            if(it>1){
              restart=T
              if(!Popfix){
                oldXi_=oldXi_s
              }
              if(ER){
                oldrho=oldrho_s
              }
              if(SB){
                oldbeta=oldbeta_s
              }
              if(SF){
                oldsigma=oldsigma_s
              }
            }
            if(it==2){
              print("Algortihm has converge but may need more data for better results ")
              restart=T
              MB=maxBit
            }
          }
          
        }
        oldMLH=MLH1
      }
      if(!restart){
        
        if(NC==1){
          test.env$Big_Xi <- N
          test.env$Big_M <-M
          test.env$q_ <-q_
        }
        if(NC>1){
          test.env$Big_Xi <- N
          test.env$Big_M <-M
          test.env$q_ <-q_
        }
        Do_BW=T
        if(!Popfix){
          oldXi_s=oldXi_
        }
        if(ER){
          oldrho_s=oldrho
        }
        if(SB){
          oldbeta_s=oldbeta
        }
        if(SF){
          oldsigma_s=oldsigma
        }
        lr=length(oldrho)
        test.env$lr<-lr
        if(Do_BW){
          if(NC==1){
            if(ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(1+(2*Klink)):(3*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      Pop=get('Pop', envir=test.env)
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      beta_=vector()
                      rho_=vector()
                      sigma_=vector()
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        beta_[x:xx]=beta[ix]
                        rho_[x:xx]=rho[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      
                      builder=build_HMM_matrix_t(n,rho_,beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair) #
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:3*Klink])
                    rho=sol[1:Klink]
                    beta_=sol[(1+Klink):(2*Klink)]
                    sigma_=sol[(1+(2*Klink)):(3*Klink)]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[(1+(2*Klink)):(3*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[(1+(3*Klink)):(4*Klink)]
                      rho_=vector()
                      beta_=vector()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        rho_=rho[ix]
                        beta_=beta[ix]
                        sigma_=sigma[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,rho_,beta=beta_,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(4*Klink),1])
                    rho=sol[1:Klink]
                    beta_=sol[(1+Klink):(2*Klink)]
                    sigma_=sol[(1+(2*Klink)):(3*Klink)]
                    Xi_=sol[(1+(3*Klink)):(4*Klink)]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma',envir = test.env)
                      Self=get('Self',envir = test.env)
                      Klink=get('Klink',envir = test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      rho_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        rho_=rho[ix]
                        beta_=beta[ix]
                      }
                      builder=build_HMM_matrix_t(n,rho_,beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*KLink),1])
                    rho=sol[1:Klink]
                    beta_=sol[(1+Klink):(2*Klink)]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_)))
                    oldrho=rho
                    oldbeta=beta_
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma',envir = test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BETA=get('BETA',envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(1+(2*Klink)):(3*Klink)]
                      Xi=vector()
                      rho_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        rho_[x:xx]=rho[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(3*(Klink)),1])
                    rho=sol[1:Klink]
                    beta_=sol[(1+(1*Klink)):(2*Klink)]
                    Xi_=sol[(1+(2*Klink)):(3*Klink)]
                    diff=max(abs(c(rho- oldrho,oldbeta-beta_,Xi_-oldXi_)))
                    oldrho=rho
                    oldbeta=beta_
                    oldXi_=Xi_
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(1+Klink):(2*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      sigma_=vector()
                      rho_=vector()
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        
                        rho_[x:xx]=rho[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldsigma), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*Klink),1])
                    rho=sol[1:Klink]
                    sigma_=sol[(1+Klink):(2*Klink)]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho,oldsigma-sigma_)))
                    oldrho=rho
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      Boxs=get('Boxs', envir=test.env)
                      sigma=param[(1+(1*Klink)):(2*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(1+(2*Klink)):(3*Klink)]
                      Xi=vector()
                      sigma_=vector()
                      rho_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        sigma_[x:xx]=sigma[ix]
                        rho_[x:xx]=rho[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      
                      Self=get('Self', envir=test.env)
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(3*(Klink)),1])
                    rho=sol[1:Klink]
                    sigma_=sol[(1+(1*Klink)):(2*Klink)]
                    Xi_=sol[(1+(2*Klink)):(3*Klink)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_,oldsigma-sigma_)))
                    oldrho=rho
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      rho_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        rho_[x:xx]=rho[ix]
                      }
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:Klink,1])
                    rho=sol[1:Klink]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(rho- oldrho)))
                    oldrho=rho
                  }
                  
                  if(!Popfix){
                    
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      rho=param[1:Klink]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(1+Klink):(2*Klink)]
                      Xi=vector()
                      rho_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:length(Xi_)){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        rho_[x:xx]=rho[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,Sigma=Self,scale=window_scaling,FS=FS,sigma = sigma,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldrho,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*Klink),1])
                    rho=sol[1:Klink]
                    Xi_=sol[(1+Klink):(2*Klink)]
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
                    oldrho=rho
                    oldXi_=Xi_
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
            }
            if(!ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(1+Klink):(2*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      beta_=vector()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*Klink),1])
                    beta_=sol[1:Klink]
                    sigma_=sol[(1+Klink):(2*Klink)]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_)))
                    oldbeta=beta_
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[(1+Klink):(2*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[(1+(2*Klink)):(3*Klink)]
                      Xi=vector()
                      sigma_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        sigma_[x:xx]=sigma[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(3*(Klink)),1])
                    beta_=sol[1:Klink]
                    sigma_=sol[(1+(1*Klink)):(2*Klink)]
                    Xi_=sol[(1+(2*Klink)):(3*Klink)]
                    diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:Klink,1])
                    beta_=sol[1:Klink]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(oldbeta-beta_)))
                    oldbeta=beta_
                    
                  }
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(1+Klink):(2*Klink)]
                      Xi=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair) #
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*(Klink)),1])
                    beta_=sol[1:Klink]
                    Xi_=sol[(1+Klink):(2*Klink)]
                    diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1:Klink]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:Klink,1])
                    sigma_=sol[1:Klink]
                    print(paste(" new Complete likelihood : ", LH ))
                    diff=max(abs(c(oldsigma-sigma_)))
                    oldsigma=sigma_
                    
                  }
                  
                  if(!Popfix){
                    
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1:Klink]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[(1+Klink):(2*Klink)]
                      Xi=vector()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma,oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*(Klink)),1])
                    sigma_=sol[1:Klink]
                    Xi_=sol[(1+Klink):(2*Klink)]
                    
                    diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
                    oldXi_=Xi_
                    oldsigma=sigma_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(!Popfix){
                    function_to_minimize<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[1:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      Self=get('Self', envir=test.env)
                      
                      builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                      Q=builder[[1]]
                      Q=t(Q)
                      A=as.vector(Q)
                      keep=which(A>0&as.vector(Big_Xi)>0)
                      A=A[keep]
                      Big_Xi=as.vector(Big_Xi)
                      Big_Xi=Big_Xi[keep]
                      Big_M=get('Big_M', envir=test.env)
                      Tc=builder[[3]]
                      g=matrix(0,nrow=length(Tc),ncol=2)
                      if(!FS){
                        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                      }
                      if(FS){
                        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                      }
                      x=as.vector(g)
                      keep=which(x>0)
                      x=x[keep]
                      m=as.vector(Big_M)
                      m=m[keep]
                      q_=get('q_', envir=test.env)
                      nu=builder[[2]]
                      if(BW){
                        LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
                      }
                      if(!BW){
                        LH=-sum(log(A)*Big_Xi)
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldXi_)),M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((Klink)),1])
                    Xi_=sol[1:length(sol)]
                    
                    diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
                    oldXi_=Xi_
                    
                    print(paste(" new Complete likelihood : ",  -LH ))
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
            }
          }
          if(NC>1){
            if(ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      beta_=vector()
                      sigma_=vector()
                      rho_=list()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                        beta_[x:xx]=beta[ix]
                        
                      }
                      LH=0
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        
                        
                        builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                      }
                      return(LH)
                    }
                    sol=BBoptim(c(unlist(oldrho),oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((NC+2)*Klink),1])
                    rho=sol[1:(NC*Klink)]
                    beta_=sol[(NC*Klink+1):((NC+1)*Klink)]
                    sigma_=sol[((NC+1)*Klink+1):((NC+2)*Klink)]
                    diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_,oldsigma-sigma_)))
                    
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    
                    oldbeta=beta_
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      Xi_=param[(Klink*(Nc+2)+1):(Klink*(Nc+3))]
                      Xi=vector()
                      sigma_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        sigma_[x:xx]=sigma[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(unlist(oldrho),oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.numeric(sol[1:((NC+2)*Klink),1])
                    rho=sol[1:(NC*Klink)]
                    beta_=sol[(NC*Klink+1):((NC+1)*Klink)]
                    sigma_=sol[((NC+1)*Klink+1):((NC+2)*Klink)]
                    Xi_=sol[((NC+2)*Klink+1):((NC+3)*Klink)]
                    diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      sigma=get('sigma', envir=test.env)
                      Self=get('Self', envir=test.env)
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      Xi_=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
                      Xi=vector()
                      sigma_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        
                        builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])}
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(unlist(oldrho),oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((NC+2)*Klink),1])
                    rho=sol[1:((NC)*Klink)]
                    beta_=sol[((NC)*Klink+1):((NC+1)*Klink)]
                    diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_)))
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    oldbeta=beta_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      Xi_=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
                      Xi=vector()
                      sigma_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(unlist(oldrho),oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(Klink*(NC+2)),1])
                    rho=sol[1:(Klink*NC)]
                    beta_=sol[(Klink*NC+1):(NC+1)*Klink]
                    Xi_=sol[(Klink*(NC+1)+1):(NC+2)*Klink]
                    diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_,Xi_-oldXi_)))
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    oldbeta=beta_
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      beta=get('beta', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(Klink*(Nc)+1):(Klink*(Nc+1))]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        
                        builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(unlist(oldrho),oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((NC+1)*Klink),1])
                    rho=sol[1:(NC*Klink)]
                    sigma_=sol[((NC)*Klink+1):((NC+1)*Klink)]
                    diff=max(abs(c(rho-unlist(oldrho),oldsigma-sigma_)))
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      
                      beta=get('beta', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(Klink*(Nc)+1):(Klink*(Nc+1))]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      Xi_=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
                      Xi=vector()
                      sigma_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(unlist(oldrho),oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((NC+2)*(Klink)),1])
                    rho=sol[1:(NC*Klink)]
                    sigma_=sol[(Klink*NC+1):((NC+1)*Klink)]
                    Xi_=sol[(Klink*(NC+1)+1):((NC+2)*Klink)]
                    diff=max(abs(c(rho-unlist(oldrho),Xi_-oldXi_,oldsigma-sigma_)))
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      beta=get('beta', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      pop_vect=get('pop_vect', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        for( chr in 1:NC){
                          builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                          Q=builder[[1]]
                          Q=t(Q)
                          A=as.vector(Q)
                          keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                          A=A[keep]
                          Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                          Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                          Tc=builder[[3]]
                          g=matrix(0,nrow=length(Tc),ncol=2)
                          if(!FS){
                            g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                            g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          }
                          if(FS){
                            g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                            g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          }
                          x=as.vector(g)
                          keep=which(x>0)
                          x=x[keep]
                          m=as.vector(Big_M[[chr]])
                          m=m[keep]
                          q_=get('q_', envir=test.env)
                          nu=builder[[2]]
                          if(BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                          }
                          if(!BW){
                            LH=LH-sum(log(A)*Big_Xi[[chr]])
                          }
                          
                        }
                        return(LH)
                      }
                      sol= BBoptim(c(unlist(oldrho)),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                      LH=as.numeric(as.matrix(sol[[2]]))
                      sol=as.matrix(sol[[1]])
                      sol=as.numeric(sol[1:(NC*Klink),1])
                      rho=sol[1:(NC*Klink)]
                      diff=max(abs(c(rho- unlist(oldrho))))
                      for(chr in 1:NC){
                        xx=1+(chr-1)*Klink
                        yy=chr*Klink
                        oldrho[[chr]]=rho[xx:yy]
                      }
                    }
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      
                      beta=get('beta', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      rho=param[1:(Klink*Nc)]
                      rho=rho*sum(Boxr)
                      rho=rho-(Boxr[1])
                      rho=10^(rho)
                      Rho=get('Rho', envir=test.env)
                      rho=rho*Rho
                      BoxB=get('BoxB', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      rho_=list()
                      Xi_=param[(Klink*(Nc)+1):(Klink*(Nc+1))]
                      Xi=vector()
                      sigma_=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        rho_[[chr]]=vector()
                        xx=0
                        for(ix in 1:Klink){
                          x=xx+1
                          xx = xx + pop_vect[ix]
                          xr=(chr-1)*Klink+ix
                          rho_[[chr]][x:xx]=rho[xr]
                          
                        }
                        builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(unlist(oldrho),oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((NC+1)*(Klink)),1])
                    rho=sol[1:(NC*Klink)]
                    Xi_=sol[(NC*Klink+1):((NC+1)*(Klink))]
                    diff=max(abs(c(rho- unlist(oldrho),Xi_-oldXi_)))
                    for(chr in 1:NC){
                      xx=1+(chr-1)*Klink
                      yy=chr*Klink
                      oldrho[[chr]]=rho[xx:yy]
                    }
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
            }
            if(!ER){
              if(SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Big_M=get('Big_M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[(Klink+1):(2*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      beta_=vector()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                        beta_[x:xx]=beta[ix]
                        
                      }
                      LH=0
                      for(chr in 1:NC){
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*Klink),1])
                    beta_=sol[1:Klink]
                    sigma_=sol[(1+Klink*NC):(2*Klink)]
                    diff=max(abs(c(oldbeta- beta,oldsigma-sigma)))
                    oldbeta=beta_
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC',envir = test.env)
                      Klink=get('Klink',envir = test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      sigma=param[(Klink+1):(2*Klink)]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[((2*Klink+1)):(3*Klink)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      
                      beta_=vector()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                        beta_[x:xx]=beta[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Self=get('Self', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Big_M=get('Big_M', envir=test.env)
                      
                      for(chr in 1:NC){
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(3*(Klink)),1])
                    beta_=sol[1:Klink]
                    sigma_=sol[(1+Klink):(2*Klink)]
                    Xi_=sol[((2*Klink)+1):(3*Klink)]
                    diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir = test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      M=get('M', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      for(chr in 1:NC){
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:Klink,1])
                    beta_=sol[1:Klink]
                    diff=max(abs(c(oldbeta-beta_)))
                    oldbeta=beta_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC',envir = test.env)
                      rho_=Rho
                      BoxB=get('BoxB', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[(1+Klink):(2*Klink)]
                      Xi=vector()
                      beta_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        beta_[x:xx]=beta[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Big_M=get('Big_M', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      
                      for(chr in 1:NC){
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*(Klink)),1])
                    beta_=sol[1:Klink]
                    Xi_=sol[(1+Klink):(2*(Klink))]
                    diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
                    oldbeta=beta_
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                  }
                }
              }
              if(!SB){
                if(SF){
                  if(Popfix){
                    function_to_minimize <-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1:Klink]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      sigma_=vector()
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      for(chr in 1:NC){
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:Klink,1])
                    sigma_=sol[1:Klink]
                    diff=max(abs(c(oldsigma-sigma_)))
                    oldsigma=sigma_
                  }
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      NC=get('NC',envir = test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Klink=get('Klink', envir=test.env)
                      Boxs=get('Boxs', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=param[1:Klink]
                      sigma=sigma*(Boxs[2]-Boxs[1])
                      sigma=sigma+Boxs[1]
                      
                      Xi_=param[(Klink+1):(2*Klink)]
                      Xi=vector()
                      sigma_=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                        sigma_[x:xx]=sigma[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      
                      Big_M=get('Big_M', envir=test.env)
                      for(chr in 1:NC){
                        builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:(2*(Klink)),1])
                    sigma_=sol[1:Klink]
                    Xi_=sol[(1+Klink):(2*Klink)]
                    diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
                    oldXi_=Xi_
                    oldsigma=sigma_
                    print(sigma_)
                    print(paste("Xi:",oldXi_))
                  }
                }
                if(!SF){
                  if(!Popfix){
                    function_to_minimize_optim<-function(param){
                      Boxr=get('Boxr', envir=test.env)
                      mu=get('mu', envir=test.env)
                      npair=get('npair', envir=test.env)
                      Big_Window=get('Big_Window', envir=test.env)
                      mu_b=get('mu_b', envir=test.env)
                      FS=get('FS', envir=test.env)
                      Rho=get('Rho', envir=test.env)
                      NC=get('NC',envir=test.env)
                      rho_=Rho
                      beta=get('beta', envir=test.env)
                      BoxP=get('BoxP', envir=test.env)
                      Xi_=param[1:length(param)]
                      Xi=vector()
                      pop_vect=get('pop_vect', envir=test.env)
                      xx=0
                      for(ix in 1:Klink){
                        x=xx+1
                        xx = xx + pop_vect[ix]
                        Xi[x:xx]=Xi_[ix]
                      }
                      Xi=Xi*sum(BoxP)
                      Xi=Xi-(BoxP[1])
                      Xi=10^Xi
                      L=get('L', envir=test.env)
                      n=get('k', envir=test.env)
                      Big_Xi=get('Big_Xi', envir=test.env)
                      Beta=get('Beta', envir=test.env)
                      window_scaling=get('window_scaling', envir=test.env)
                      BW=get('BW', envir=test.env)
                      
                      Pop=get('Pop', envir=test.env)
                      LH=0
                      Big_M=get('Big_M', envir=test.env)
                      Self=get('Self', envir=test.env)
                      sigma=get('sigma', envir=test.env)
                      cut_edge=get('cut_edge', envir=test.env)
                      for(chr in 1:NC){
                        builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                        Q=builder[[1]]
                        Q=t(Q)
                        A=as.vector(Q)
                        keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                        A=A[keep]
                        Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                        Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                        Tc=builder[[3]]
                        g=matrix(0,nrow=length(Tc),ncol=2)
                        if(!FS){
                          g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                          g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                        }
                        if(FS){
                          g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                          g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                        }
                        x=as.vector(g)
                        keep=which(x>0)
                        x=x[keep]
                        m=as.vector(Big_M[[chr]])
                        m=m[keep]
                        q_=get('q_', envir=test.env)
                        nu=builder[[2]]
                        if(BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                        }
                        if(!BW){
                          LH=LH-sum(log(A)*Big_Xi[[chr]])
                        }
                        
                      }
                      return(LH)
                    }
                    sol= BBoptim(c(oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
                    LH=as.numeric(as.matrix(sol[[2]]))
                    sol=as.matrix(sol[[1]])
                    sol=as.numeric(sol[1:((Klink)),1])
                    Xi_=sol[1:length(sol)]
                    diff=max(abs(Xi_-oldXi_))
                    oldXi_=Xi_
                    print(paste("Xi:",oldXi_))
                    
                  }
                }
              }
            }
          }
        }
      }
      if(diff_o>=0.003){
        diff=(max(diff_o,diff))
      }
      diff_conv=c(diff_conv,diff_o)
      end_time <- Sys.time()
      print(end_time-start_time)
    }
    res_t <- list();
    if(Popfix){
      if(NC==1){
        rho_=oldrho*sum(Boxr)
        rho_=rho_-(Boxr[1])
        rho_=10^(rho_)
        rho_=rho_*Rho
        rho=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          rho[x:xx]=rho_[ix]
        }
        
      }
      if(NC>1){
        rho_=list()
        rho=list()
        for(chr in 1:NC){
          rho_[[chr]]=oldrho[[chr]]*sum(Boxr)
          rho_[[chr]]=rho_[[chr]]-(Boxr[1])
          rho_[[chr]]=10^(rho_[[chr]])
          rho_[[chr]]=rho_[[chr]]*Rho[chr]
          rho[[chr]]=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            rho[[chr]][x:xx]=rho_[[chr]][ix]
          }
        }
      }
      if(SB){
        beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
        beta=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          beta[x:xx]=beta_[ix]
        }
      }
      if(SF){
        sigma_=oldsigma*(Boxs[2]-Boxs[1])
        sigma_=sigma_+Boxs[1]
        sigma=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          sigma[x:xx]= sigma_[ix]
        }
      }
      if(NC==1){
        builder=build_HMM_matrix_t(k,(rho),beta=beta,L=L,Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair )
      }
      if(NC>1){
        builder=list()
        for(chr in 1:NC){
          builder[[chr]]=build_HMM_matrix_t(k,(rho[chr]),beta=beta,L=L[chr],Pop=Pop,Xi=NA,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS ,Big_Window=Big_Window,npair=npair)
        }
      }
    }
    if(!Popfix){
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        oldXi[x:xx]=oldXi_[ix]
      }
      Xi_=oldXi*sum(BoxP)
      Xi_=Xi_-(BoxP[1])
      Xi_=10^Xi_
      
      if(NC==1){
        rho_=oldrho*sum(Boxr)
        rho_=rho_-(Boxr[1])
        rho_=10^(rho_)
        rho_=rho_*Rho
        rho=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          rho[x:xx]=rho_[ix]
        }
        
      }
      if(NC>1){
        rho_=list()
        rho=list()
        for(chr in 1:NC){
          rho_[[chr]]=oldrho[[chr]]*sum(Boxr)
          rho_[[chr]]=rho_[[chr]]-(Boxr[1])
          rho_[[chr]]=10^(rho_[[chr]])
          rho_[[chr]]=rho_[[chr]]*Rho[chr]
          rho[[chr]]=vector()
          xx=0
          for(ix in 1:Klink){
            x=xx+1
            xx = xx + pop_vect[ix]
            rho[[chr]][x:xx]=rho_[[chr]][ix]
          }
        }
      }
      if(SB){
        beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
        beta=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          beta[x:xx]=beta_[ix]
        }
      }
      if(SF){
        sigma_=oldsigma*(Boxs[2]-Boxs[1])
        sigma_=sigma_+Boxs[1]
        sigma=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          sigma[x:xx]= sigma_[ix]
        }
      }
      if(NC==1){
        builder=build_HMM_matrix_t(k,(rho),beta=beta,L=L,Pop=Pop,Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
      }
      if(NC>1){
        builder=list()
        for(chr in 1:NC){
          builder[[chr]]=build_HMM_matrix_t(k,(rho[[chr]]),beta=beta,L=L[chr],Pop=Pop,Xi_,Beta,scale=window_scaling,sigma =sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
        }
      }
    }
    if(NC==1){
      Q = builder[[1]]
      nu= builder[[2]]
      Tc=builder[[3]]
      g=matrix(0,nrow=length(Tc),ncol=2)
      if(!FS){
        g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
        g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
      }
      if(FS){
        g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
        g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
      }
      
      Tc_r=builder[[4]]
    }
    if(NC>1){
      Q=list()
      nu=list()
      Tc=list()
      Tc_r=list()
      g=list()
      for(chr in 1:NC){
        Q[[chr]] = builder[[chr]][[1]]
        nu[[chr]]= builder[[chr]][[2]]
        Tc[[chr]]=builder[[chr]][[3]]
        Tc_r[[chr]]=builder[[chr]][[4]]
        g[[chr]]=matrix(0,nrow=length(Tc[[chr]]),ncol=2)
        if(!FS){
          g[[chr]][,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]])
          g[[chr]][,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]])          }
        if(FS){
          g[[chr]][,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]]/3))
          g[[chr]][,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc[[chr]]/3))
        }
      }
    }
    res_t$LH<-LH
    res_t$Tc<-Tc_r
    Tc_or=Tc_r
    Tc_o=Tc
    if(!Popfix){
      Xi_t=oldXi_
      Xit=vector()
      pop_vect=get('pop_vect', envir=test.env)
      xx=0
      for(ix in 1:length(Xi_t)){
        x=xx+1
        xx = xx + pop_vect[ix]
        Xit[x:xx]=Xi_t[ix]
      }
      Xit=Xit*sum(BoxP)
      Xit=Xit-(BoxP[1])
      Xit=10^Xit
    }
    if(NC>1){
      rho_t=list()
      rhot=list()
      for(chr in 1:NC){
        rho_t[[chr]]=oldrho[[chr]]*sum(Boxr)
        rho_t[[chr]]=rho_t[[chr]]-(Boxr[1])
        rho_t[[chr]]=10^(rho_t[[chr]])
        rho_t[[chr]]=rho_t[[chr]]*Rho[chr]
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          rhot[[chr]][x:xx]=rho_t[[chr]][ix]
        }
      }
    }
    if(NC==1){
      rho_t=oldrho*sum(Boxr)
      rho_t=rho_t-(Boxr[1])
      rho_t=10^(rho_t)
      rho_t=rho_t*Rho
      xx=0
      rhot=vector()
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        rhot[x:xx]=rho_t[ix]
      }
    }
    if(SB){
      beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
      beta=vector()
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        beta[x:xx]=beta_[ix]
      }
    }
    if(SF){
      sigma_=oldsigma*(Boxs[2]-Boxs[1])
      sigma_=sigma_+Boxs[1]
      sigma=vector()
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        sigma[x:xx]= sigma_[ix]
      }
    }
    if(!Popfix){
      res_t$Xi=Xit
    }
    res_t$mu<-mu
    res_t$beta=beta
    res_t$sigma=sigma
    if(NC==1){
      rhot=rhot/(2*L)
      
    }
    if(NC>1){
      for(chr in 1:NC){
        rhot[[chr]]=rhot[[chr]]/(2*L[chr])
      }
    }
    res_t$rho=rhot
    res_t$Beta<- Beta
    res_t$Self<- Self
    old_list[[length(old_list)+1]]=res_t
    if(!SB){
      oldBeta=Beta
    }
    if(SB){
      oldBeta=Beta
      Beta=beta
    }
    if(SF){
      oldSelf=Self
      Self=sigma
    }
    if(!SF){
      oldSelf=Self
    }
    if(Beta==oldBeta&Self==oldSelf){
      mb=maxBit
      print("no need  for more iteration")
    }
    if(NC==1){
      gamma_t=mean(rhot)/mu
    }
    if(NC>1){
      gamma_t=vector()
      for(chr in 1:NC){
        gamma_t[chr]=mean(rhot[[chr]])/mu[chr]
      }
    }
    print(gamma_t)
    print(gamma)
    if(any(gamma_t<=0.5*gamma|gamma_t>=2*gamma)&redo_R){
      maxBit=maxBit+1
      gamma=gamma_t
    }
    
  }
  if(NC==1){
    
    if(!restart){
      LH=0
      test=Build_zip_Matrix_mailund(Q,g,Os[[1]][[2]],nu)
      for(i in 1:length(Os)){
        fo=forward_zip_mailund(Os[[i]][[1]],g,nu,test[[1]])
        LH=LH+(sum(fo[[2]]))
      }
      res<-list()
      res$LH=LH
      res$Tc=Tc_or
      Q=t(Q)
      #    res$Q <- Q;
      #    res$g <- g;
      if(NC==1){
        rho_=oldrho*sum(Boxr)
        rho_=rho_-(Boxr[1])
        rho_=10^(rho_)
        rho_=rho_*Rho
        rho=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          rho[x:xx]=rho_[ix]/(2*L)
        }
        
      }
      if(SB){
        beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
        beta=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          beta[x:xx]=beta_[ix]
        }
      }
      if(SF){
        sigma_=oldsigma*(Boxs[2]-Boxs[1])
        sigma_=sigma_+Boxs[1]
        sigma=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          sigma[x:xx]= sigma_[ix]
        }
      }
      if(!Popfix){
        res$Xi=Xi_
      }
      #     res$Self<-oldSelf
      if(FS){
        mu=4*mu/3
      }
      res$mu=mu
      res$L<-L
      res$beta=beta
      res$sigma=sigma
      res$rho=rho
      #      res$Beta<- oldBeta
      res$Os<-Os
      #     res$old<-old_list
    }
    
    if(restart){
      res<-old_list[[length(old_list)]]
      res$Os<-Os
      res$old<-old_list
      res$L<-L
    }
  }
  if(NC>1){
    LH=0
    for(chr in 1:NC){
      test=Build_zip_Matrix_mailund(Q[[chr]],g[[chr]],Os[[chr]][[1]][[2]],nu[[chr]])
      for(i in 1:length(Os[[chr]])){
        fo=forward_zip_mailund(Os[[chr]][[i]][[1]],g[[chr]],nu[[chr]],test[[1]])
        LH=LH+(sum(fo[[2]]))
      }
      Q[[chr]]=t(Q[[chr]])
    }
    res <- list();
    res$LH=LH
    res$Tc=Tc_or[[1]]
    #  res$Q <- Q;
    #  res$g <- g;
    res$L<-L
    
    
    rho_=list()
    rho=list()
    for(chr in 1:NC){
      rho_[[chr]]=oldrho[[chr]]*sum(Boxr)
      rho_[[chr]]=rho_[[chr]]-(Boxr[1])
      rho_[[chr]]=10^(rho_[[chr]])
      rho_[[chr]]=rho_[[chr]]*Rho[chr]
      rho[[chr]]=vector()
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        rho[[chr]][x:xx]=rho_[[chr]][ix]/(2*L[chr])
      }
    }
    
    
    if(SB){
      beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
      beta=vector()
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        beta[x:xx]=beta_[ix]
      }
    }
    if(SF){
      sigma_=oldsigma*(Boxs[2]-Boxs[1])
      sigma_=sigma_+Boxs[1]
      sigma=vector()
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        sigma[x:xx]= sigma_[ix]
      }
    }
    if(!Popfix){
      res$Xi=Xi_
    }
    res$beta=beta
    res$sigma=sigma
    res$rho=rho
    # res$Self<-oldSelf
    if(FS){
      mu=4*mu/3
    }
    res$mu=mu
    #  res$Beta<- oldBeta
    res$Os<-Os
    #  res$old<-old_list
  }
  return(res)
}


##############
# Build simu #
##############

Build_simu<-function(NC,n,param){
  O=list()
  for(chr in 1:NC){
    parameters=param[[chr]]
    M=parameters[1]
    theta=parameters[2]
    rho=parameters[3]
    L=parameters[4]
    beta=parameters[5]
    sigma=parameters[6]
    chi=parameters[7:(6+n)]
    time=parameters[(7+n):(6+n+n)]
    setwd("~/Documents/testcpp/scrm_full")
    command_simu=paste("./scrm",M,1,"-t",theta,"-r",rho,L,sep=" ")
    for(x in 1:length(chi)){
      command_simu=paste(command_simu,"-eN",time[x],chi[x],sep=" ")
    }
    command_simu=paste(command_simu," -B 1",beta," -S 1",sigma," -p 9 > Check_simu.txt",sep=" ")
    system(command_simu)
    path="~/Documents/testcpp/scrm_full/Check_simu.txt"
    DNAseqfile=Get_data(path)
    O_total=Seqlist2data(DNAseqfile,L,M,1)
    O[[NC]]=O_total[[1]]
  }
  return(O)
}

################
# Plot results #
################

Plot_esmc_results<-function(results,mu,path=NA,WP=T,NC=1,x=c(1,10^7),y=c(2,8),title="Demographic history"){
  if(is.na(path)){
    path="plot_esmc.pdf"
  }
  if(NC==1){
    if(WP){
      pdf(path)
      plot(x,c(1,1), log=c("x"), ylim =y ,type="n", xlab= paste("Generations ago",sep=" "), ylab="population size (log10)")
      Pop_=results$mu/(2*mu)
      lines((results$Tc*Pop_), log10((as.numeric(results$Xi))*Pop_), type="s", col="red")
      title(title,adj = 0)
      dev.off()
    }
    if(!WP){
      plot(x,c(1,1), log=c("x"), ylim =y ,type="n", xlab= paste("Generations ago",sep=" "), ylab="population size (log10)")
      Pop_=results$mu/(2*mu)
      lines((results$Tc*Pop_), log10((as.numeric(results$Xi))*Pop_), type="s", col="red")
      title(title,adj = 0)
    }
  }
  if(NC>1){
    if(WP){
      pdf(path)
      plot(x,c(1,1), log=c("x"), ylim =y ,type="n", xlab= paste("Generations ago",sep=" "), ylab="population size (log10)")
      for(i in 1:NC){
        Pop_=results[[i]]$mu/(2*mu)
        lines((results[[i]]$Tc*Pop_), log10((as.numeric(results[[i]]$Xi))*Pop_), type="s", col="red")
      }
      title(title,adj = 0)
      dev.off()
    }
    if(!WP){
      plot(x,c(1,1), log=c("x"), ylim =y ,type="n", xlab= paste("Generations ago",sep=" "), ylab="population size (log10)")
      for(i in 1:NC){
        Pop_=results[[i]]$mu/(2*mu)
        lines((results[[i]]$Tc*Pop_), log10((as.numeric(results[[i]]$Xi))*Pop_), type="s", col="red")
      }
      title(title,adj = 0)
    }
  }
}

##################
# Final Function #
##################

get_first_coal_time_msprime<-function(path,mut=T){
  DNAfile=Get_data(path)
  theta=NA
  DNA=matrix(0,ncol=length(DNAfile),nrow=2)
  start=1
  if(mut==T){
    start=2
    theta=as.numeric(substr(DNAfile[[1]][2],6,nchar(DNAfile[[1]][2])))
    if(is.na(theta)){
      browser()
    }
  }
  
  #browser()
  for(i in start:length(DNAfile)){
    data=DNAfile[[i]][1]
    pos_time=which(strsplit(data, "")[[1]]=="e")
    DNA[1,i]=as.numeric(substr(data,pos_time+1,nchar(data)))
    DNA[2,i]=round(as.numeric(substr(data,1,pos_time-5)))
    
  }
  DNA=DNA[,-1]
  pos_0=which(DNA[2,]==0)
  if(length(pos_0)>0){
    DNA=DNA[,-pos_0]
  }
  #diff_check=which(c(DNA[1:(dim(DNA)[2]-1)]-DNA[2:(dim(DNA)[2])])==0)
  
  output=list()
  output$DNA=DNA
  output$theta=theta
  return(output)
}

Get_sim_data<-function(path,L,M,nsim){
  DNAseqfile=Get_data(path)
  O_total=Seqlist2data(DNAseqfile,L,M,nsim)
  return(O_total)
}

Get_real_data<-function(path=NA,M,filename){
  O_total_=create_good_file(path,M,filename)
  O=O_total_$output
  return(O)
}

Plot_beta_sigma<-function(r,rho,BoxB=c(0.05,1),Boxs=c(0,0.99),PLOT=F){
  
  x_s=seq(Boxs[1],Boxs[2],0.01)
  y_b=rho/(r*2*(1-x_s)/(2-x_s))
  keep=which(y_b<=BoxB[2]&y_b>=BoxB[1])
  y_b=y_b[keep]
  x_s=x_s[keep]
  if(PLOT){
    plot(x_s,y_b,type="l",xlab = "Self-fertilization rate",ylab = "germination rate",main = "Possible  germination rate given self-fertilization rate values")
  }
  output=list()
  output$sigma=x_s
  output$beta=y_b
  return(output)
}

get_first_coal_time<-function(file,mut=F){
  data=Get_data(file,heavy=T)
  Output=matrix(0,nrow = 4, ncol= (length(data)))
  start=F
  end=F
  bonus=0
  for(i in 1:length(data)){
    if(length(substr(data[[i]],1,3))>0){
      if(start&mut&substr(data[[i]],1,3)[1]=="seg"){
        end=T
        Output=Output[,-c((i-start_i):(length(data)))]
      }
    }
    if(!end){
      if(start){
        if(length(which(strsplit(data[[i]],split = "")[[1]]=="]"))==1 & length(which(strsplit(data[[i]],split = "")[[1]]=="["))==1){
          Output[4,(i-start_i+bonus)]=as.numeric(substr(data[[i]],2,(gregexpr(data[[i]],pattern = "]")[[1]][1]-1)))
          for(xx in 1:length(as.numeric(gregexpr(data[[i]],pattern = ")")[[1]]))){
            id1=max(as.numeric(gregexpr(data[[i]],pattern = ",")[[1]])[which(as.numeric(gregexpr(data[[i]],pattern = ",")[[1]])<as.numeric(gregexpr(data[[i]],pattern = ")")[[1]][xx]))])
            if(xx==1){
              Output[1,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(gregexpr(data[[i]],pattern = ":")[[1]][min(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])>id1))]-1),(gregexpr(data[[i]],pattern = ":")[[1]][min(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])>id1))]-1)))
              Output[2,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(gregexpr(data[[i]],pattern = ":")[[1]][max(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])<id1))]-1),(gregexpr(data[[i]],pattern = ":")[[1]][max(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])<id1))]-1)))
              Output[3,(i-start_i+bonus)]=2*as.numeric(substr(data[[i]],(gregexpr(data[[i]],pattern = ":")[[1]][max(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])<id1))]+1),(id1-1)))
            }
            if(xx>1){
              if(id1>as.numeric(gregexpr(data[[i]],pattern = ")")[[1]][(xx-1)])){
                genealogy=data[[i]]
                if(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1))!=")"){
                  if(as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]+1),(id1-1)))<(0.5*Output[3,(i-start_i+bonus)])){
                    Output[1,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(gregexpr(data[[i]],pattern = ":")[[1]][min(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])>id1))]-1),(gregexpr(data[[i]],pattern = ":")[[1]][min(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])>id1))]-1)))
                    Output[2,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(gregexpr(data[[i]],pattern = ":")[[1]][max(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])<id1))]-1),(gregexpr(data[[i]],pattern = ":")[[1]][max(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])<id1))]-1)))
                    Output[3,(i-start_i+bonus)]=2*as.numeric(substr(data[[i]],(gregexpr(data[[i]],pattern = ":")[[1]][max(which(as.numeric(gregexpr(data[[i]],pattern = ":")[[1]])<id1))]+1),(id1-1)))
                  }
                }
              }
            }
            
          }
        }else{
          stop("problem in data")
          nb_rec=length(which(strsplit(data[[i]],split = "")[[1]]=="]"))
          pos=which(strsplit(data[[i]],split = "")[[1]]=="]")
          Output=cbind(Output,matrix(0,nrow = 4, ncol= (nb_rec-1)))
          data_s=strsplit(data[[i]],split = "")[[1]]
          for(ii in 1:nb_rec){
            bonus=bonus+1
            Output[4,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(max(which(data_s[1:(pos[ii]-1)]=="["))+1),(pos[ii]-1)))
            Output[1,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(gregexpr(paste(data_s[pos[ii]:length(data_s)],collapse =""),pattern = ":")[[1]][1]-1),(gregexpr(paste(data_s[pos[ii]:length(data_s)],collapse =""),pattern = ":")[[1]][1]-1)))
            Output[2,(i-start_i+bonus)]=as.numeric(substr(data[[i]],(gregexpr(paste(data_s[pos[ii]:length(data_s)],collapse =""),pattern = ":")[[1]][2]-1),(gregexpr(paste(data_s[pos[ii]:length(data_s)],collapse =""),pattern = ":")[[1]][2]-1)))
            Output[3,(i-start_i+bonus)]=2*as.numeric(substr(data[[i]],(gregexpr(paste(data_s[pos[ii]:length(data_s)],collapse =""),pattern = ":")[[1]][1]+1),(gregexpr(paste(data_s[pos[ii]:length(data_s)],collapse =""),pattern = ",")[[1]][1]-1)))
          }
        }
      }
      if(length(data[[i]])>0){
        if(data[[i]]=="//"){
          start=TRUE
          start_i=i
        }
      }
    }
    
  }
  return(Output)
}

get_genealogy<-function(file,M,mut=F){
  data=Get_data(file,heavy = T)
  
  Output=list()
  coal_time=matrix(0,nrow = (M), ncol= (length(data)))
  id_split=matrix(0,nrow = (M-1), ncol= (length(data)))
  id_create=matrix(0,nrow = (M-1), ncol= (length(data)))
  start=F
  end=F
  bonus=0
  for(i in 1:length(data)){
    
    if(length(substr(data[[i]],1,3))>0){
      if(start&mut&substr(data[[i]],1,3)[1]=="seg"){
        end=T
        coal_time=coal_time[,-c((i-start_i):(length(data)))]
        id_split=id_split[,-c((i-start_i):(length(data)))]
        id_create=id_create[,-c((i-start_i):(length(data)))]
      }
    }
    if(!end){
      if(start){
        genealogy=data[[i]]
        if(length(which(strsplit(genealogy,split = "")[[1]]=="]"))==1 & length(which(strsplit(genealogy,split = "")[[1]]=="["))==1){
          coal_time[M,(i-start_i+bonus)]=as.numeric(substr(genealogy,2,(gregexpr(genealogy,pattern = "]")[[1]][1]-1)))
          for(j in 1:(M-1)){
            for(xx in 1:length(as.numeric(gregexpr(genealogy,pattern = ")")[[1]]))){
              id1=max(as.numeric(gregexpr(genealogy,pattern = ",")[[1]])[which(as.numeric(gregexpr(genealogy,pattern = ",")[[1]])<as.numeric(gregexpr(genealogy,pattern = ")")[[1]][xx]))])
              if(xx==1){
                id_split[(M-j),(i-start_i+bonus)]=min(c(as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1))),as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1)))))
                id_create[(M-j),(i-start_i+bonus)]=max(c(as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1))),as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1)))))
                coal_time[(M-j),(i-start_i+bonus)]=2*as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]+1),(id1-1)))
                x_k=xx
                id1_s=id1
              }
              if(xx>1){
                if(id1>as.numeric(gregexpr(genealogy,pattern = ")")[[1]][(xx-1)])){
                  if(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1))!=")"){
                    if(as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]+1),(id1-1)))<(0.5*coal_time[(M-j),(i-start_i+bonus)])){
                      id_split[(M-j),(i-start_i+bonus)]=min(c(as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1))),as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1)))))
                      id_create[(M-j),(i-start_i+bonus)]=max(c(as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][min(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])>id1))]-1))),as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1),(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]-1)))))
                      coal_time[(M-j),(i-start_i+bonus)]=2*as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1))]+1),(id1-1)))
                      x_k=xx
                      id1_s=id1
                    }
                  }
                  
                }
              }
              
            }
            if(j<(M-1)){
              pos_e=c()
              if(x_k<xx){
                pos_e=c(pos_e,as.numeric((gregexpr(genealogy,pattern = ")")[[1]][x_k+1]-1)))
              }
              
              if(!is.na(as.numeric((gregexpr(genealogy,pattern = ",")[[1]][min(na.rm =T,which(as.numeric(gregexpr(genealogy,pattern = ",")[[1]])>as.numeric(gregexpr(genealogy,pattern = ")")[[1]][x_k])))])))){
                pos_e=c(pos_e,(as.numeric((gregexpr(genealogy,pattern = ",")[[1]][min(na.rm =T,which(as.numeric(gregexpr(genealogy,pattern = ",")[[1]])>as.numeric(gregexpr(genealogy,pattern = ")")[[1]][x_k])))]))-1))
              }
              pos_e=min(pos_e)
              time=(0.5*coal_time[(M-j),(i-start_i+bonus)])+as.numeric(substr(genealogy,(gregexpr(genealogy,pattern = ")")[[1]][x_k]+2),pos_e))
              new_ind=paste(id_split[(M-j),(i-start_i+bonus)],":",time,sep="")
              genealogy=paste(substr(genealogy,1,as.numeric(gregexpr(genealogy,pattern = ":")[[1]][max(which(as.numeric(gregexpr(genealogy,pattern = ":")[[1]])<id1_s))]-3)),new_ind,substr(genealogy,(pos_e+1),nchar(genealogy)),sep = "")
            }
            
          }
        }else{
          stop('Problem in data')
        }
        
        
        
      }
      if(length(data[[i]])>0){
        if(data[[i]]=="//"){
          start=TRUE
          start_i=i
        }
      }
    }
    
    
    
  }
  Output$Coal_time=coal_time
  Output$id_split=id_split
  Output$id_create=id_create
  return(Output)
}

get_theta<-function(file){
  data=Get_data(file,heavy=T)
  theta_s=T
  while(theta_s){
    for(i in 1:length(data)){
      if(length(substr(data[[i]],1,3))>0){
        if(substr(data[[i]],1,3)[1]=="seg"){
          theta=as.numeric(strsplit(data[[i]]," ")[2])
          theta_s=F
        }
      }
    }
  }
  return(theta)
}

build_N<-function(genealogy,Tc,Newick=T){
  Output=matrix(0,ncol=length(Tc),nrow=length(Tc))
  if(Newick){
    for(i in 1:dim(genealogy)[2]){
      state=max(which(Tc<genealogy[3,i]))
      if(i==1){
        former_state=state
        length_seq=genealogy[4,i]
        if(length_seq>1){
          Output[former_state,state]=Output[former_state,state]+1
          
          former_state=state
          Output[former_state,state]=Output[former_state,state]+(length_seq-1)
          
        }
      }else{
        length_seq=genealogy[4,i]
        if(length_seq>1){
          Output[former_state,state]=Output[former_state,state]+1
          former_state=state
          Output[former_state,state]=Output[former_state,state]+(length_seq-1)
          
        } else {
          Output[former_state,state]=Output[former_state,state]+1
          former_state=state
        }
      }
      
    }
  }else{
    for(i in 1:dim(genealogy)[2]){
      state=max(which(Tc<genealogy[1,i]))
      if(i==1){
        former_state=state
        length_seq=genealogy[2,i]
        if(length_seq>1){
          Output[former_state,state]=Output[former_state,state]+1
          
          former_state=state
          Output[former_state,state]=Output[former_state,state]+(length_seq-1)
          
        }
      }else{
        length_seq=genealogy[2,i]
        if(length_seq>1){
          Output[former_state,state]=Output[former_state,state]+1
          former_state=state
          Output[former_state,state]=Output[former_state,state]+(length_seq-1)
          
        } else {
          Output[former_state,state]=Output[former_state,state]+1
          former_state=state
        }
      }
      
    }
  }
  
  return(Output)
}

build_M<-function(Mutation,genealogy,Tc){
  Output=matrix(0,ncol=2,nrow=length(Tc))
  L=0
  for(i in 2:dim(genealogy)[2]){
    pos1=1+genealogy[4,(i-1)]+L
    L=L+genealogy[4,i]
    state=max(which(Tc<genealogy[3,i]))
    Output[state,1]=Output[state,1]+genealogy[4,i]
    genealogy[4,i]=pos1
  }
  state=max(which(Tc<genealogy[3,1]))
  q=rep(0,length(Tc))
  q[state]=1
  
  Output[state,1]=Output[state,1]+genealogy[4,1]
  genealogy[4,1]=1
  for(i in 1:dim(Mutation)[2]){
    pos=Mutation[3,i]
    state=max(which(Tc<genealogy[3,max(which(genealogy[4,]<pos))]))
    Output[state,1]=Output[state,1]-1
    Output[state,2]=Output[state,2]+1
  }
  res=list()
  res$M=Output
  res$q=q
  return(res)
}

Optimize_N<-function(file,Rho,theta=NA,L,n=40,ER=T,Pop=T,SB=FALSE,SF=FALSE,BoxB=c(0.1,1),BoxP=c(3,3),Boxr=c(1,1),Boxs=c(0,0.97),pop_vect=NA,window_scaling=c(1,0),sigma=0,beta=1,SCALED=F,Big_Window=F,NC=1,npair=2,FS=F,BW=F,mut=F,Correct_window=F){
  cut_edge=F
  mu=theta/(2*L)
  if(!is.na(theta)){
    gamma=Rho/theta
    print(gamma)
  }
  
  b=get_first_coal_time(file,mut=mut)
  #b[3,]=b[3,]#/2
  
  Popfix=!Pop
  scale_T=1
  if(Correct_window){
    theta_=get_theta(file)
    scale_T=theta_/theta
    b[3,]=b[3,]/scale_T
    mu=theta_/(2*L)
  }
  if(is.na(mu)){
    stop()
  }
  if(as.numeric(Big_Window)==0){
    Vect=0:(n-1)
    Tc= window_scaling[2] -(0.5*(2-sigma)*log(1-(Vect/n))/((beta^2)*window_scaling[1]))
  }
  if(as.numeric(Big_Window)==1){
    alpha_t=0.1
    tmax=15
    Vect=1:(n-1)
    alpha_t=alpha_t/(npair)
    Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
  }
  if(as.numeric(Big_Window)==2){
    tmax=100
    alpha_t=0.1
    Vect=1:(n-1)
    npair=10
    alpha_t=alpha_t/npair
    Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
  }
  #print(scale_T)
  #Tc=Tc*scale_T
  N=build_N(b,Tc)
  for(xxx in 1:dim(N)[1]){
    #N[which(N[,xxx]==0),xxx]<-0.0001
  }
  print("N is built")
  
  # browser()
  if(BW){
    a=Get_sim_data(file,L,2,1)
    res=build_M(a[[1]],b,Tc)
    print("M is built")
    M=res$M
    #M[which(M[,1]==0),1]<-0.0001
    #M[which(M[,2]==0),2]<-0.0001
    print(sum(M))
    q_=res$q
  }
  if(SCALED){
    corrector_N=rowSums(N)
    N=diag(1/corrector_N)%*%N
    if(BW){
      corrector_M=rowSums(M)
      M=diag(1/corrector_M)%*%M
    }
  }
  #Tc=Tc/scale_T
  scale_T=1
  Rho=gamma*2*L*mu
  print(Rho/(2*L))
  test.env <- new.env()
  test.env$L <- L
  test.env$k <- n
  test.env$Rho <- Rho
  test.env$window_scaling <- window_scaling
  test.env$Pop<-!Pop
  test.env$NC<-NC
  test.env$FS<-FS
  test.env$Big_Window <- Big_Window
  test.env$npair <- npair
  test.env$Beta <- beta
  test.env$Self <- sigma
  test.env$BoxB <- BoxB
  test.env$Boxs <- Boxs
  test.env$Big_Xi <- N
  test.env$mu_b <- mu
  test.env$BW <- BW
  test.env$scale_T <- scale_T
  if(BW){
    test.env$Big_M <- M
    test.env$q_ <- q_
  }
  lr=length(Rho)
  test.env$lr<-lr
  if(any(!is.na(pop_vect))){
    Klink=length(pop_vect)
  }
  if((all(is.na(pop_vect))|sum(pop_vect)!=n)){
    Klink=0.5*n
    pop_vect=rep(2, Klink)
    print("Default pop vector")
  }
  test.env$pop_vect <- pop_vect
  if(SB){
    oldbeta=(sqrt(beta)-BoxB[1])/(BoxB[2]-BoxB[1])
  }
  if(SF){
    oldsigma=(sigma-Boxs[1])/(Boxs[2]-Boxs[1])
  }
  if(!SB){
    test.env$beta <- beta
  }
  if(!SF){
    test.env$sigma <- sigma
  }
  if(Pop){
    oldXi_=rep((BoxP[1]/sum(BoxP)),Klink)
    oldXi=vector()
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      oldXi[x:xx]=oldXi_[ix]
    }
    Xi_=oldXi*sum(BoxP)
    Xi_=Xi_-(BoxP[1])
    Xi_=10^Xi_
  }
  
  if(!ER){
    oldrho=0
    if(NC>1){
      oldrho=rep(0,NC)
    }
    
  }
  if(ER){
    oldrho=(Boxr[1]/sum(Boxr))
    if(NC>1){
      oldrho=rep((Boxr[1]/sum(Boxr)),NC)
    }
    
  }
  
  if(NC==1){
    if(ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              FS=get('FS', envir=test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[3]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              builder=build_HMM_matrix(n,rho_,beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge) #
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              
              
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:3,1])
            rho=sol[1]
            beta_=sol[2]
            sigma_=sol[3]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_)))
            oldrho=rho
            oldbeta=beta_
            oldsigma=sigma_
            
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[3]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[4:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,rho_,beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(3+(Klink)),1])
            rho=sol[1]
            beta_=sol[2]
            sigma_=sol[3]
            Xi_=sol[4:length(sol)]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldrho=rho
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma',envir = test.env)
              Self=get('Self',envir = test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,rho_,beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:2,1])
            rho=sol[1]
            beta_=sol[2]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho,oldbeta-beta_)))
            oldrho=rho
            oldbeta=beta_
            
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma',envir = test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[2]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BETA=get('BETA',envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[3:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2+(Klink)),1])
            rho=sol[1]
            beta_=sol[2]
            Xi_=sol[3:length(sol)]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,Xi_-oldXi_)))
            oldrho=rho
            oldbeta=beta_
            oldXi_=Xi_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldsigma), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:2,1])
            rho=sol[1]
            sigma_=sol[2]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho,oldsigma-sigma_)))
            oldrho=rho
            oldsigma=sigma_
            
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[3:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2+(Klink)),1])
            rho=sol[1]
            sigma_=sol[2]
            Xi_=sol[3:length(sol)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_,oldsigma-sigma_)))
            oldrho=rho
            oldXi_=Xi_
            oldsigma=sigma_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma', envir=test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldrho),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1,1])
            rho=sol[1]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho)))
            oldrho=rho
            
          }
          
          if(!Popfix){
            
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma', envir=test.env)
              rho_=param[1]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[2:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+(Klink)),1])
            rho=sol[1]
            Xi_=sol[2:length(sol)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
            oldrho=rho
            oldXi_=Xi_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
    }
    if(!ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:2,1])
            beta_=sol[1]
            sigma_=sol[2]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_)))
            oldbeta=beta_
            oldsigma=sigma_
            
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[3:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2+(Klink)),1])
            beta_=sol[1]
            sigma_=sol[2]
            Xi_=sol[3:length(sol)]
            diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              sigma=get('sigma', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1,1])
            beta_=sol[1]
            print(paste(" new Complete likelihood : ", LH ))
            
            diff=max(abs(c(oldbeta-beta_)))
            oldbeta=beta_
            
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              sigma=get('sigma', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[2:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge) #
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+(Klink)),1])
            beta_=sol[1]
            Xi_=sol[2:length(sol)]
            
            diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1,1])
            sigma_=sol[1]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(oldsigma-sigma_)))
            oldsigma=sigma_
            
          }
          
          if(!Popfix){
            
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              Rho=get('Rho', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[2:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(Big_M)>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              
              return(LH)
            }
            sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+(Klink)),1])
            sigma_=sol[1]
            Xi_=sol[2:length(sol)]
            diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
            oldXi_=Xi_
            oldsigma=sigma_
            print(sigma_)
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[1:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix(n,(rho_),beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              
              scale_T=get('scale_T', envir=test.env)
              Tc=Tc*scale_T
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              if(BW){
                g=matrix(0,nrow=length(Tc),ncol=2)
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              return(LH)
            }
            sol= BBoptim(c(oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((Klink)),1])
            Xi_=sol[1:length(sol)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
            oldXi_=Xi_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        
        
      }
    }
  }
  if(NC>1){
    if(ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[lr+2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              LH=0
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(length(oldrho)+2),1])
            rho=sol[1:length(oldrho)]
            beta_=sol[length(oldrho)+1]
            sigma_=sol[length(oldrho)+2]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_)))
            oldrho=rho
            oldbeta=beta_
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[lr+2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[3+lr:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2+length(oldrho)+(Klink)),1])
            rho=sol[1:length(oldrho)]
            beta_=sol[length(oldrho)+1]
            sigma_=sol[length(oldrho)+2]
            Xi_=sol[(length(oldrho)+3):length(sol)]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldrho=rho
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              LH=0
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(length(oldrho)+1),1])
            rho=sol[1:length(oldrho)]
            beta_=sol[length(oldrho)+1]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_)))
            oldrho=rho
            oldbeta=beta_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(lr+1)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(lr+2):length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              LH=0
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for( chr in 1:chr){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+length(oldrho)+(Klink)),1])
            rho=sol[1:length(oldrho)]
            beta_=sol[1+length(oldrho)]
            Xi_=sol[(2+length(oldrho)):length(sol)]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,Xi_-oldXi_)))
            oldrho=rho
            oldbeta=beta_
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
          }
        }
      }
      
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              sigma=param[1+lr]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              LH=0
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(length(oldrho)+1),1])
            rho=sol[1:length(oldrho)]
            sigma_=sol[length(oldrho)+1]
            diff=max(abs(c(rho- oldrho,oldsigma-sigma_)))
            oldrho=rho
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(2+lr):length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1+lr]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              
              LH=0
              for( chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+length(oldrho)+(Klink)),1])
            rho=sol[1:length(oldrho)]
            sigma_=sol[1+length(oldrho)]
            Xi_=sol[(2+length(oldrho)):length(sol)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_,oldsigma-sigma_)))
            oldrho=rho
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              LH=0
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for( chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(length(oldrho)),1])
            rho=sol[1:length(oldrho)]
            diff=max(abs(c(rho- oldrho)))
            oldrho=rho
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              lr=get('lr', envir=test.env)
              rho_=param[1:lr]
              rho_=rho_*sum(Boxr)
              rho_=rho_-(Boxr[1])
              rho_=10^(rho_)
              Rho=get('Rho', envir=test.env)
              rho_=rho_*Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(1+lr):length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              LH=0
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for( chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(length(oldrho)+(Klink)),1])
            rho=sol[1:length(oldrho)]
            Xi_=sol[(1+length(oldrho)):length(sol)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
            oldrho=rho
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
          }
        }
      }
    }
    if(!ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              
              LH=0
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:2,1])
            beta_=sol[1]
            sigma_=sol[2]
            diff=max(abs(c(oldbeta- beta,oldsigma-sigma)))
            oldbeta=sol[1]
            oldsigma=sol[2]
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC',envir = test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[2]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[3:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2+(Klink)),1])
            beta_=sol[1]
            sigma_=sol[2]
            Xi_=sol[3:length(sol)]
            diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir = test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              LH=0
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1,1])
            beta_=sol[1]
            diff=max(abs(c(oldbeta-beta_)))
            oldbeta=beta_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC',envir = test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[2:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+(Klink)),1])
            beta_=sol[1]
            Xi_=sol[2:length(sol)]
            diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
          }
        }
      }
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldsigma),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1,1])
            sigma_=sol[1]
            diff=max(abs(c(oldsigma-sigma_)))
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              Rho=get('Rho', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir = test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[2:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(1+(Klink)),1])
            sigma_=sol[1]
            Xi_=sol[2:length(sol)]
            diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
            oldXi_=Xi_
            oldsigma=sigma_
            print(sigma_)
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC',envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[1:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              BW=get('BW', envir=test.env)
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                q_=get('q_', envir=test.env)
              }
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                scale_T=get('scale_T', envir=test.env)
                Tc=Tc*scale_T
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                if(BW){
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  x=as.vector(g)
                  keep=which(x>0&as.vector(as.vector(Big_M[[chr]]))>0)
                  x=x[keep]
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  nu=builder[[2]]
                  LH=-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((Klink)),1])
            Xi_=sol[1:length(sol)]
            diff=max(abs(Xi_-oldXi_))
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
            
          }
        }
      }
    }
  }
  
  if(SB){
    beta=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
  }
  if(SF){
    sigma=oldsigma*(Boxs[2]-Boxs[1])
    sigma=sigma+Boxs[1]
  }
  if(ER){
    rho_=oldrho*sum(Boxr)
    rho_=rho_-(Boxr[1])
    rho_=10^(rho_)
    rho_=rho_*Rho
    rho_=rho_/(2*L)
  }else{
    rho_=Rho/(2*L)
  }
  if(Pop){
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      oldXi[x:xx]=oldXi_[ix]
    }
    Xi_=oldXi*sum(BoxP)
    Xi_=Xi_-(BoxP[1])
    Xi_=10^Xi_
  }
  res<-list()
  res$beta=beta
  res$sigma=sigma
  res$rho=rho_
  if(Pop){
    res$Xi=Xi_
  }
  Beta=get('Beta', envir=test.env)
  Self=get('Self', envir=test.env)
  res$scale_T=scale_T
  res$mu=mu
  res$Tc=Tc
  res$N=N
  return(res)
}

get_Mat<-function(file,Rho,theta=NA,L,n=40,BoxB=c(0.1,1),BoxP=c(3,3),Boxr=c(1,1),Boxs=c(0,0.97),pop_vect=NA,window_scaling=c(1,0),sigma=0,beta=1,SCALED=F,Big_Window=F,NC=1,npair=2,mut=F,Correct_window=F,build_M=F){
  output=list()
  cut_edge=F
  mu=theta/(2*L)
  b=get_first_coal_time(file,mut=mut)
  b[3,]=b[3,]#/2
  Popfix=!Pop
  scale_T=1
  if(Correct_window){
    theta_=get_theta(file)
    scale_T=theta_/theta
    mu=theta_/(2*L)
  }
  if(is.na(mu)){
    stop()
  }
  if(as.numeric(Big_Window)==0){
    Vect=0:(n-1)
    Tc= window_scaling[2] -(0.5*(2-sigma)*log(1-(Vect/n))/((beta^2)*window_scaling[1]))
  }
  if(as.numeric(Big_Window)==1){
    alpha_t=0.1
    tmax=15
    Vect=1:(n-1)
    alpha_t=alpha_t/(npair)
    Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
  }
  if(as.numeric(Big_Window)==2){
    tmax=100
    alpha_t=0.1
    Vect=1:(n-1)
    npair=10
    alpha_t=alpha_t/npair
    Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
  }
  print(scale_T)
  Tc=Tc*scale_T
  N=build_N(b,Tc)
  for(xxx in 1:dim(N)[1]){
    #N[which(N[,xxx]==0),xxx]<-0.0001
  }
  print("N is built")
  if(build_M){
    a=Get_sim_data(file,L,2,1)
    res=build_M(a[[1]],b,Tc)
    print("M is built")
    M=res$M
    #M[which(M[,1]==0),1]<-0.0001
    #M[which(M[,2]==0),2]<-0.0001
    print(sum(M))
    q_=res$q
  }
  if(SCALED){
    corrector_N=rowSums(N)
    N=diag(1/corrector_N)%*%N
    if(build_M){
      corrector_M=rowSums(M)
      M=diag(1/corrector_M)%*%M
    }
  }
  output$N=N
  if(build_M){
    output$M=M
    output$q=q_
  }else{
    output$M=0
    output$q=0
  }
  return(output)
  
}

Optimize_N_t<-function(file,Rho,theta=NA,L,n=30,ER=T,Pop=T,SB=FALSE,SF=FALSE,BoxB=c(0.1,1),BoxP=c(3,3),Boxr=c(1,1),Boxs=c(0,0.97),pop_vect=NA,window_scaling=c(1,0),sigma=0,beta=1,SCALED=F,Big_Window=F,NC=1,npair=2,FS=F,BW=F,mut=F,Correct_window=F,Newick=T){
  cut_edge=F
  Share=T
  if(!is.na(theta)){
    gamma=theta/Rho
  }
  if(NC==1){
    mu=theta/(2*L)
    if(Newick){
      
      b=get_first_coal_time(file,mut=mut)
      b[3,]=b[3,]#/2
      Popfix=!Pop
      scale_T=1
      if(Correct_window){
        theta_=get_theta(file)
        scale_T=theta_/theta
        mu=theta_/(2*L)
      }
      if(is.na(mu)){
        stop()
      }
      if(as.numeric(Big_Window)==0){
        Vect=0:(n-1)
        Tc= window_scaling[2] -(0.5*(2-sigma)*log(1-(Vect/n))/((beta^2)*window_scaling[1]))
      }
      if(as.numeric(Big_Window)==1){
        alpha_t=0.1
        tmax=15
        Vect=1:(n-1)
        alpha_t=alpha_t/(npair)
        Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
      }
      if(as.numeric(Big_Window)==2){
        tmax=100
        alpha_t=0.1
        Vect=1:(n-1)
        npair=10
        alpha_t=alpha_t/npair
        Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
      }
      #print(scale_T)
      Tc=Tc*scale_T
      N=build_N(b,Tc)
    }else{
      b=get_first_coal_time_msprime(file,mut=mut)
      #b[1,]=b[1,]#/2
      Popfix=!Pop
      scale_T=1
      if(Correct_window){
        theta_=b$theta
        scale_T=theta_/theta
        mu=theta_/(2*L)
      }
      if(is.na(mu)){
        stop()
      }
      if(as.numeric(Big_Window)==0){
        Vect=0:(n-1)
        Tc= window_scaling[2] -(0.5*(2-sigma)*log(1-(Vect/n))/((beta^2)*window_scaling[1]))
      }
      if(as.numeric(Big_Window)==1){
        alpha_t=0.1
        tmax=15
        Vect=1:(n-1)
        alpha_t=alpha_t/(npair)
        Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
      }
      if(as.numeric(Big_Window)==2){
        tmax=100
        alpha_t=0.1
        Vect=1:(n-1)
        npair=10
        alpha_t=alpha_t/npair
        Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
      }
      #print(scale_T)
      Tc=Tc*scale_T
      N=build_N(b$DNA,Tc,Newick)
    }
  }else{
    if(Share){
      mu=theta/(2*L)
      N=matrix(0,n,n)
      mu_t=c()
      for(chr in 1:NC){
        if(Newick){
          
          b=get_first_coal_time(file[[chr]],mut=mut)
          b[3,]=b[3,]#/2
          Popfix=!Pop
          scale_T=1
          if(Correct_window){
            theta_=get_theta(file)
            scale_T=theta_/theta
            mu_t=c(mu_t,theta_/(2*L))
          }
          if(is.na(mu)){
            stop()
          }
          if(as.numeric(Big_Window)==0){
            Vect=0:(n-1)
            Tc= window_scaling[2] -(0.5*(2-sigma)*log(1-(Vect/n))/((beta^2)*window_scaling[1]))
          }
          if(as.numeric(Big_Window)==1){
            alpha_t=0.1
            tmax=15
            Vect=1:(n-1)
            alpha_t=alpha_t/(npair)
            Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
          }
          if(as.numeric(Big_Window)==2){
            tmax=100
            alpha_t=0.1
            Vect=1:(n-1)
            npair=10
            alpha_t=alpha_t/npair
            Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
          }
          #print(scale_T)
          Tc=Tc*scale_T
          N=N+build_N(b,Tc)
        }else{
          b=get_first_coal_time_msprime(file[[chr]],mut=mut)
          b$DNA[1,]=b$DNA[1,]/2
          Popfix=!Pop
          scale_T=1
          if(Correct_window){
            theta_=b$theta
            #scale_T=theta_/theta
            mu_t=c(mu_t,theta_/(2*L))
          }
          if(is.na(mu)){
            stop()
          }
          if(as.numeric(Big_Window)==0){
            Vect=0:(n-1)
            Tc= window_scaling[2] -(0.5*(2-sigma)*log(1-(Vect/n))/((beta^2)*window_scaling[1]))
          }
          if(as.numeric(Big_Window)==1){
            alpha_t=0.1
            tmax=15
            Vect=1:(n-1)
            alpha_t=alpha_t/(npair)
            Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
          }
          if(as.numeric(Big_Window)==2){
            tmax=100
            alpha_t=0.1
            Vect=1:(n-1)
            npair=10
            alpha_t=alpha_t/npair
            Tc= c(0,window_scaling[2] + (0.5*(2-sigma)*(alpha_t*exp((Vect/n)*log(1+(tmax/alpha_t))-1))/((beta^2)*window_scaling[1])))
          }
          print(scale_T)
          scale_T=mean(b$DNA[1,])/(mean(Tc))
          Tc=Tc*scale_T
          print(scale_T)
          N=N+build_N(b$DNA,Tc,Newick)
        }
      }
      if(Correct_window){
        mu=mean(mu_t)}
      if(is.na(mu)){
        browser()
      }
      rm(mu_t)
      NC=1
    }else{
    }
  }
  
  
  
  print("N is built")
  
  # browser()
  if(BW){
    a=Get_sim_data(file,L,2,1)
    res=build_M(a[[1]],b,Tc)
    print("M is built")
    M=res$M
    #M[which(M[,1]==0),1]<-0.0001
    #M[which(M[,2]==0),2]<-0.0001
    print(sum(M))
    q_=res$q
  }
  if(SCALED){
    corrector_N=rowSums(N)
    N=diag(1/corrector_N)%*%N
    if(BW){
      corrector_M=rowSums(M)
      M=diag(1/corrector_M)%*%M
    }
  }
  Tc=Tc/scale_T
  scale_T=1
  Rho=(2*L)*mu*gamma
  test.env <- new.env()
  test.env$L <- L
  test.env$k <- n
  test.env$Rho <- Rho
  test.env$window_scaling <- window_scaling
  test.env$Pop<-!Pop
  test.env$NC<-NC
  test.env$FS<-FS
  test.env$Big_Window <- Big_Window
  test.env$npair <- npair
  
  
  test.env$Self <- sigma
  test.env$BoxB <- BoxB
  test.env$Boxs <- Boxs
  test.env$Big_Xi <- N
  test.env$mu_b <- mu
  test.env$BW <- BW
  if(BW){
    test.env$Big_M <- M
    test.env$q_ <- q_
  }
  lr=length(Rho)
  test.env$lr<-lr
  if(any(!is.na(pop_vect))){
    Klink=length(pop_vect)
  }
  if((all(is.na(pop_vect))|sum(pop_vect)!=n)){
    Klink=0.5*n
    pop_vect=rep(2, Klink)
    print("Default pop vector")
  }
  test.env$pop_vect <- pop_vect
  if(SB){
    BoxB[1]=max(sqrt(0.01),sqrt(BoxB[1]))
    BoxB[2]=min(sqrt(1),sqrt(BoxB[2]))
    beta=max((BoxB[1]^2),beta)
    beta=min(beta,(BoxB[2]^2))
    Beta=beta
    oldbeta=rep(beta,length(pop_vect))
    
  }
  if(SF){
    sigma=min(Boxs[2],sigma)
    sigma=max(sigma,Boxs[1])
    Self=sigma
    sigma=rep(sigma,length(pop_vect))
    oldsigma=(sigma-Boxs[1])/(Boxs[2]-Boxs[1])
  }
  if(Pop){
    oldXi_=rep((BoxP[1]/sum(BoxP)),Klink)
    oldXi=vector()
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      oldXi[x:xx]=oldXi_[ix]
    }
    Xi_=oldXi*sum(BoxP)
    Xi_=Xi_-(BoxP[1])
    Xi_=10^Xi_
  }
  if(!SB){
    Beta=beta
    test.env$beta <- beta
  }
  if(!SF){
    Self=sigma
    oldSelf=Self
    test.env$sigma <- sigma
  }
  if(!ER){
    Boxr=c(0,0)
  }
  if(ER){
    oldrho=rep((Boxr[1]/sum(Boxr)),length(pop_vect))
  }
  if(NC==1){
    if(ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(1+(2*Klink)):(3*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              Pop=get('Pop', envir=test.env)
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              beta_=vector()
              rho_=vector()
              sigma_=vector()
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                beta_[x:xx]=beta[ix]
                rho_[x:xx]=rho[ix]
                sigma_[x:xx]=sigma[ix]
              }
              
              builder=build_HMM_matrix_t(n,rho_,beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair) #
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:3*Klink])
            rho=sol[1:Klink]
            beta_=sol[(1+Klink):(2*Klink)]
            sigma_=sol[(1+(2*Klink)):(3*Klink)]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_)))
            oldrho=rho
            oldbeta=beta_
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[(1+(2*Klink)):(3*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[(1+(3*Klink)):(4*Klink)]
              rho_=vector()
              beta_=vector()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                rho_=rho[ix]
                beta_=beta[ix]
                sigma_=sigma[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix_t(n,rho_,beta=beta_,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(4*Klink),1])
            rho=sol[1:Klink]
            beta_=sol[(1+Klink):(2*Klink)]
            sigma_=sol[(1+(2*Klink)):(3*Klink)]
            Xi_=sol[(1+(3*Klink)):(4*Klink)]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldrho=rho
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma',envir = test.env)
              Self=get('Self',envir = test.env)
              Klink=get('Klink',envir = test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              M=get('M', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              rho_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                rho_=rho[ix]
                beta_=beta[ix]
              }
              builder=build_HMM_matrix_t(n,rho_,beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*KLink),1])
            rho=sol[1:Klink]
            beta_=sol[(1+Klink):(2*Klink)]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho,oldbeta-beta_)))
            oldrho=rho
            oldbeta=beta_
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma',envir = test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(1+Klink):(2*Klink)]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BETA=get('BETA',envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(1+(2*Klink)):(3*Klink)]
              Xi=vector()
              rho_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                rho_[x:xx]=rho[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldbeta,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(3*(Klink)),1])
            rho=sol[1:Klink]
            beta_=sol[(1+(1*Klink)):(2*Klink)]
            Xi_=sol[(1+(2*Klink)):(3*Klink)]
            diff=max(abs(c(rho- oldrho,oldbeta-beta_,Xi_-oldXi_)))
            oldrho=rho
            oldbeta=beta_
            oldXi_=Xi_
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(1+Klink):(2*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              sigma_=vector()
              rho_=vector()
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                
                rho_[x:xx]=rho[ix]
                sigma_[x:xx]=sigma[ix]
              }
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldsigma), function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*Klink),1])
            rho=sol[1:Klink]
            sigma_=sol[(1+Klink):(2*Klink)]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho,oldsigma-sigma_)))
            oldrho=rho
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              Boxs=get('Boxs', envir=test.env)
              sigma=param[(1+(1*Klink)):(2*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(1+(2*Klink)):(3*Klink)]
              Xi=vector()
              sigma_=vector()
              rho_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                sigma_[x:xx]=sigma[ix]
                rho_[x:xx]=rho[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              Self=get('Self', envir=test.env)
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(3*(Klink)),1])
            rho=sol[1:Klink]
            sigma_=sol[(1+(1*Klink)):(2*Klink)]
            Xi_=sol[(1+(2*Klink)):(3*Klink)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_,oldsigma-sigma_)))
            oldrho=rho
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              rho_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                rho_[x:xx]=rho[ix]
              }
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:Klink,1])
            rho=sol[1:Klink]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(rho- oldrho)))
            oldrho=rho
          }
          
          if(!Popfix){
            
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Klink=get('Klink', envir=test.env)
              sigma=get('sigma', envir=test.env)
              rho=param[1:Klink]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(1+Klink):(2*Klink)]
              Xi=vector()
              rho_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:length(Xi_)){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                rho_[x:xx]=rho[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,Sigma=Self,scale=window_scaling,FS=FS,sigma = sigma,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              
              if(BW){
                
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldrho,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldrho,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*Klink),1])
            rho=sol[1:Klink]
            Xi_=sol[(1+Klink):(2*Klink)]
            diff=max(abs(c(rho- oldrho,Xi_-oldXi_)))
            oldrho=rho
            oldXi_=Xi_
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
    }
    if(!ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(1+Klink):(2*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              beta_=vector()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
                beta_[x:xx]=beta[ix]
              }
              builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*Klink),1])
            beta_=sol[1:Klink]
            sigma_=sol[(1+Klink):(2*Klink)]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_)))
            oldbeta=beta_
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Klink=get('Klink', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[(1+Klink):(2*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[(1+(2*Klink)):(3*Klink)]
              Xi=vector()
              sigma_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                sigma_[x:xx]=sigma[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(3*(Klink)),1])
            beta_=sol[1:Klink]
            sigma_=sol[(1+(1*Klink)):(2*Klink)]
            Xi_=sol[(1+(2*Klink)):(3*Klink)]
            diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              sigma=get('sigma', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              M=get('M', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                beta_[x:xx]=beta[ix]
              }
              builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:Klink,1])
            beta_=sol[1:Klink]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(oldbeta-beta_)))
            oldbeta=beta_
            
          }
          if(!Popfix){
            function_to_minimize<-function(param){
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              sigma=get('sigma', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Klink=get('Klink', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(1+Klink):(2*Klink)]
              Xi=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta_,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair) #
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldbeta,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*(Klink)),1])
            beta_=sol[1:Klink]
            Xi_=sol[(1+Klink):(2*Klink)]
            diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Klink=get('Klink', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1:Klink]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
              }
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:Klink,1])
            sigma_=sol[1:Klink]
            print(paste(" new Complete likelihood : ", LH ))
            diff=max(abs(c(oldsigma-sigma_)))
            oldsigma=sigma_
            
          }
          
          if(!Popfix){
            
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              Rho=get('Rho', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              Klink=get('Klink', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1:Klink]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[(1+Klink):(2*Klink)]
              Xi=vector()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                sigma_[x:xx]=sigma[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldsigma,oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*(Klink)),1])
            sigma_=sol[1:Klink]
            Xi_=sol[(1+Klink):(2*Klink)]
            
            diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
            oldXi_=Xi_
            oldsigma=sigma_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(!Popfix){
            function_to_minimize<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              sigma=get('sigma', envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[1:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              Self=get('Self', envir=test.env)
              
              builder=build_HMM_matrix_t(n,(rho_),beta=beta,Pop = Pop,Xi=Xi,L=L,Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
              Q=builder[[1]]
              Q=t(Q)
              A=as.vector(Q)
              keep=which(A>0&as.vector(Big_Xi)>0)
              A=A[keep]
              Big_Xi=as.vector(Big_Xi)
              Big_Xi=Big_Xi[keep]
              
              Tc=builder[[3]]
              g=matrix(0,nrow=length(Tc),ncol=2)
              if(!FS){
                g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
              }
              if(FS){
                g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
              }
              
              if(BW){
                Big_M=get('Big_M', envir=test.env)
                x=as.vector(g)
                keep=which(x>0&as.vector(as.vector(Big_M))>0)
                x=x[keep]
                
                m=as.vector(Big_M)
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                LH=-sum(log(A)*Big_Xi)-sum(log(x)*m)-sum(log(nu)*q_)
              }
              if(!BW){
                LH=-sum(log(A)*Big_Xi)
              }
              return(LH)
            }
            sol= BBoptim(c(oldXi_),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=15+length(c(oldXi_)),M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((Klink)),1])
            Xi_=sol[1:length(sol)]
            
            oldXi_=Xi_
            
            print(paste(" new Complete likelihood : ",  -LH ))
            print(paste("Xi:",oldXi_))
          }
        }
      }
    }
  }
  if(NC>1){
    if(ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              beta_=vector()
              sigma_=vector()
              rho_=list()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
                beta_[x:xx]=beta[ix]
                
              }
              LH=0
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                
                
                builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_Xi=get('Big_Xi', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
              }
              return(LH)
            }
            sol=BBoptim(c(unlist(oldrho),oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((NC+2)*Klink),1])
            rho=sol[1:(NC*Klink)]
            beta_=sol[(NC*Klink+1):((NC+1)*Klink)]
            sigma_=sol[((NC+1)*Klink+1):((NC+2)*Klink)]
            diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_,oldsigma-sigma_)))
            
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            
            oldbeta=beta_
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              Xi_=param[(Klink*(Nc+2)+1):(Klink*(Nc+3))]
              Xi=vector()
              sigma_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                sigma_[x:xx]=sigma[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
              }
              return(LH)
            }
            sol= BBoptim(c(unlist(oldrho),oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.numeric(sol[1:((NC+2)*Klink),1])
            rho=sol[1:(NC*Klink)]
            beta_=sol[(NC*Klink+1):((NC+1)*Klink)]
            sigma_=sol[((NC+1)*Klink+1):((NC+2)*Klink)]
            Xi_=sol[((NC+2)*Klink+1):((NC+3)*Klink)]
            diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              sigma=get('sigma', envir=test.env)
              Self=get('Self', envir=test.env)
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              Xi_=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
              Xi=vector()
              sigma_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                
                builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])}
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(unlist(oldrho),oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((NC+2)*Klink),1])
            rho=sol[1:((NC)*Klink)]
            beta_=sol[((NC)*Klink+1):((NC+1)*Klink)]
            diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_)))
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            oldbeta=beta_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              beta=((param[(Klink*Nc+1):(Klink*(Nc+1))]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              Xi_=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
              Xi=vector()
              sigma_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(unlist(oldrho),oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(Klink*(NC+2)),1])
            rho=sol[1:(Klink*NC)]
            beta_=sol[(Klink*NC+1):(NC+1)*Klink]
            Xi_=sol[(Klink*(NC+1)+1):(NC+2)*Klink]
            diff=max(abs(c(rho-unlist(oldrho),oldbeta-beta_,Xi_-oldXi_)))
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            oldbeta=beta_
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
          }
        }
      }
      
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              beta=get('beta', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(Klink*(Nc)+1):(Klink*(Nc+1))]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
              }
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                
                builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(unlist(oldrho),oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((NC+1)*Klink),1])
            rho=sol[1:(NC*Klink)]
            sigma_=sol[((NC)*Klink+1):((NC+1)*Klink)]
            diff=max(abs(c(rho-unlist(oldrho),oldsigma-sigma_)))
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              
              beta=get('beta', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(Klink*(Nc)+1):(Klink*(Nc+1))]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              Xi_=param[(Klink*(Nc+1)+1):(Klink*(Nc+2))]
              Xi=vector()
              sigma_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                sigma_[x:xx]=sigma[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(unlist(oldrho),oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((NC+2)*(Klink)),1])
            rho=sol[1:(NC*Klink)]
            sigma_=sol[(Klink*NC+1):((NC+1)*Klink)]
            Xi_=sol[(Klink*(NC+1)+1):((NC+2)*Klink)]
            diff=max(abs(c(rho-unlist(oldrho),Xi_-oldXi_,oldsigma-sigma_)))
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              beta=get('beta', envir=test.env)
              sigma=get('sigma', envir=test.env)
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              pop_vect=get('pop_vect', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                for( chr in 1:NC){
                  builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                  Q=builder[[1]]
                  Q=t(Q)
                  A=as.vector(Q)
                  keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                  A=A[keep]
                  Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                  Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                  Tc=builder[[3]]
                  g=matrix(0,nrow=length(Tc),ncol=2)
                  if(!FS){
                    g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                    g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  }
                  if(FS){
                    g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                    g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  }
                  x=as.vector(g)
                  keep=which(x>0)
                  x=x[keep]
                  
                  if(BW){
                    Big_M=get('Big_M', envir=test.env)
                    m=as.vector(Big_M[[chr]])
                    m=m[keep]
                    q_=get('q_', envir=test.env)
                    nu=builder[[2]]
                    LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                  }
                  if(!BW){
                    LH=LH-sum(log(A)*Big_Xi[[chr]])
                  }
                  
                }
                return(LH)
              }
              sol= BBoptim(c(unlist(oldrho)),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
              LH=as.numeric(as.matrix(sol[[2]]))
              sol=as.matrix(sol[[1]])
              sol=as.numeric(sol[1:(NC*Klink),1])
              rho=sol[1:(NC*Klink)]
              diff=max(abs(c(rho- unlist(oldrho))))
              for(chr in 1:NC){
                xx=1+(chr-1)*Klink
                yy=chr*Klink
                oldrho[[chr]]=rho[xx:yy]
              }
            }
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              
              beta=get('beta', envir=test.env)
              sigma=get('sigma', envir=test.env)
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC', envir=test.env)
              Klink=get('Klink', envir=test.env)
              rho=param[1:(Klink*Nc)]
              rho=rho*sum(Boxr)
              rho=rho-(Boxr[1])
              rho=10^(rho)
              Rho=get('Rho', envir=test.env)
              rho=rho*Rho
              BoxB=get('BoxB', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              rho_=list()
              Xi_=param[(Klink*(Nc)+1):(Klink*(Nc+1))]
              Xi=vector()
              sigma_=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                rho_[[chr]]=vector()
                xx=0
                for(ix in 1:Klink){
                  x=xx+1
                  xx = xx + pop_vect[ix]
                  xr=(chr-1)*Klink+ix
                  rho_[[chr]][x:xx]=rho[xr]
                  
                }
                builder=build_HMM_matrix_t(n,(rho_[[chr]]),beta=beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(unlist(oldrho),oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((NC+1)*(Klink)),1])
            rho=sol[1:(NC*Klink)]
            Xi_=sol[(NC*Klink+1):((NC+1)*(Klink))]
            diff=max(abs(c(rho- unlist(oldrho),Xi_-oldXi_)))
            for(chr in 1:NC){
              xx=1+(chr-1)*Klink
              yy=chr*Klink
              oldrho[[chr]]=rho[xx:yy]
            }
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
          }
        }
      }
    }
    if(!ER){
      if(SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Klink=get('Klink', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              
              Beta=get('Beta', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[(Klink+1):(2*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              beta_=vector()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
                beta_[x:xx]=beta[ix]
                
              }
              LH=0
              for(chr in 1:NC){
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*Klink),1])
            beta_=sol[1:Klink]
            sigma_=sol[(1+Klink*NC):(2*Klink)]
            diff=max(abs(c(oldbeta- beta,oldsigma-sigma)))
            oldbeta=beta_
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC',envir = test.env)
              Klink=get('Klink',envir = test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              sigma=param[(Klink+1):(2*Klink)]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[((2*Klink+1)):(3*Klink)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              
              beta_=vector()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
                beta_[x:xx]=beta[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Self=get('Self', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(3*(Klink)),1])
            beta_=sol[1:Klink]
            sigma_=sol[(1+Klink):(2*Klink)]
            Xi_=sol[((2*Klink)+1):(3*Klink)]
            diff=max(abs(c(oldbeta-beta_,oldsigma-sigma_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            oldsigma=sigma_
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir = test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Klink=get('Klink', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              M=get('M', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                beta_[x:xx]=beta[ix]
              }
              for(chr in 1:NC){
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                m=as.vector(Big_M[[chr]])
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                if(BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta),function_to_minimize,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:Klink,1])
            beta_=sol[1:Klink]
            diff=max(abs(c(oldbeta-beta_)))
            oldbeta=beta_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC',envir = test.env)
              rho_=Rho
              BoxB=get('BoxB', envir=test.env)
              Klink=get('Klink', envir=test.env)
              beta=((param[1:Klink]*(BoxB[2]-BoxB[1]))+BoxB[1])^2
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[(1+Klink):(2*Klink)]
              Xi=vector()
              beta_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                beta_[x:xx]=beta[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              
              for(chr in 1:NC){
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta_,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta_+((1-beta_)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta_+((1-beta_)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
              }
              return(LH)
            }
            sol= BBoptim(c(oldbeta,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*(Klink)),1])
            beta_=sol[1:Klink]
            Xi_=sol[(1+Klink):(2*(Klink))]
            diff=max(abs(c(oldbeta-beta_,Xi_-oldXi_)))
            oldbeta=beta_
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
          }
        }
      }
      if(!SB){
        if(SF){
          if(Popfix){
            function_to_minimize <-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir=test.env)
              Rho=get('Rho', envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              Klink=get('Klink', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1:Klink]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              sigma_=vector()
              Pop=get('Pop', envir=test.env)
              LH=0
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                sigma_[x:xx]=sigma[ix]
              }
              for(chr in 1:NC){
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta,Pop = Pop,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                m=as.vector(Big_M[[chr]])
                m=m[keep]
                q_=get('q_', envir=test.env)
                nu=builder[[2]]
                if(BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldsigma),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:Klink,1])
            sigma_=sol[1:Klink]
            diff=max(abs(c(oldsigma-sigma_)))
            oldsigma=sigma_
          }
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              Rho=get('Rho', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              NC=get('NC',envir = test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Klink=get('Klink', envir=test.env)
              Boxs=get('Boxs', envir=test.env)
              Self=get('Self', envir=test.env)
              sigma=param[1:Klink]
              sigma=sigma*(Boxs[2]-Boxs[1])
              sigma=sigma+Boxs[1]
              
              Xi_=param[(Klink+1):(2*Klink)]
              Xi=vector()
              sigma_=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
                sigma_[x:xx]=sigma[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              
              for(chr in 1:NC){
                builder=build_HMM_matrix_t(n,(rho_[chr]),beta=beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma_,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
              }
              return(LH)
            }
            sol= BBoptim(c(oldsigma,oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:(2*(Klink)),1])
            sigma_=sol[1:Klink]
            Xi_=sol[(1+Klink):(2*Klink)]
            diff=max(abs(c(Xi_-oldXi_,oldsigma-sigma_)))
            oldXi_=Xi_
            oldsigma=sigma_
            print(sigma_)
            print(paste("Xi:",oldXi_))
          }
        }
        if(!SF){
          if(!Popfix){
            function_to_minimize_optim<-function(param){
              Boxr=get('Boxr', envir=test.env)
              mu=get('mu', envir=test.env)
              npair=get('npair', envir=test.env)
              Big_Window=get('Big_Window', envir=test.env)
              mu_b=get('mu_b', envir=test.env)
              FS=get('FS', envir=test.env)
              Rho=get('Rho', envir=test.env)
              NC=get('NC',envir=test.env)
              rho_=Rho
              beta=get('beta', envir=test.env)
              BoxP=get('BoxP', envir=test.env)
              Xi_=param[1:length(param)]
              Xi=vector()
              pop_vect=get('pop_vect', envir=test.env)
              xx=0
              for(ix in 1:Klink){
                x=xx+1
                xx = xx + pop_vect[ix]
                Xi[x:xx]=Xi_[ix]
              }
              Xi=Xi*sum(BoxP)
              Xi=Xi-(BoxP[1])
              Xi=10^Xi
              L=get('L', envir=test.env)
              n=get('k', envir=test.env)
              Big_Xi=get('Big_Xi', envir=test.env)
              Beta=get('Beta', envir=test.env)
              window_scaling=get('window_scaling', envir=test.env)
              BW=get('BW', envir=test.env)
              
              Pop=get('Pop', envir=test.env)
              LH=0
              
              Self=get('Self', envir=test.env)
              sigma=get('sigma', envir=test.env)
              cut_edge=get('cut_edge', envir=test.env)
              for(chr in 1:NC){
                builder=build_HMM_matrix(n,(rho_[chr]),beta,Pop = Pop,Xi=Xi,L=L[chr],Beta=Beta,scale=window_scaling,sigma = sigma,Sigma = Self,FS=FS,Big_Window=Big_Window,npair=npair,cut_edge =cut_edge)
                Q=builder[[1]]
                Q=t(Q)
                A=as.vector(Q)
                keep=which(A>0&as.vector(Big_Xi[[chr]])>0)
                A=A[keep]
                Big_Xi[[chr]]=as.vector(Big_Xi[[chr]])
                Big_Xi[[chr]]=Big_Xi[[chr]][keep]
                Tc=builder[[3]]
                g=matrix(0,nrow=length(Tc),ncol=2)
                if(!FS){
                  g[,2]=1-exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                  g[,1]=exp(-mu*(beta+((1-beta)*mu_b))*2*Tc)
                }
                if(FS){
                  g[,2]= 0.75 - (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                  g[,1]= 0.25 + (0.75*exp(-4*mu*(beta+((1-beta)*mu_b))*2*Tc/3))
                }
                x=as.vector(g)
                keep=which(x>0)
                x=x[keep]
                
                if(BW){
                  Big_M=get('Big_M', envir=test.env)
                  m=as.vector(Big_M[[chr]])
                  m=m[keep]
                  q_=get('q_', envir=test.env)
                  nu=builder[[2]]
                  LH=LH-sum(log(A)*Big_Xi[[chr]])-sum(log(x)*m)-sum(log(nu)*q_[[chr]])
                }
                if(!BW){
                  LH=LH-sum(log(A)*Big_Xi[[chr]])
                }
                
              }
              return(LH)
            }
            sol= BBoptim(c(oldXi_),function_to_minimize_optim,lower=0,upper=1,method=c(2),control = list(maxit=30,M=c(20)))
            LH=as.numeric(as.matrix(sol[[2]]))
            sol=as.matrix(sol[[1]])
            sol=as.numeric(sol[1:((Klink)),1])
            Xi_=sol[1:length(sol)]
            diff=max(abs(Xi_-oldXi_))
            oldXi_=Xi_
            print(paste("Xi:",oldXi_))
            
          }
        }
      }
    }
  }
  if(SB){
    beta_=((oldbeta*(BoxB[2]-BoxB[1]))+BoxB[1] )^2
    beta=vector()
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      beta[x:xx]=beta_[ix]
    }
  }
  if(SF){
    sigma_=oldsigma*(Boxs[2]-Boxs[1])
    sigma_=sigma_+Boxs[1]
    sigma=vector()
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      sigma[x:xx]= sigma_[ix]
    }
  }
  if(NC==1){
    if(ER){
      rho_=oldrho*sum(Boxr)
      rho_=rho_-(Boxr[1])
      rho_=10^(rho_)
      rho_=rho_*Rho
      rho=vector()
      xx=0
      for(ix in 1:Klink){
        x=xx+1
        xx = xx + pop_vect[ix]
        rho[x:xx]=rho_[ix]
      }
    }else{
      rho=Rho
      rho_=Rho
    }
  }
  if(NC>1){
    if(ER){
      rho_=list()
      rho=list()
      for(chr in 1:NC){
        rho_[[chr]]=oldrho[[chr]]*sum(Boxr)
        rho_[[chr]]=rho_[[chr]]-(Boxr[1])
        rho_[[chr]]=10^(rho_[[chr]])
        rho_[[chr]]=rho_[[chr]]*Rho[chr]
        rho[[chr]]=vector()
        xx=0
        for(ix in 1:Klink){
          x=xx+1
          xx = xx + pop_vect[ix]
          rho[[chr]][x:xx]=rho_[[chr]][ix]
        }
      }
    }else{
      rho=Rho
      rho_=Rho
    }
  }
  if(Pop){
    xx=0
    for(ix in 1:Klink){
      x=xx+1
      xx = xx + pop_vect[ix]
      oldXi[x:xx]=oldXi_[ix]
    }
    Xi_=oldXi*sum(BoxP)
    Xi_=Xi_-(BoxP[1])
    Xi_=10^Xi_
  }
  res<-list()
  res$beta=beta
  res$sigma=sigma
  
  
  if(Pop){
    res$Xi=Xi_
  }
  Beta=get('Beta', envir=test.env)
  Self=get('Self', envir=test.env)
  if(NC==1){
    rho_=rho_/(2*L)
    builder=build_HMM_matrix_t(n,(rho_*2*L),beta,Pop = Pop,Xi=Xi_,L=L,Beta=Beta,scale=window_scaling,sigma=sigma,Sigma=Self,FS=FS,Big_Window=Big_Window,npair=npair) # (n=20,rho,beta=1,L,Pop=T,Xi=NA,Beta=1,scale=c(1,0),sigma=0,Sigma=0,FS=F,Big_Window=F,tmax=15,alpha_t=0.1,npair=2,cut_edge=F)
    Tc=builder[[4]]
  }
  if(NC>1){
    Tc=list()
    rho_=rho_/(2*L[chr])
    for(chr in 1:NC){
      builder=build_HMM_matrix_t(n,(rho_[chr]*2*L[chr]),beta,Pop = Pop,Xi=Xi_,L=L,Beta=Beta,scale=window_scaling,sigma=sigma,Sigma=Self,FS=FS,Big_Window=Big_Window,npair=npair)
      Tc[[chr]]=builder[[4]]
    }
    
    
    
    
  }
  res$rho=rho_
  res$Tc=Tc
  res$mu=mu
  res$N=N
  return(res)
}

