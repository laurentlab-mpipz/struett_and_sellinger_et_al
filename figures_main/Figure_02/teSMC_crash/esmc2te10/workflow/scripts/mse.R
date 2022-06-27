
# save image for debugging
#save.image(); stop("\n_____________\n saved image\n=============\n\n")
# setwd("../../.")
# load(".RData")

result_list = readRDS(snakemake@input$full_data)

total = result_list[[2]]
mu_ex = result_list[[4]]
s_ex = result_list[[3]]
Tc = result_list[[5]]

# this is the script of Thibaut to calculate the mse
nsim = 1
mu = as.numeric(as.character(snakemake@wildcards[["mutation_rate"]]))

mat_save_p=matrix(NA,nrow=(1*1*nsim),ncol=30)
mat_save_t=matrix(NA,nrow=(1*1*nsim),ncol=30)
mat_save_s=matrix(NA,nrow=(1*1*nsim),ncol=30)
count_p=0
count_t=0
count_s=0

x = 1 # loop if nsim > 1

test_<-total[[1+((x-1)*4)]]
Pop_=mu_ex[[1+((x-1)*4)]]/(mu)
count_p=count_p+1
count_t=count_t+1
count_s=count_s+1

mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
mat_save_t[count_t,]=(Tc[[1+((x-1)*4)]]*Pop_)
mat_save_s[count_s,]=s_ex[[1+((x-1)*4)]]

# as I run the different inference models separately, I do not need this part
# test_<-total[[2+((x-1)*4)]]
# Pop_=mu_ex[[2+((x-1)*4)]]/(mu)
# count_p=count_p+1
# count_t=count_t+1
# count_s=count_s+1
# mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
# mat_save_t[count_t,]=(Tc[[2+((x-1)*4)]]*Pop_)
# mat_save_s[count_s,]=s_ex[[2+((x-1)*4)]]
# 
# test_<-total[[3+((x-1)*4)]]
# Pop_=mu_ex[[3+((x-1)*4)]]/(mu)
# count_p=count_p+1
# count_t=count_t+1
# count_s=count_s+1
# mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
# mat_save_t[count_t,]=(Tc[[3+((x-1)*4)]]*Pop_)
# mat_save_s[count_s,]=s_ex[[3+((x-1)*4)]]
# 
# test_<-total[[(4*x)]]
# Pop_=mu_ex[[(4*x)]]/(mu)
# count_p=count_p+1
# count_t=count_t+1
# count_s=count_s+1
# mat_save_p[count_p,]=(log10((test_)*Pop_)*0.5)
# mat_save_t[count_t,]=(Tc[[(4*x)]]*Pop_)
# mat_save_s[count_s,]=s_ex[[(4*x)]]

# before we created the matrix, now we calculate the 

message("the true simulated values are very hard coded")

count_p=0
count_t=0
count_s=0
MSE_pop=matrix(0,nrow = nsim,ncol = 1)
MSE_selfing=matrix(0,nrow = nsim,ncol = 1)
for(x_rep in 1:nsim){
    count_p=count_p+1
    count_t=count_t+1
    count_s=count_s+1
            
    nb_point=10^4 # Number of point where MSE shoulb be calculated
    x_time=seq(log10((as.numeric(mat_save_t[count_t, 2]))),
        log10(max(as.numeric(mat_save_t[count_t,]))),
        (log10(max(as.numeric(mat_save_t[count_t,]))) - 
            log10((as.numeric(mat_save_t[count_t,2]))) )/nb_point)
    x_time=10^x_time
    #Pop_size= # simulated population size at each x_time time
    current_selfing_value=0.99
    ancestral_selfing_value=0.1


    # get right demography
    which_demography = as.numeric(snakemake@wildcards$demography)+1
    demography = matrix(snakemake@config$population_sizes_backward_in_time,
                          byrow = T, ncol = 2)[which_demography, 1]


    Pop_size=demography

    # time_of_change:
    which_demography = as.numeric(snakemake@wildcards$demography)+1
    time_of_change=matrix(snakemake@config$selfing_rates_backward_in_time,
                          byrow = T, ncol = 4)[which_demography,4]
    sigma_t=rep(current_selfing_value,length(x_time))
    if (any(x_time>time_of_change)) {
      sigma_t[which(x_time>time_of_change)]=ancestral_selfing_value
    }


    Estimated_pop=numeric(length(x_time))
    Estimated_selfing=numeric(length(x_time))
    count_x=0
        for(xx in x_time){
            count_x=count_x+1
            Estimated_pop[count_x]=     mat_save_p[count_p,max(which((as.numeric(mat_save_t[count_t,]))<=xx))]
            Estimated_selfing[count_x]= mat_save_s[count_s,max(which((as.numeric(mat_save_t[count_t,]))<=xx))]
      }
    MSE_pop[x_rep,1]=log10(mean((Pop_size- (10^Estimated_pop))^2))
    MSE_selfing[x_rep,1]=log10(mean((sigma_t- (10^Estimated_selfing))^2))
}

# setwd(" ") # Location to put MSE
filename_MSE_pop=snakemake@output$demgr
write.csv(MSE_pop,filename_MSE_pop)
      
filename_MSE_selfing= snakemake@output$sigma
write.csv(MSE_selfing,filename_MSE_selfing)

