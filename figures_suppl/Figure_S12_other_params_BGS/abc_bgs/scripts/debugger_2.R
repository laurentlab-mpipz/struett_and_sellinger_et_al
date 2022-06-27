
c = readRDS("/home/struett/dadadada.RDS")

a <- c[[1]]
b <- c[[2]]

colnames(a) <- colnames(b)

df.estim.mode = a
params_pod = b

get_rb(df.estim.mode, params_pod)
