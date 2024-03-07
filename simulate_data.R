rm(list=ls())

library(Rcpp)
sourceCpp("Transfil.cpp")

bednet_settings <- list("time" = c(200,300),
                        "coverage" = c(0.3,1))

MDA_settings <- list("time" = c(100,112,124),
                     "coverage" = c(0.65,0.65,0.65))

pop_size <- 1000
x <- new(Data, pop_size)
x$set_bednet(bednet_settings)
x$set_MDA(MDA_settings)
x$set_parameter("k", 0.3)
x$set_parameter("vth", 120)
x$set_parameter("importation_rate", 0.005)
x$burn_in()
x$timestep(30*12) # in months

plot(x$mf_prev, type="l", ylab="Mf prevalence", xlab = "Time (months)")
