lik_reloaded.old <- readRDS("lik.old.RData")       #object from new upd function
lik2_reloaded.old <- readRDS("lik2.old.RData")    #object from original upd function

sim_reloaded.old <- readRDS("sim.RData")           #object from new upd function
sim2_reloaded.old <- readRDS("sim2.RData")         #object from original upd function


setequal(lik_reloaded.old,lik2_reloaded.old) #checks to see if the two objects are equal
setequal(sim_reloaded.old,sim2_reloaded.old) #checks to see if the two objects are equal






