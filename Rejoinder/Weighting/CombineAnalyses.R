load(".//New-BartoModel//new.RData")
saveRDS(ests.mat.new, file = "barto-0-unweighted.rds")

load(".//New-BartoModelw//new.RData")
saveRDS(ests.mat.new, file = "barto-1-cfweighted.rds")

load(".//New-BartoModelw2S//new.RData")
saveRDS(ests.mat.new, file = "barto-2-2sweighted.rds")

load(".//BB//BB-stat13//.RData")
saveRDS(est.mat, file = "barto-3-bb-linked.rds")

load(".//BB-Unlinked//.RData")
saveRDS(est.mat, file = "barto-4-bb-unlinked.rds")

load(".//New-CorrectModel//correct.RData")
saveRDS(ests.mat.new, file = "barto-c-correct.rds")

save.image("CombinedAnalysis.RData")
load("CombinedAnalysis.RData")
rm(list=ls())
ests0 <- readRDS("barto-0-unweighted.rds")
ests1 <- readRDS("barto-1-cfweighted.rds")
ests2 <- readRDS("barto-2-2sweighted.rds")
ests3 <- readRDS("barto-3-bb-linked.rds")
ests4 <- readRDS("barto-4-bb-unlinked.rds")
estsc <- readRDS("barto-c-correct.rds")
save.image("CombinedAnalysis.RData")





