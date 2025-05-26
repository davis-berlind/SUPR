results = Reduce("rbind", lapply(list.files()[1:10], function(file) read.csv(file)))
saveRDS(results, "~/MICH/results.rds")