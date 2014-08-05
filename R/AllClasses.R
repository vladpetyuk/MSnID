

.PROTON_MASS <- 1.007276466812

setClass(Class="MSnID",
            representation(
                workDir="character", # working directory
                psms="data.table") # peptide-to-spectra matches
)




# partitions and parameter transforms to be added to filter class
setClass(Class="MSnIDFilter",
            representation(
                filterList="list",
                validParNames="character")
)

