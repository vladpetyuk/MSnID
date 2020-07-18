
rec_insert <- function(x, y, s){
    if(length(y) == 0)
        return(x)
    else{
        x <- paste0(substr(x, 0, y[1]), s, substr(x, y[1] + 1, nchar(x)))
        y <- y[-1]
        rec_insert(x, y, s)
    }
}



.add_mod_symbol <- function(obj, mod_mass, symbol){
    
    #
    # get the positions of PTMs of the specified mass
    mod_pos <- obj$modification
    # Bunch of temp vars because I am not ready to use %>% pipe in this package
    t1 <- strsplit(mod_pos, mod_mass)
    t2 <- lapply(t1, Filter, f=function(x) grepl("^\\s",x))
    t3 <- lapply(t2, sub, pattern = "^\\s\\((\\d+).*",replacement = "\\1")
    t4 <- lapply(t3, as.numeric)
    t5 <- lapply(t4, `+`, 2) # account for X. at N-term
    t6 <- lapply(t5, function(..) .. + seq_along(..) - 1) # account for shifts
    mod_pos <- t6
    
    # map to the peptide sequences
    pep <- obj$peptide_mod
    if(is.null(pep)){
        pep <- obj$peptide
    }else{
        # how account for already placed mods?
        # get positions of non-AA symbols
        prior_mod_pos <- lapply(pep, gregexpr, pattern = "[^ARNDCQEGHILKMFPSTWYV.-]")
        prior_mod_pos <- lapply(prior_mod_pos, `[[`, 1)
        prior_mod_pos <- lapply(prior_mod_pos, as.numeric)
        prior_mod_pos <- lapply(prior_mod_pos, function(..) ifelse(.. < 0, NA, ..))
        prior_mod_pos <- lapply(prior_mod_pos, function(..) .. - seq_along(..)) # back to AA ind
        
        for(i in seq_along(mod_pos)){
            y <- mod_pos[[i]]
            d <- prior_mod_pos[[i]]
            for(j in seq_along(y)){
                shift <- sum(y[j] >= d, na.rm = TRUE)
                y[j] <- y[j] + shift
            }
            mod_pos[[i]] <- y
        }
    }
    
    # insert the current ones
    obj$peptide_mod <- sapply(seq_along(pep), function(..) {
        rec_insert(pep[..], mod_pos[[..]], symbol)})
    return(obj)
    
}




