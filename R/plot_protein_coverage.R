


utils::globalVariables(c("Last_AA_First",
                         "First_AA_First",
                         "ProtLen",
                         "Length",
                         "n",
                         "ymin",
                         "ymax"))



.plot_protein_coverage <- function(object, accession, save_plot = FALSE, ...){
    
    x <- psms(object) %>%
        select(-c(First_AA, Last_AA)) %>%
        dplyr::rename(Last_AA = Last_AA_First,
                      First_AA = First_AA_First)
    
    prot_len <- filter(x, accession == !!accession) %>% 
        distinct(ProtLen) %>% 
        as.numeric()
    
    prot_name <- accession
    
    prot <- x %>%
        filter(accession == !!accession) %>%
        # dplyr::rename(`First_residue` = `First residue`) %>%
        # dplyr::rename(`Last_residue` = `Last residue`) %>%
        group_by(First_AA, Last_AA) %>%
        tally() %>%
        ungroup() %>%
        mutate(Length = Last_AA - First_AA + 1) %>%
        arrange(First_AA, -Length, -n)
    
    
    # setting staggered ymin
    min_y <- 0.1
    step_y <- 0.033
    width_y <- 0.025
    
    prot$ymin <- min_y
    if(nrow(prot) > 1){
        for(i in 2:nrow(prot)){
            current_y <- min_y
            while(TRUE){
                # is there a conflict
                max_last_residue <- prot %>%
                    dplyr::slice(1:(i-1)) %>%
                    dplyr::filter(ymin == current_y) %>%
                    pull(Last_AA) %>%
                    max()
                if(max_last_residue + 0 >= prot[i,"First_AA",drop=TRUE]){
                    current_y <- current_y + step_y
                }else{
                    break()
                }
            }
            prot[i,"ymin"] <- current_y
        }
    }
    
    
    prot$ymax <- prot$ymin + width_y
    p <-
        ggplot(data = prot) +
        geom_rect(aes(xmin = 0, xmax = prot_len + 2, ymin = -0.04, ymax = +0.04)) +
        geom_rect(aes(xmin=First_AA, xmax=Last_AA, ymin=ymin, ymax=ymax, fill=n),
                  color="white", size=1) +
        scale_fill_viridis_c(option = "D")
    
    p <- p +    
        ylab(NULL) +
        xlab("residue") +
        theme_classic() +
        theme(axis.ticks.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.line.y = element_blank(),
              axis.line.x = element_blank()) +
        scale_x_continuous(breaks = seq(0,prot_len,20)) +
        theme(axis.text.x = element_text(angle=90),
              plot.title = element_text(hjust = 0.5, size=16)) +
        ggtitle(prot_name)
    
    if(max(prot$ymax) < 0.5)
        p <- p + ylim(-0.04, 0.5)
    
    if(save_plot){
        file_name <- gsub(", ", "_", prot_name)
        file_name <- gsub("\\|", "_", file_name) # OneDrive does not allow pipes
        ggsave(filename = paste0(file_name,".png"), plot = p)
    }else{
        return(p)
    }
}





