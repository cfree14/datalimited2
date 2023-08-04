
#' Retrieve resilience estimates from FishBase
#'
#' Retrieves estimates of resilience for fish and invertebrate species from FishBase
#' and SeaLifeBase. Resilience describes the ability of a species population to recover
#' after a perturbance and is reported as: very low, low, medium, or high.
#'
#' @param species A character vector of species
#' @return A dataframe containing species and their resilience values
#' @examples
#' # Retrieve resilience for three species
#' resilience(species=c("Gadus morhua", "Mola mola", "Centropristis striata"))
#' @export
resilience <- function(species){

  # Build FB/SLB taxa key
  taxa_key_fb <- rfishbase::fishbase %>% # (1-old) rfishbase::load_taxa(server="fishbase") (2-old)# server="https://fishbase.ropensci.org"
    as.data.frame() %>%
    mutate(type="fish") %>%
    select(type, everything()) %>%
    setNames(tolower(colnames(.))) %>%
    mutate(sciname=paste(genus, species)) %>%
    mutate_all(as.character)
  taxa_key_slb <- rfishbase::sealifebase %>% # (1-old) rfishbase::load_taxa(server="sealifebase") (2-old) rfishbase::sealifebase %>%
    as.data.frame() %>%
    mutate(type="invert") %>%
    select(type, everything()) %>%
    setNames(tolower(colnames(.))) %>%
    mutate(sciname=paste(genus, species)) %>%
    mutate_all(as.character)
  taxa_key <-  taxa_key_fb %>%
    bind_rows(taxa_key_slb) %>%
    setNames(tolower(names(.))) %>%
    select(type, class, order, family, genus, species, sciname) %>%
    unique()

  # Check whether species are in FishBase
  spp <- species
  spp_do <- spp[spp %in% taxa_key$sciname]
  spp_not_in_fb <- spp[!spp %in% taxa_key$sciname]
  if(length(spp_not_in_fb)>0){
    print(paste0("The following species are not in FishBase: ", paste(spp_not_in_fb, collapse=", ")))
  }

  # Divide into fish and inverts
  taxa_do <- taxa_key %>%
    filter(sciname %in% spp_do)
  spp_fin <- taxa_do$sciname[taxa_do$type=="fish"]
  spp_inv <- taxa_do$sciname[taxa_do$type=="invert"]

  # Finfish
  #######################################

  # If there are finfish
  if(length(spp_fin)>0){

    # Get resilience info from Fishbase
    fin_orig <- rfishbase::stocks(spp_fin, server="fishbase")
    lh_fin <- fin_orig %>%
      select(Species, Resilience) %>%
      filter(!is.na(Resilience)) %>%
      rename(species=Species, resilience=Resilience) %>%
      mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")))

    # Calculate mode (default is to resolve ties with lower value)
    lh_fin_sum <- as.data.frame.matrix(table(lh_fin$species, lh_fin$resilience))
    lh_fin_sum$resilience <- apply(lh_fin_sum, 1, function(x) colnames(lh_fin_sum)[which.max(x)])

    # Format for export
    lh_fin1 <- lh_fin_sum %>%
      mutate(species=row.names(lh_fin_sum)) %>%
      select(species, resilience)

    # Invertebrates
    #######################################

    # If there are inverts
    if(length(spp_inv)>0){

      # Get resilience info from Fishbase
      inv_orig <- rfishbase::stocks(spp_inv, server="sealifebase")
      lh_inv <- inv_orig %>%
        select(Species, Resilience) %>%
        filter(!is.na(Resilience)) %>%
        rename(species=Species, resilience=Resilience) %>%
        mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")))

      # Calculate mode (default is to resolve ties with lower value)
      lh_inv_sum <- as.data.frame.matrix(table(lh_inv$species, lh_inv$resilience))
      lh_inv_sum$resilience <- apply(lh_inv_sum, 1, function(x) colnames(lh_inv_sum)[which.max(x)])

      # Format for export
      lh_inv1 <- lh_inv_sum %>%
        mutate(species=row.names(lh_inv_sum)) %>%
        select(species, resilience)

    }

  }

  # Merge results
  #######################################

  # Merge finfish and inverts
  if(length(spp_fin)>0 & length(spp_inv)>0){ # fish and inverts
    out <- rbind(lh_fin1, lh_inv1) %>% arrange(species)
  }
  if(length(spp_fin)>0 & length(spp_inv)==0){ # only fish
    out <- lh_fin1
  }
  if(length(spp_fin)==0 & length(spp_inv)>0){ # only inverts
    out <- lh_inv1
  }

  # Make resilience a factor
  out$resilience <- factor(out$resilience, levels=c("Very low", "Low", "Medium", "High"))

  # Return
  return(out)

}
