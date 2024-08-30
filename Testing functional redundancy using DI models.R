library(DImodels)
library(dplyr)

# Condition 1 delta_ik = delta_jk, i.e., between_FG interactions
# Condition 2 delta_ij = 0, i.e., don't interaction with each other
# Condition 3 beta_i = beta_j

#' @param model A DI model object used to test functional redundancy in. Note currently only DI models with theta = 1 are supported
#' @param groups A character vector of the same length as the number of species (and in same order) 
#'               in the experiment specifying the functional groupings of the species. For example,
#'               to test for functional redundancy between grasses in a system with four species, 
#'               groups could be c("Gr", "Gr", "p3", "p4"). The species to be tested for functional
#'               redundancy should be given the same characters.
#' @param selection Selection method to be used in the automated model selection process. 
#'                  Options are "Ftest", "AIC", "AICc", "BIC" and "BICc". 
#'                  The default is selection = "Ftest".
#' @param verbose A logical value indicating whether to print detailed anova tables. Default is TRUE
#' 
#' Additional parameters passed to ..., these accept all parameters accepted by the DI() function.
#' This is useful if you wish to test for functional redundancy but don't have a fitted model.
#' 
#' If not specifying a model object, these parameter are mandatory
#' @param y The column name of the response vector, which must be in quotes, for example, y = "yield".
#' @param prop A vector of s column names identifying the species proportions in each community in 
#'             the dataset. For example, if the species proportions columns are labelled p1 to p4, 
#'             then prop = c("p1","p2","p3","p4"). Alternatively, the column numbers can be specified,
#'             for example, prop = 4:7, where species proportions are in the 4th to 7th columns.
#' @param data Specify the dataset, for example, data = Switzerland. The dataset name should not 
#'             appear in quotes.
#'
#' Additional optional parameters
#' @param DImodel This argument is chosen (over custom_formula) to fit an automated version of a DI model. 
#'                The chosen tag determines the form of the species interactions to be included in the model. 
#'                The tags (or options) are:
#'                  + STR (no identity or interaction effects, only an intercept is fitted, 
#'                    plus the experiment structural variables block, density and treat, if specified).
#'                Each of the following includes the species proportions as specified in prop, the 
#'                interaction variables according to the tag, plus block, density and treat if specified.
#'                  + ID (no interaction terms),
#'                  + AV (a single average species pairwise interaction variable), 
#'                  + FG (functional group interaction variables, the FG argument 
#'                    must be specified to use this option),
#'                  + ADD (the additive species interaction variables),
#'                  + FULL (all pairwise interaction variables).
#'                The DImodel tag should appear in quotes, for example, DImodel = "STR".
#' @param extra_formula n conjunction with DImodel, additional terms can be added using extra_formula. 
#'                      A ~ must be included before specifying the terms. For example, if DImodel = "AV" 
#'                      has been specified, adding extra_formula = ~ I(AV**2) will add a quadratic average 
#'                      pairwise interaction variable or extra_formula = ~ treatment:AV will add an 
#'                      interaction between the average pairwise species interaction variable and the treatment. 
#'                      Any variable included directly in extra_formula must already be contained in the dataset 
#'                      (interaction variables can be created using the function DI_data, if required).
#' @param FG If species are classified by g functional groups, this argument takes a text list (of length s) 
#'           of the functional group to which each species belongs. For example, for four grassland species 
#'           with two grasses and two legumes: FG could be FG = c("G","G","L","L"), where G stands for grass 
#'           and L stands for legume.
#'           The FG argument is required if DImodel = "FG" is specified.
#' @param ID 	This argument takes a text list (of length s) dsecirbing groupings for the identity effects 
#'            of the species. For example, if there are four species and you wish to group the identity 
#'            effects all four species into a single term: ID could be ID = c("ID1","ID1","ID1","ID1"), 
#'            where "ID1" is the name of the ID group. Similarly if the we wish to have two identity effect 
#'            groups where identity effect of species 1 and 3, and species 2 and 4 are grouped together: 
#'            ID could be ID = c("ID1","ID2","ID1","ID2"), where "ID1" and "ID2" are the names of the ID 
#'            groups. These ideas expand to any number of species and any number or combination of groups. 
#'            Finally, the ID groups do not have to be named "ID1" and "ID2", the user can specify any name 
#'            for the groups.
#'              + If the ID argument is not specified, each species will be assumed to have a 
#'                separate identity effect.
#'              + Specify an grouping for the ID does not affect the interaction terms. The interactions are 
#'                still calculated using the individual species proportions.
#'              + The ID argument is defunct when using the custom_formula argument, since species identity 
#'                effects must be included directly in the custom_formula argument.
#' @param block The name of the block variable (if present), which must be in quotes, for example, 
#'              block = "block". If no blocking variable, omit this argument.
#' @param density The name of the density variable (if present), which must be in quotes, 
#'                for example, density = "density". If no density variable, omit this argument.
#' @param treat The name of a column in the dataset containing the value of a treatment factor or covariate. 
#'              The treatment name must be included in quotes, for example, treat = "nitrogen". If the
#'              treatment is a factor, the variable must already be specified as a factor prior to using DI.
#'                + When used in conjunction with DImodel, the treatment will be included in the model as an 
#'                  additive factor or covariate, for example, specifying treat = nitrogen, DImodel = ID will 
#'                  fit the model p1 + p2 + ... + ps + nitrogen. Additional treatments, or interactions between 
#'                  the treatment and other model terms can be included via the extra_formula argument.
#'                + The treat argument is defunct when using the custom_formula argument, and any treatment must 
#'                  be included directly in the custom_formula argument.
#' 
test_functional_redundancy <- function(model, groups,
                                       selection = c("F-test", "AIC", "BIC", "AICc", "BICc"),
                                       verbose = TRUE,
                                       ...){
  selection <- match.arg(selection)
  DI_args <- list(...)
  DI_args$estimate_theta <- FALSE
  DI_args$theta <- 1
  if(is.null(DI_args$block)) DI_args$block <- NA
  if(is.null(DI_args$density)) DI_args$density <- NA
  if(is.null(DI_args$treat)) DI_args$treat <- NA
  # browser()
  if(missing(model)){
    if(is.null(DI_args$y) || is.null(DI_args$prop) || is.null(DI_args$data)){
      stop(paste0("If not specifying a model object, it is mandatory to specify", 
                  "the default arguments of the DI function. ",
                  "Please specify the appropriate values for `y`, `prop`, and `data`."))
    }
    if(is.null(DI_args$DImodel)){
      # FULL model is not feasible if there are more than 10 species
      if(length(DI_args$prop) > 10) {
        warning("The number of pairwise interactions would be too high for ecosystems with over 10 species ",
                "the full pairwise model will not be fit and the functional group model would be ",
                "assumed to be the complex model in this situation. ",
                "Manually specify a value in the `DImodel` parameter incase this is not desirable.")
        DI_args$no_full <- TRUE
        DI_args$DImodel <- "AV"
      } else {
        DI_args$DImodel <- "FULL" 
      }
    }
    
    model <- suppressWarnings(
      suppressMessages(       
             DI(data = DI_args$data, y = DI_args$y, prop = DI_args$prop,
                block = DI_args$block, density = DI_args$density, 
                treat = DI_args$treat,
                estimate_theta = DI_args$estimate_theta, theta = DI_args$theta,
                extra_formula = DI_args$extra_formula,
                DImodel = DI_args$DImodel)
        )
      )
  } 
  
  # Ensure model is a DI models  
  if(!inherits(model, "DI")){
    stop("Specify a model of class `DI` fit using the `DI()` or `autoDI()` functions")
  }
  
  # Ensure theta equals 1
  if(!is.na(attr(model, "theta_val")) && attr(model, "theta_val") != 1){
    stop("Currently functional redundancy can only tested for DI models where theta = 1")
  }
  
  if(attr(model, "DImodel") != "FULL" && !isTRUE(DI_args$no_full) && selection == "F-test"){
    warning("The specified model for testing functional redudandcy against is not nested with ",
            "the models fit internally in the function. Therefore, F-tests can not be used for ",
            "model selection, using 'AICc'. Choose a different information criteria using `selection` parameter.")
    selection <- "AICc"
  }
  # Attributes of the full model
  data <- model$original_data
  mod_args <- attributes(model)
  extra_formula <- eval(model$DIcall$extra_formula) 
  mod_args$block <- if(!is.null(model$DIcall$block)) eval(model$DIcall$block) else NA
  mod_args$density <- if(!is.null(model$DIcall$density)) eval(model$DIcall$density) else NA
  mod_args$treat <- if(!is.null(model$DIcall$treat)) eval(model$DIcall$treat) else NA
  
  if(length(groups) != length(mod_args$prop)){
    stop("The number of elements in the `groups` argument should be the same as the number of species in the model.")
  }
  
  # Model fulfilling first condition, i.e. delta_ik = delta_jk 
  # This is the FG model
  # Use theta estimate from comp_model
  theta <-  1 #if(!is.na(mod_args$theta_val)) as.numeric(mod_args$theta_val) else DI_args$theta
  cond1_mod <- suppressWarnings(
    suppressMessages(
          DI(data = data, y = mod_args$y, prop = mod_args$prop,
             block = mod_args$block, density = mod_args$density,
             treat = mod_args$treat,
             estimate_theta = FALSE, theta = theta,
             extra_formula = extra_formula,
             DImodel = "FG", FG = groups) 
      )  
    )
  
  # Model additionally fulfilling second condition, i.e. where delta_ij = 0
  ## Manually pull out the FG interactions
  FG_ints <- DI_data_FG(prop = mod_args$prop, FG = groups, data = data, 
                        theta = theta)$FG
  ## Set within FG interactions to 0 by removing them
  bfg_ints <- paste0("FG_", colnames(FG_ints)[!startsWith(colnames(FG_ints), "wfg")])
  ## Add FG terms to data
  mod2_data <- cbind(data, `colnames<-`(FG_ints, paste0("FG_", colnames(FG_ints))))
  
  ## Append extra formula
  extra_formula_more <- as.formula(paste0("~ ",
                                          paste0(c(extra_formula[[length(extra_formula)]],
                                                   bfg_ints),
                                                 collapse = " + ")))
  # Same as cond1_mod, but now within FG terms are set to 0
  cond2_mod <- suppressWarnings(
      suppressMessages(
        DI(data = mod2_data, y = mod_args$y, prop = mod_args$prop, 
           DImodel = "ID", block = mod_args$block, 
           density = mod_args$density, treat = mod_args$treat,
           estimate_theta = FALSE, theta = theta,
           extra_formula = extra_formula_more)
      )
    )
  
  # Model additionally fulfilling third condition, i.e. beta_i = beta_j
  # Same as cond2_mod, but now the ID effects are grouped
  cond3_mod <- suppressWarnings(
      suppressMessages(
        DI(data = mod2_data, y = mod_args$y, prop = mod_args$prop, 
           DImodel = "ID", block = mod_args$block, 
           density = mod_args$density, treat = mod_args$treat,
           estimate_theta = FALSE, theta = theta,
           ID = groups, extra_formula = extra_formula_more)
      )
    )
  
  mods_list <- list(cond3_mod, cond2_mod, cond1_mod, model)
  desc <- if(isTRUE(DI_args$no_full)) c(rep(TRUE, 3), FALSE) else rep(TRUE, 4)
  
  selected <- redundancy_table(mods_list, groups = groups, selection = selection,
                               descriptions = desc, verbose = verbose)
  selected_model <- selected$selected_model
  mod_code <- selected$code
  # browser()
  # Restimate theta for when functionally redundancy is present
  if(mod_code %in% c("Model a", "Model b", "Model c") && isTRUE(mod_args$theta_flag)){
    if(mod_code != "Model c") {
      IDs <- attr(selected_model, "ID")
      selected_model <- suppressWarnings(
        suppressMessages(
          custom_theta(obj = selected_model, ID = IDs, 
                                     FG = groups, FG_in_model = bfg_ints,
                                     block = mod_args$block, treat = mod_args$treat, 
                                     density = mod_args$density,
                                     extra_formula = extra_formula_more)
        )
      )
    } else {
     selected_model <- suppressWarnings(
       suppressMessages(
         DI(data = data, y = mod_args$y, prop = mod_args$prop,
            block = mod_args$block, density = mod_args$density,
            treat = mod_args$treat,
            estimate_theta = TRUE, 
            extra_formula = extra_formula,
            DImodel = "FG", FG = groups) 
       )  
     )
    }
  }
  return(selected_model)
}

# Helper functions run these and skip to line 462
AICc <- DImodels:::AICc
BICc <- DImodels:::BICc

redundancy_table <- function(model_list, groups,
                             selection = c("F-test", "AIC", "BIC", "AICc", "BICc"),
                             descriptions = c(TRUE, TRUE, TRUE, TRUE),
                             verbose = TRUE){
  # Partial matching for selection metric
  selection <- match.arg(selection)
  
  # Describing the conditions
  if(isTRUE(verbose)){
    message("For two species i and j to be functionally redundant the following conditions must be met.\n\n",
            "(1) They should interact with other species (k) in the same way (i.e., \U03B4ik = \U03B4jk)\n",
            "(2) They should not interact with each other (i.e., \U03B4ij = 0)\n",
            "(3) Species should have same the identity effects (i.e., \U03B2i = \U03B2j)",
            "\nThese concepts extend to test for functional redundancy between more than 2 species.")
    cat(paste0("\033[0;", 32, "m",
               paste0("\nTesting functional redundancy as per the following grouping\nc(",
                      paste('"', paste(groups, collapse='", "'), '"', sep=''), ")"),
               "\033[0m","\n"))
    message("\nModel selection\n")
  }
  
  
  # Hard coded descriptions for the models
  # desc <- c("Model with grouped ID effects (i.e., \U03B2i = \U03B2j), between FG interactions (i.e., \U03B4ik = \U03B4jk ) and within FG interactions set to 0 (i.e., \U03B4ij = 0)",
  #           "Model with grouped ID effects (i.e., \U03B2i = \U03B2j) and between FG interactions (i.e., \U03B4ik = \U03B4jk )",
  #           "Model with grouped ID effects (i.e., \U03B2i = \U03B2j) but all interactions",
  #           "Full model")
  mod_d_desc <- switch (attr(model_list[[4]], "DImodel"),
                        "FULL" = "Full model with separate identity and interaction effects",
                        "AV" = "Model with a single average interaction term (AV model)",
                        "ADD" = "Model with additive species interactions (ADD model)",
                        "ID" = "Model with only species identity effects and no interactions (ID model)",
                        "FG" = "Model with functional group interaction terms (FG model)",
                        "Custom user model"
  )
  
  desc <- c("Model satisfying conditions 1, 2, and 3",
            "Model satisfying conditions 1 and 2",
            "Model satisfying only condition 1 (FG model)",
            mod_d_desc)
  
  final_table <- data.frame(selection = sapply(model_list, function(x) selection)) #,
  #treat = "none", #sapply(model_list, function(x) attr(x, "treat")),
  #theta = sapply(model_list, function(x) attr(x, "theta_flag")))
  rnames <- paste0("Model ", letters[1:length(desc)])
  rownames(final_table) <- rnames
  desc <- desc[descriptions]
  final_table <- slice(final_table, c(1:nrow(final_table))[descriptions])
  model_list <- model_list[descriptions]
  if(selection == "F-test"){
    # Family lock change later if more families are added
    Test <- "F" # Could also be Chisq for non-gaussian families
    final_table$Description <-  desc
    anovas <- eval(parse(text = paste("anova(",
                                      paste("model_list[[", 1:length(model_list), "]]",
                                            sep = "", collapse = ","),
                                      ",test ='", Test, "')", sep = "")
    ))
    rownames(anovas) <- rownames(final_table)
    p_values <- anovas$`Pr(>F)` # Change to ChiSq if more families are added
    p_less <- which(p_values < .05)
    p_value_selected <- ifelse(length(p_less) == 0, 1, max(p_less))
    selected <- rnames[p_value_selected]
    
    anova_data <- as.data.frame(anovas)
    anova_data["Test"] <- c("", "a vs b", "b vs c", "c vs d")[descriptions]
    anova_data <- anova_data %>% 
      mutate(across(all_of(colnames(.)[1:6]), function(x) round(x, 4))) %>% 
      select(1, 2, 3, 4, 7, 5, 6) %>% 
      mutate(across(everything(), ~ifelse(is.na(.x), "", .x)))
    
    # Add signif stars
    stars <- symnum(as.numeric(anova_data$`Pr(>F)`), corr = FALSE, na = FALSE, 
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    symbols = c("***", "**", "*", ".", " "))
    anova_data <- mutate(anova_data, "Pr(>F)" = paste(`Pr(>F)`, format(stars)))
    
    if(isTRUE(verbose)){
      message("Selection using F-tests.\n")
      print(final_table, right = FALSE)
      cat(paste0("\033[0;", 32, "m",
                 "\nFor functional redundancy to be present, none of these tests should be significant.",
                 "\033[0m","\n"))
      print(anova_data)
      sleg <- attr(stars, "legend")
      cat("---\nSignif. codes:  ", sleg, sep = "", fill = getOption("width") + 
            4 + max(nchar(sleg, "bytes") - nchar(sleg)))
      # printCoefmat(anovas, na.print = "")
    }
  } else {
    final_table[, selection] <- sapply(model_list, get(selection))
    final_table$Description <-  desc
    selected <- rnames[which.min(final_table[[selection]])]
    
    if(isTRUE(verbose)){
      message(paste("Selection by ", selection, 
                    sep = ""))
      #final_table <- final_table[, c(1, 4, 2, 3, 5)]
      cat(paste0("\033[0;", 32, "m",
                 paste0("\nFor functional redundancy to be present, `Model a` should have the lowest ", selection, "."),
                 "\033[0m","\n"))
      
      print(final_table, right = FALSE)
      cat(paste0("\033[0;", 31, "m",
                 paste0("\nWarning: Model with the lowest ",
                        selection, " will be selected, even if the difference is very small. Please inspect other models to see differences in ",
                        selection, "."),
                 "\033[0m","\n"))
    }
  }
  
  cat(paste0("\033[0;", 32, "m",paste0("\nSelected model: ", final_table$Description[which(rnames == selected)]),"\033[0m","\n"))
  final_message <- if(selected == "Model a") "There is functional redundancy between the specified species" 
  else "There is no functional redundancy between the specified species"
  cat(paste0("\033[0;", 32, "m",final_message,"\033[0m","\n"))
  # mods <<- model_list
  return(list("selected_model" = model_list[[which(rnames == selected)]],
              "code" = selected))
}

data("Switzerland")
# Use a fitted model object for testing functional redundancy
mod <- DI(y = "yield", prop = c("p1", "p2", "p3", "p4"), 
          DImodel = "FULL", data = Switzerland, estimate_theta = F)

# Test for functional redundancy between legumes and grasses
t1 <- test_functional_redundancy(model = mod, groups = c("Gr", "Gr", "Le1", "Le2"))
# The results of the test are printed at the last line
# The grasses are functionally redundant. 
# The final selected model is returned, thus t1 contains the chosen model object
t1

# Can also test for functional redundancy between legumes
# The results of the test also indicate the reason for the absence of functional redundancy. 
# In this example, since the model a vs b test is significant, it means that the legumes have different 
# identity effects, thus violating condition 3. 
t2 <- test_functional_redundancy(model = mod, groups = c("Gr1", "Gr2", "Le", "Le"))
t2

# It is also possible to fit a model on the fly and test for functional redundancy
# Now we need to specify the data, prop and y arguments, with any additional arguments
# to the DI function
# Again, we tested for functionally redundancy between grasses and the model suggests
# no functional redundancy
t3 <- test_functional_redundancy(data = Switzerland, prop = 4:7, y = "yield",
                                 groups = c("Gr", "Gr", "Le1", "Le2"))

# However, currently only models where theta = 1 can be fit
t4 <- test_functional_redundancy(data = Switzerland, prop = 4:7, y = "yield",
                                 groups = c("Gr", "Gr", "Le1", "Le2"),
                                 theta = 0.5)
# If theta is specified with a value other than 1, it'll be forced to be 1
t4

# Additionally, if a model where theta not equals 1 is specified a error will be thrown
# Model with theta
theta_mod <- DI(y = "yield", prop = c("p1", "p2", "p3", "p4"), 
                DImodel = "FULL", data = Switzerland, 
                estimate_theta = TRUE)
test_functional_redundancy(model = theta_mod, groups = c("Gr", "Gr", "Le1", "Le2"))

# By default, the function will fit the full pairwise model as the complex model (model d) 
# and use that for testing redundancy, i.e., the DImodel parameter will he a default value 
# of "FULL", but the user can specify different custom models, either
# by specifying them in the model paramter or using a different value for the DImodel parameter
# Both examples are shown here

# For further examples we have fixed the groups paramter to be c("Gr", "Gr", "L1", "L2"),
# implying that we are testing for functional redundancy between grasses.
diff_mod <- DI(y = "yield", prop = c("p1", "p2", "p3", "p4"), 
               DImodel = "AV", data = Switzerland)
t5 <- test_functional_redundancy(model = diff_mod, groups = c("Gr", "Gr", "L1", "L2"))

t6 <- test_functional_redundancy(data = Switzerland, prop = 4:7, y = "yield",
                                 groups = c("Gr", "Gr", "L1", "L2"),
                                 DImodel = "ADD")
# Notice that when specifying a custom model, this model might not be hierarchical in nature and
# hence F-test can't be performed, so model comparison is done using information criteria 
# (AICc by default). A different information criteria can be if the user desires by using the selection parameter
t5 <- test_functional_redundancy(model = diff_mod, 
                                 selection = "BIC", 
                                 groups = c("Gr", "Gr", "L1", "L2"))

t6 <- test_functional_redundancy(data = Switzerland, prop = 4:7, y = "yield",
                                 groups = c("Gr", "Gr", "L1", "L2"),
                                 DImodel = "ADD",
                                 selection = "BICc")

# If it is not possible to fit the full model, then I've designed the function to skip
# fitting it and only compare the three heirarchical models (Model a, b and c).
# This cutoff is for datasets with more than 10 species
data("Bell")
# Notice there is no Model d now
t7 <- test_functional_redundancy(data = Bell, prop = 4:75, y = "response",
                                 groups = c(rep("G1", 70), "G2", "G3"),
                                 estimate_theta = F, selection = "F-test")

# If this is not desirable, then the user has options
# They could either specify a different model in the DImodel paramter
# Note: Information criteria will be used for model selection in this case
t8 <- test_functional_redundancy(data = Bell, prop = 4:75, y = "response",
                                 groups = c(rep("G1", 70), "G2", "G3"),
                                 DImodel = "AV",
                                 estimate_theta = F, selection = "F-test")

# Or they could fit a model themselves and pass it to the function
full_mod <- DI(data = Bell, y = "response", prop = 4:75, DImodel = "FULL")
# This is not sensible from an ecological point but just as an example
# It'll take some time to fit this model
t9 <- test_functional_redundancy(model = full_mod, 
                                 groups = c(rep("G1", 70), "G2", "G3"))

# Specify verbose = FALSE to not print the tables
t9 <- test_functional_redundancy(model = full_mod, verbose = FALSE,
                                 groups = c(rep("G1", 70), "G2", "G3"))
