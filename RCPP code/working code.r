# Load necessary libraries
library(ape)
library(TreeTools)
library(treebalance)
library(rlist)
library(tictoc)
library(ggplot2)
library(dplyr)
library(fpCompare)
library(tictoc)
library(tidyverse)
library(treesliceR)
library(phytools)
library(future.apply)
library(furrr)
library(purrr)
library(microbenchmark)


read_convert <- function(file){
  
  # Checks the format of file/tree and converts to phylo objext
  suppressWarnings({
    tree <- try(ape::read.tree(file), silent = TRUE) # Newick format
    if ("try-error" %in% class(tree)){
      tree <- try(ape::read.nexus(file), silent = TRUE) # Nexus format
      if ("try-error" %in% class(tree)){
        tree <- try(ape::read.tree(text=file), silent = TRUE) # String in newick format
        if ("try-error" %in% class(tree)){
          if (inherits(file,"phylo")){ # Already phylo object
            tree <- file
          }else if(!(inherits(file,"phylo"))){ # file is none of the above
            return(print("Tree must be in Newick or NEXUS format, or be a phylo object.")) 
          }
        }
      }
    }
  })  
  
  # If no branch lengths, assign all to be one
  if (length(tree$edge.length) == 0){
    tree$edge.length <- rep(1, times = length(tree$edge[,1]))
  }
  
  return(tree)
}

# Function to calculate total abundance descending from each branch of a phylo object,
abundance_phylo <- function(tree, abundance_data) {
  
  node_labels <- append(tree$tip.label, tree$node.label) # all node labels
  
  # Initialize a dictionary to store the total abundance for each branch
  total_abundance <- list()
  
  # Define a recursive function to calculate total abundance
  calculate_abundance <- function(node) {
    
    # Initialize the total abundance for the current node
    node_abundance <- abundance_data[which(abundance_data[,1] == node_labels[node]), 2]
    
    # Loop through each child of the current node
    for (child in tree$edge[tree$edge[, 1] == node, 2]) {
      
      # If the child is not a leaf, recursively calculate the total abundance of its subtree
      if (!(child %in% c(1:length(tree$tip.label)))) {
        child_abundance <- calculate_abundance(child)
        
      }
      # If the child is a leaf, use its abundance as the total abundance of its subtree
      else {
        # Store child/leaf abundance
        child_abundance <- abundance_data[which(abundance_data[,1] == node_labels[child]), 2]
        # Total abundance of leaf is its own abundance
        total_abundance[[as.character(child)]] <<- abundance_data[which(abundance_data[,1] == node_labels[child]), 2]
      }
      
      # Add the total abundance of the child subtree to the current node's total abundance
      node_abundance <- node_abundance + child_abundance
    }
    
    # Store the total abundance for the current node
    total_abundance[[as.character(node)]] <<- node_abundance
    
    # Return the total abundance for the current node
    return(node_abundance)
  }
  
  # Start the recursive calculation from the root node
  calculate_abundance(tree$edge[1, 1])
  
  
  # Return the dictionary of total abundances
  return(total_abundance)
}

# Calculates all ancestor nodes of a given node, including itself
getAllAncestors <- function(tree,node){
  getAncestor <- function(tree, node){ # Returns immediate ancestor of given node
    i <- which(tree$edge[, 2] == node)
    return(tree$edge[i, 1])
  }
  root_node <- length(tree$tip.label) + 1 # Number assigned to root node
  ancestors <- c(node) # Include node as an ancestor of itself
  anc_node <- root_node + 1 # Assign so while loop will start
  if (node != root_node){
    while (anc_node > root_node){
      anc_node <- getAncestor(tree, node) # Ancestor of node
      ancestors <- append(ancestors, anc_node) # Append ancestor to list
      node <- anc_node # Move to ancestor
    }
  }
  return(ancestors)
}

# Calculates all descendant branches of a node, recording distance from node or "x"
get_descendant_branches <- function(tree, node){
  root <- FALSE
  if (node == (length(tree$tip.label)+1)){ # Is the root node
    root <- TRUE
    branches <- which(tree$edge[,1] == node) # Select row index of direct descendant branches
    if (length(branches) == 0) { # If node is leaf return empty matrix
      return(data.frame("Start_node" = numeric(), "End_node" = numeric(),
                        "Branch_length" = numeric(), "x" = numeric()))}
    else{
      descendants <- c() # Empty array to store descendants
      # Select row index in tree$edge for all branches that descend from node
      for (i in 1:length(branches)){
        descendants <- append(descendants, 
                              which(descendant_edges(edge = branches[i],
                                                    tree$edge[,1], 
                                                    tree$edge[,2])))
      }
    }  
  }else{ # Not root node
    # Calculate descendant branches from branch connecting node to parent
    # then delete this branch after
    branch <- which(tree$edge[,2] == node) # Select row index of parent branch
    if (length(which(node == c(1:length(tree$tip.label))))) { # If node is leaf return empty matrix
      return(data.frame("Start_node" = numeric(), "End_node" = numeric(),
                        "Branch_length" = numeric(), "x" = numeric()))}
    else{
      descendants <- which(descendant_edges(tree$edge[,1], tree$edge[,2], edge = branch))
    }
  }
    sub_tree <- tree$edge[descendants,] # Select subtree from node
    x <- sapply(sub_tree[,2], distance, top_node=node, tree=tree) # Record distance between node and end of each branch
    
    # x is the distance from node to the end of the branch 
    # i.e. x and branch length are only equal for direct descendant branches
    df_branch_info <- data.frame("Start_node" = sub_tree[,1], "End_node" = sub_tree[,2],
                                 "Branch_length" = tree$edge.length[descendants], 
                                 "x" = x) # Store information

  if (root == FALSE){ # Delete parent branch
    df_branch_info <- df_branch_info[-1,]
  }
  return(df_branch_info)
}

# Calculates T_i and S_i for a given node i in a phylogenetic tree
compute_T_i_S_i <- function(tree, node, abundances) {
  
  # Choose branches that are immediate descendants of node and record data in dataframe
  df_descen_info <- data.frame("Start_node" = tree$edge[tree$edge[, 1] == node, 1],
                               "End_node" = tree$edge[tree$edge[, 1] == node, 2],
                               "Branch_lengths" = tree$edge.length[which(tree$edge[, 1] == node)])
  
  # Sum abundances of all branches in df_descen_info, delete the branch/es 
  # corresponding with shortest branch length, repeat until all branches
  # have been summed over (i.e. deleted from dataframe)
  Length <- length(unique(df_descen_info[,"Branch_lengths"])) # Preassign length
  df_list <- vector(mode = "list", length = Length) # List for dataframes
  T_i <- 0 # Initialise T_i sum
  prev_x <- 0 # Keep track of previous value of x
  # Sum over each unique branch length
  for (j in 1:Length) {
    abundance_sum <- 0 # Initialise sum
    if (nrow(df_descen_info) != 0){ # Checks if any branches are left
      for (i in 1:nrow(df_descen_info)) { # Sum over all branches in current region
        if (!(is.na(df_descen_info["End_node"][i,]))){
          end_node <- df_descen_info["End_node"][i,] # Select end node of current branch
          abundance_sum <- abundance_sum + abundances[[as.character(end_node)]] # Sum abundances
        }
      }
    }
    
    # Record all info and delete any branches "passed"
    if (nrow(df_descen_info) != 0){ # If there are branches left
      # Append abundance and corresponding value of x
      df_list[j] <- list(data.frame("S_i" = abundance_sum,
                                    "x" = min(df_descen_info["Branch_lengths"])))
      
      T_i <- T_i + (abundance_sum * (min(df_descen_info["Branch_lengths"]) - prev_x)) # Sum T_i
      prev_x <- min(df_descen_info["Branch_lengths"]) # Update previous x
      
      df_descen_info <- df_descen_info[-(which(df_descen_info["Branch_lengths"] %==% min(df_descen_info["Branch_lengths"]))),] # Remove branch/s corresponding to x value
    }
  }
  # Combine dataframes
  DF_S_i <- Reduce(rbind, df_list)
  return(list(T_i, DF_S_i))
}

# Calculate S_i_a for a given node i and ancestor a in a phylogenetic tree
calculate_S_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i) {
  
  temp <- function(a,b){return(a[[b]])}
  
  # Select all descendant branches from ancestor
  if (node != curr_ancestor){ # If not considering itself as ancestor  
    df_node_info <- get_descendant_branches(tree, curr_ancestor) # Get descendant branches from ancestor
    df_node_info <- df_node_info[df_node_info$x > h,] # Select branches that pass node
    
    # Want all branch lengths and x values to be measure from node, not ancestor
    # Select index values of nodes that have descendant branches that start before
    # given node, then adjust x and branch lengths accordingly
    curr_index <- which(df_node_info[,"x"] - df_node_info[,"Branch_length"] < h)
    # Correct branch lengths
    df_node_info[curr_index,"Branch_length"] <- df_node_info[curr_index,"x"] - h
    # Correct x
    df_node_info$x <- df_node_info$x - h
  }else{ # Considering itself as ancestor
    df_node_info <- get_descendant_branches(tree, node) # No corrections needed
  }
  
  df_node_info_ed <- df_node_info # Create dataframe to delete values from
  
  # A branch is in the current region iff its x and branch length are equal, 
  # otherwise it doesn't start at the beginning of current region.
  # Sum over all branches whose x and branch length are equal, deleting the 
  # branch/es corresponding with smallest branch length/x, then transform
  # all x and branch lengths by - length of previous region 
  # i.e. shifting then to be measured from the started of new region
  Length <- length(df_node_info_ed$Branch_length)
  df_list <- vector(mode = "list", length = Length) # List for dataframes
  abund_list <- vector(mode = "list", length = Length) # Empty dictionary for abundances
  x <- 0 # Keep track of distance from node
  # Sum over each branch as all could have different branch lengths
  for (j in 1:Length) {
    # Only do until reach the end of nodes longest direct descendant branch, l_i
    if (x != l_i){
      # Row index of branches where x == branch length
      index_curr_branches <- which(df_node_info_ed[,"x"] %==% df_node_info_ed[,"Branch_length"])
      
      if (nrow(df_node_info_ed) != 0){ # If there are branches left to sum over
        # Store all branch abundances present for this value of x
        abund_list[[j]] <- sapply(as.character(end_node <- df_node_info_ed[index_curr_branches,"End_node"]),
                                                temp, a = abundances)
        
        abundance_sum <- sum(abund_list[[j]]) # Sum abundances
        
        # Append dataframe to list
        df_list[j] <- list(data.frame("S_i" = abundance_sum,
                                      "x" = (min(df_node_info_ed[index_curr_branches,]["x"]) + x)))
        
        prev_x <- min(df_node_info_ed[index_curr_branches,]["x"]) # Update previous x
        x <- x + prev_x # Update x
        
        # Indices of rows to be deleted 
        index_to_be_deleted <- which(df_node_info_ed[,"x"] %==% min(df_node_info_ed[index_curr_branches,]["x"]))
        df_node_info_ed$x <- df_node_info_ed$x - prev_x # "Reset" x = 0 level
        # Shorten current branch lengths by previous x
        df_node_info_ed[index_curr_branches,]["Branch_length"] <- df_node_info_ed[index_curr_branches,]["Branch_length"] - prev_x
        # Remove branch/s
        df_node_info_ed <- df_node_info_ed[-index_to_be_deleted,]
        
      }
    }else{
      break
    }
  }
  # Delete unused entries
  df_list <- Filter(Negate(is.null), df_list)
  abund_list <- Filter(Negate(is.null), abund_list)
  # Combine dataframes
  DF_S_i <- Reduce(rbind, df_list)
  return(list(DF_S_i, abund_list))
}

# Calculates E/J/M or all given i and a
calculate_EJM_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i, 
                              index_letter, q, individual){
  
  term <- function(a,b,log_base){-(a/b) * log(a/b, base = log_base)}
  
  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h, l_i) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list

  if (individual == FALSE){ # Create empty dataframe/s
    df_E_i_a <- data.frame("E_i_a" = numeric(), "x" = numeric()) 
    df_J_i_a <- data.frame("J_i_a" = numeric(), "x" = numeric())
    df_M_i_a <- data.frame("M_i_a" = numeric(), "x" = numeric())
  }else if (individual == TRUE){
    df_In_i_a <- data.frame("E_i_a" = numeric(), "x" = numeric()) 
  }


  for (k in 1:length(abund_list)){
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    abund_vec <- unlist(abund_list[k], use.names = FALSE)
    
    if (individual == FALSE){
      m <- log(length(abund_vec)) # Calculate out-degree term
      e <- sum(sapply(abund_vec, term, b = S, log_base=exp(1))) # Calculate diversity term
      if (length(abund_vec) > 1 ){ # Removes case of only one branch in region
        j <- sum(sapply(abund_vec, term, b = S, log_base=length(abund_vec))) # Calculate balance term
      }else if (length(abund_vec) == 1){
        j <- 1
      }
      # Store values and corresponding x
      df_M_i_a <- rbind(df_M_i_a, data.frame("M_i_a" = m, "x" = df_S_i_a[k, "x"]))
      df_E_i_a <- rbind(df_E_i_a, data.frame("E_i_a" = e, "x" = df_S_i_a[k, "x"]))
      df_J_i_a <- rbind(df_J_i_a, data.frame("J_i_a" = j, "x" = df_S_i_a[k, "x"]))
    }else if (individual == TRUE){
      if ((index_letter == "D")&(q == 0)){
        ind_val <- log(length(abund_vec))
      }else if ((index_letter == "D")&(q == 1)){
        ind_val <- sum(sapply(abund_vec, term, b = S, log_base=exp(1)))
      }else if (index_letter == "J"){
        if (length(abund_vec) > 1 ){ # Removes case of only one branch in region
          ind_val <- sum(sapply(abund_vec, term, b = S, log_base=length(abund_vec))) # Calculate balance term
        }else if (length(abund_vec) == 1){
          ind_val <- 1
        }
      }
      df_In_i_a <- rbind(df_In_i_a, data.frame("In_i_a" = ind_val, "x" = df_S_i_a[k, "x"]))
    }
  }
  if (individual == FALSE){
    return(list("1DN" = df_E_i_a, "1JN" = df_J_i_a, "0DN" = df_M_i_a))
  }else if (individual == TRUE){
    return(df_In_i_a)
  }
}

# Calculate S_i_a by first attaching all of trees branches to root node/creates star tree
calculate_S_i_a_star <- function(tree, node, abundance_data) {
  
  temp <- function(a,b){return(a[[b]])}
  
  abundances <- abundance_phylo(tree, abundance_data) # Node/branch abundances
  
  df_node_info <- get_descendant_branches(tree, node) # Get descendant branches from ancestor
  # Set all x to be equal to branch length, i.e. all branches are attached to
  # root node
  df_node_info$x <- df_node_info$Branch_length
  df_node_info_ed <- df_node_info # Create node to delete data from
  
  # Sum abundances of all branches in df_node_info_ed, delete the branch/es 
  # corresponding with shortest branch length/x, repeat until all branches
  # have been summed over (i.e. deleted from dataframe)
  Length <- length(unique(df_node_info$x)) # Preassign length
  df_list <- vector(mode = "list", length = Length) # List for dataframes
  abund_list <- vector(mode = "list", length = Length) # Empty dictionary for abundances
  names(abund_list) <- as.character(c(1:Length)) # Change names to characters
  # Sum over each unique branch length
  for (j in 1:Length) {
    if (nrow(df_node_info_ed) != 0){ # Checks if any branches are left
      end_nodes <- df_node_info_ed$End_node # End nodes of current branches
      # Store all branch abundances present for this value of x
      abund_list[[as.character(j)]] <- sapply(as.character(end_nodes), temp, a = abundances)
      
      abundance_sum <- sum(abund_list[[as.character(j)]]) # Sum abundances 
      
      # Append abundance and corresponding value of x
      df_list[j] <- list(data.frame("S_i" = abundance_sum,
                                    "x" = (min(df_node_info_ed$x))))
      # Remove branch/s corresponding to x value
      df_node_info_ed <- df_node_info_ed[-(which(df_node_info_ed$x %==% min(df_node_info_ed$x))),]
    }
  }
  # Combine dataframes
  DF_S_i <- Reduce(rbind, df_list)
  return(list(DF_S_i, abund_list))
}

# Calculates integral for each ancestor
calculate_integral <- function(tree, node, curr_ancestor, S_i, abundances, index_letter,
                               q, individual){
  
  # Set l_i, longest direct descendant branch of node
  if (length(S_i["x"] > 1)){
    l_i <- max(S_i["x"]) # Select longest branch length of a immediate descendant
  }else{
    l_i <- S_i["x"] # If only one value, select that
  }
  
  # If node is a leaf/has size zero
  if (l_i == 0){
    if (individual == TRUE){
      return(0)
    }else{
      return(c(0,0,0))
    }
  }
  
  h <- distance(tree,curr_ancestor,node) # Distance between current node and ancestor
  
  # Assign distance to parent
  if (curr_ancestor == (length(tree$tip.label) + 1)){ # If ancestor is the root node
    d_parent <- l_i  # By definition its Inf, but integral is only nonzero to l_i
  }else{
    d_parent <- tree$edge.length[which(tree$edge[,2] == curr_ancestor)] + h
  }
  
  # Ancestor integral
  # Calculates integral as sums of areas
  int_2 <- 0 # Initialise sum
  prev_x <- h # Keep track of x, h as integral is from h
  if (l_i <= h){ # Nodes furthest away child is closer than distance to ancestor
    int_2 <- 0
    if (individual == TRUE){
      return(0)
    }else{
      return(c(0,0,0))
    }
  }else{
    current_row <- which(S_i[,"x"] > h)[1] # Start sum over first x that "reaches" ancestor
    while (current_row <= length(S_i[,"x"])) { # Sum over all rows of S_i
      x_s <- S_i[current_row,2] # Current value of x for S_i
      value_s <- S_i[current_row,1] # Current value of S_i for given x
      if (x_s > d_parent){ # Region reaches past ancestor's parent
        int_2 <- int_2 + (value_s * (d_parent - prev_x)) # Integral
        prev_x <- x_s # Update previous x
        break # Exit while loop
      }else { # x_s < d_parent
        int_2 <- int_2 + (value_s * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row <- current_row + 1 # Move to next row
      }
    }
  }
  
  # Function for index integral
  # Calculates integral as sums of areas
  integral <- function(df){
    prev_x <- 0 # Keep track of x
    current_row_s <- 1 # Start sum over rows of S_i
    current_row_i <- 1 # Start sum over rows of index value
    int_1 <- 0 # Initialise integral sum
    
    # qD/1J integral
    # Sum over all rows of S_i
    while ((current_row_s <= length(S_i[,"x"]))&(current_row_i <= length(df$x))){
      x_s <- S_i[current_row_s,2] # Current value of x for S_i
      x_i <- df[current_row_i, 2] # Current value of x for qD/1J
      value_s <- S_i[current_row_s,1] # Current value of S_i for given x
      value_i <- df[current_row_i,1] # Current value of qD/1J for given x
      if (x_s < x_i){
        int_1 <- int_1 + (value_s * value_i * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row_s <- current_row_s + 1 # Move to next row
      }else if (isTRUE(all.equal(x_s, x_i))){
        int_1 <- int_1 + (value_s * value_i * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row_s <- current_row_s + 1 # Move to next row
        current_row_i <- current_row_i + 1 # Move to next row
      }else{ # x_i > x_s
        if (x_i < l_i){
          int_1 <- int_1 + (value_s * value_i * (x_i - prev_x)) # Sum integral
          prev_x <- x_i # Update previous x
        }else{ # x_i > l_i
          int_1 <- int_1 + (value_s * value_i * (l_i - prev_x)) # Sum integral
          prev_x <- x_i # Update previous x
        }
        current_row_i <- current_row_i + 1 # Move to next row
      }
    }
    return(int_1) 
  }
  
  if (int_2 != 0){ # Ancestor integral is not zero, calculate other integral/s
    Sum_i_a <- calculate_DJ_i_a(tree, node, abundances, curr_ancestor, h, l_i, 
                                index_letter, q, individual)
    
    if (individual == TRUE){
      int <- integral(Sum_i_a)*int_2
    }else{
      int <- sapply(Sum_i_a, integral)*int_2
    }
  }else if (int_2 == 0){ # Ancestor integral is zero, return zero
    if (individual == TRUE){
      int <- 0
    }else{
      int <- c(0,0,0)
    }
  }
  
  return(int)
}

# Calculates E_i, J_i and M_i 
calculate_EJM_i <- function(tree, node, abundances, index_letter, q, individual){
  
  ancestors <- getAllAncestors(tree,node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe

  EJM_i <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                  abundances=abundances, index_letter = index_letter, q = q,
                  individual=individual)

  
  if (individual == FALSE){
    ejm_i <- apply(EJM_i,1,sum)
    if (T_i != 0){ # h > 0
      E_i <- (1/T_i) * ejm_i[1]
      J_i <- (1/T_i) *ejm_i[2]
      M_i <- (1/T_i) *ejm_i[3]
    }else if (T_i == 0){ # h = 0
      E_i <- 0
      J_i <- 1
      M_i <- 0
    }
    v <- c(E_i, J_i, M_i) 
  }else if (individual == TRUE){
    ejm_i <- sum(EJM_i)
    if (T_i != 0){ # h > 0
      v <- (1/T_i) * ejm_i
    }else if ((T_i == 0)&(index_letter == "J")){ # h = 0
      v <- 1
    }else if (T_i == 0){
      v <- 0
    }
  }

  
  return(v)
}

calculate_DJ_i <- function(tree, node, abundances, index_letter, q, individual){
  
  ancestors <- c(TreeTools::ListAncestors(tree$edge[,1], tree$edge[,2], node), node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe
  
  if (T_i != 0){ # Has descendant branch with branch length greater than zero
    # Calculate integral value/s
    DJ_i <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                  abundances=abundances, index_letter = index_letter, q = q,
                  individual=individual)
    
    # Add all ancestor contributions and normalise
    if (individual == FALSE){
      dj_i <- apply(DJ_i,1,sum)
      D1_i <- (1/T_i) * dj_i[1]
      J1_i <- (1/T_i) *dj_i[2]
      D0_i <- (1/T_i) *dj_i[3]
      v <- c(D1_i, J1_i, D0_i) 
    }else if (individual == TRUE){
      dj_i <- sum(DJ_i)
      v <- (1/T_i) * dj_i
    }
  }else{
    if (individual == FALSE){
      v <- c(0, 0, 0) 
    }else if (individual == TRUE){
      v <- 0
    }
  }
  
  return(v)
}

calculate_DJ_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i, 
                             index_letter, q, individual) {

  index_letter <- toupper(index_letter) # Capitalise input

  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h, l_i) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list

  # Function to calculate the index/indices sum/s
  sum_function <- function(abund_vec, S, individual, index_letter){
    if (individual == TRUE) { # Calculating one index
      # Checks desired index and calculates value
      if (index_letter == "J") { # J
        # J is defined to be 1 when only one branch is present
        if (length(abund_vec) != 1) { # There is more than one branch in region
          ind_val <- sum(-(abund_vec / S) * log((abund_vec / S),
                                                base=length(abund_vec)))
        }else if (length(abund_vec) == 1) { # There is one branch in region
          ind_val <- 1
        }
      }else { # D
        # Checks q and calculates corresponding index
        if (q == 1) {
          ind_val <- sum(-(abund_vec / S) * log((abund_vec / S)))
        }else if (q == 0) {
          ind_val <- log(length(abund_vec))
        }
      }
      return(ind_val)
    }else if (individual == FALSE) { # Calculate all indices
      D1 <- sum(-(abund_vec / S) * log((abund_vec / S)))
      D0 <- log(length(abund_vec))
      if (length(abund_vec) != 1) { # There is more than one branch in region
        J1 <- sum(-(abund_vec / S) * log((abund_vec / S),
                                         base=length(abund_vec)))
      }else if (length(abund_vec) == 1) { # There is one branch in region
        J1 <- 1
      }
      return(c(D0, D1, J1))
    }
  }
  # Calculate index values
  values <- mapply(sum_function, abund_vec = abund_list, S = df_S_i_a$S_i, 
                   individual=individual, index_letter=index_letter)

  # Return either a dictionary of dataframes with index values or 
  # a single dataframe with desired index values
  if (individual == FALSE) { # All indices
    return(list("1DN" = data.frame("D1" = values[2,], "x" = df_S_i_a$x),
                "1JN" = data.frame("J1" = values[3,], "x" = df_S_i_a$x),
                "0DN" = data.frame("D0" = values[1,], "x" = df_S_i_a$x)))
  }else if (individual == TRUE) { # Single index
    return(data.frame("In" = values, "x" = df_S_i_a$x))
  }
}

descendant_edges <- function(parent, child, edge,
                            nEdge = length(parent)) {

  ret <- logical(nEdge)
  edgeSister <- match(parent[edge], parent[-edge])
  if ((!(is.na(edgeSister)))&(edgeSister >= edge)) {
    # Added check for linearity in tree
    # edgeSister is really 1 higher than you think, because we knocked out
    # edge "edge" in the match
    ret[edge:edgeSister] <- TRUE
        
    # Return:
    ret
  } else {
    nextEdge <- edge
    revParent <- rev(parent)
    repeat {
      if (revDescendant <- match(child[nextEdge], revParent, nomatch=FALSE)) {
        nextEdge <- 1 + nEdge - revDescendant
      } else break;
    }
    ret[edge:nextEdge] <- TRUE
        
    # Return:
    ret
    }
}

# Calculates longitudinal or star mean for D for q=0 or 1, or J for q =1, 
# or all for one mean type
long_star <- function(file, node_abundances = FALSE, mean_type, index_letter = "D", q = 1,
                      individual = FALSE){
  
  index_letter <- toupper(index_letter) # Capitalise input
  mean_type <- toupper(mean_type) # Capitalise input
  
  tree <- read_convert(file) # Convert tree to phylo object

  # Check if tree is linear
  if ((length(tree$tip.label) == 1)&(mean_type == "LONGITUDINAL")){
    # If tree is linear all longitudinal indices are 1
    if (individual == TRUE){
      index <- 1
      return(index)
    }else{
      List <- list("D0L"= 1,"D1L" = 1,"J1L" = 1)
      return(List)
    }
  }else{
    # If tree is linear but star mean is desired
    # Need to form tree properly
    if ((length(tree$tip.label) == 1)&(mean_type == "STAR")){
      if (inherits(file, "phylo") == FALSE){ # Check if phylo object given
      # If phylo object not given, form tree proprly
        tree <- form_linear_phylo(file)
      }
      # If no abundances given, it will be assumed leaf has size one
      # and index values are trivial
      # Also if tree is only root node and leaf,
      # abundances do not matter and index values are again trivial
      if ((!is.data.frame(node_abundances))){
        if (individual == TRUE){
          index <- 1
          return(index)
        }else{
          List <- list("D0S"= tree$Nnode,"D1S" = tree$Nnode,"J1S" = 1)
          return(List)
        }
      }
    }
  }

  node <- tree$edge[1,1] # Select root node
  
  # Checks if tree has abundance data
  # If it does it calculates branch/node abundance data.
  # If it doesn't, it assign leaves to be equally abundant and internal nodes 
  # to have size zero
  if (is.data.frame(node_abundances)) { # Tree has abundance data
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }else if (!(is.data.frame(node_abundances))) { # Tree doesn't have abundance data
    num_tips <- tree$edge[1,1] - 1 # Number of tips
    tree$tip.label<- as.character(c(1:num_tips)) # Assign node labels
    num_nodes <- tree$Nnode # Number of nodes
    # Assign node labels
    tree$node.label <- as.character(c((num_tips + 1):(num_tips + num_nodes)))

    # Create abundance dataframe
    node_abundances <- data.frame("names" = c(tree$node.label, tree$tip.label),
                                  "values" = rep(c(0, (1/num_tips)), times=c(tree$Nnode, num_tips)))
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }

  # Selects mean type and runs corresponding function for S_i_a
  if (mean_type == "LONGITUDINAL") {
    S_i_a_res <- calculate_S_i_a(tree, node, abundances, node, 0, sum(tree$edge.length)) # Run function
  }else if (mean_type == "STAR") {
    S_i_a_res <- calculate_S_i_a_star(tree, node, node_abundances) # Run function
  }

  # Calculate index value/s
  # For each region of x in df_S_i_a, calculates index/indices using
  # corresponding abundance list, this contains every branch abundance
  # in this region, index is calculated and summed over every region of x
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list
  # Function to calculate index/indices sum/s
  sum_function <- function(abund_vec, x, S, individual, index_letter){
    # Calculates index values
    if (individual == TRUE) { # One index
      if (index_letter == "J") { # J for q = 1
        if (length(abund_vec) != 1) { # More than one branch in region
          h <- (sum(-(abund_vec) * log((abund_vec / S),
           base = length(abund_vec))) * x)
        }else if (length(abund_vec) == 1) { # Only one branch in region
          h <- (1 * S * x)
        }
      }else{ # D
        # Check desired index
        if (q == 1) {
          h <- (sum(-(abund_vec) * log((abund_vec / S))) * x)
        }else if (q == 0) {
          h <- (S * x * log(length(abund_vec), base = exp(1)))
        }
      }
      return(h)
    }else if (individual == FALSE) { # All indices
      h1 <- (sum(-(abund_vec) * log((abund_vec / S))) * x)
      h0 <- (S * x * log(length(abund_vec), base = exp(1)))
      if (length(abund_vec) != 1) { # More than one branch in region
        j1 <- (sum(-(abund_vec) * log((abund_vec / S),
         base=length(abund_vec))) * x)
      }else if (length(abund_vec) == 1) { # Only one branch in region
        j1 <- (1 * S * x)
      }
      return(c(h0, h1, j1))
    }
  }

  # Calculate size of regions
  X <- append(df_S_i_a$x[1], diff(df_S_i_a$x))
  # Calculate index values
  values <- mapply(sum_function, abund_vec = abund_list, S = df_S_i_a$S_i,
                   x = X, individual = individual, index_letter = index_letter)
  # Calculate normalisation term
  T_S_sum <- sum(df_S_i_a$S_i * X)

  # Normalise index/indices
  if (individual == TRUE) { # Single index
    if (index_letter == "J") { # J
      H <- (sum(values) / T_S_sum)
    }else { # D
      H <- exp((sum(values) / T_S_sum))
    }
  }else if (individual == FALSE) { # All indices
    if (mean_type == "STAR") { # Star mean
      H <- list("D0S" = exp((sum(values[1,]) / T_S_sum)),
                "D1S" = exp((sum(values[2,]) / T_S_sum)),
                "J1S" = sum(values[3, ]) / T_S_sum)
    }else if(mean_type == "LONGITUDINAL") { # Longitudinal mean
      H <- list("D0L" = exp((sum(values[1,]) / T_S_sum)),
                "D1L" = exp((sum(values[2,]) / T_S_sum)),
                "J1L" = sum(values[3,]) / T_S_sum)
    }
  }

  return(H)
}

distance <- function(tree, top_node, bottom_node){
  curr_node <- bottom_node # Select bottom node as starting node
  dist <- 0 # Initialise distance sum
  while (curr_node != top_node) { # Do until reach desired top node
    # Each node only has one parent branch,
    # select that branch, add the distance between node and parent,
    # move to parent node and repeat
    dist <- dist + tree$edge.length[which(tree$edge[,2] == curr_node)] # Sum distance
    curr_node <- tree$edge[which(tree$edge[,2] == curr_node),1] # Move to parent node
  }
  return(dist)
}

# Calculates node-wise mean for D for q = 0 or 1, or for J for q = 1, or all
node <- function(file, node_abundances = FALSE, index_letter = "D", q = 1,
                 individual = FALSE){
  
  tree <- read_convert(file) # Capitalise inputs
  index_letter <- toupper(index_letter) # Capitalise inputs

  # Checks if the tree is linear
  # If tree is linear all node indices are 1
  if (length(tree$tip.label) == 1){
    if (individual == TRUE){
      index <- 1
      return(index)
    }else{
      
      List <- list("D1N"= 1,"J1N" = 1,"D0N" = 1)
      return(List)
    }
  }else{
      # Checks if tree has abundance data
  # If it does it calculates branch/node abundance data.
  # If it doesn't, it assign leaves to be equally abundant and internal nodes 
  # to have size zero
  if (is.data.frame(node_abundances)){ # Tree has abundance data
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }else if (!(is.data.frame(node_abundances))){ # Tree doesn't have abundance data
    num_tips <- tree$edge[1,1] - 1 # Number of tips
    tree$tip.label<- as.character(c(1:num_tips)) # Assign node labels
    num_nodes <- tree$Nnode # Number of nodes
    # Assign tip labels
    tree$node.label<- as.character(c((num_tips+1):(num_tips + num_nodes)))
    
    # Create abundance dataframe
    node_abundances <- data.frame("names" = c(tree$node.label, tree$tip.label),
                                  "values" = rep(c(0, (1/num_tips)), times=c(tree$Nnode, num_tips)))
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }
  
  # Calculates h_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }
  
  # Numbers of every internal node, corresonding to numbers in phylo object
  nodes <- c((length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode))
  
  # Calculates node-averaged indice/s
  if (individual == FALSE){# All indices
    # Calculates the indices values for each node
    DJ <- sapply(nodes, calculate_DJ_i, tree=tree, abundances = abundances, 
                 individual = individual, q = q, index_letter = "D")
    # Normalise the indices
    D1N <- (1/T)*sum(DJ[1,])
    J1N <- (1/T)*sum(DJ[2,])
    D0N <- (1/T)*sum(DJ[3,])
    # List of index values
    List <- list("D1N"= exp(D1N),"J1N" = J1N,"D0N" = exp(D0N))
    return(List)
  }else if (individual == TRUE){ # One index
    # Calculates the index value for each node
    DJ <- sapply(nodes, calculate_DJ_i, tree=tree, abundances = abundances, 
                 individual = individual, index_letter = index_letter, q = q)
    # Normalise
    if (index_letter == "J"){ # Index J
      index <- (1/T)*sum(DJ)
    }else if (!(index_letter == "J")){ # Index D
      index <- exp((1/T)*sum(DJ))
    }
    
    return(index)
  }
  }
}

# Returns a list of all index values, key is letter for index with first number second
# e.g. 1JN has key J1N
all_indices <- function(file, node_abundances = FALSE){
  # Calculates all indices
  node <- node(file, node_abundances, "D", 1, FALSE)
  star <- long_star(file, node_abundances, "Star", "D", 0, FALSE)
  long <- long_star(file, node_abundances, "Longitudinal", "D", 0, FALSE)
  
  # List of each index value
  values <- list("D0N"= node$D0N,"D1N" = node$D1N,"J1N" = node$J1N, "D0S" = star$D0S,
                 "D1S" = star$D1S, "J1S" = star$J1S, "D0L" = long$D0L, "D1L" = long$D1L,
                 "J1L" = long$J1L)
  return(values)
}























# Set the directory for the tree files
tree_directory <- "C:\\Users\\chiti\\Desktop\\coding\\PhD work\\tree files\\new method trees\\"

# List all tree files in the directory (supporting .nexus, .nwk, .tree, etc.)
tree_files <- list.files(tree_directory, full.names = TRUE, pattern = "\\.(nexus|nwk|tree|tre|txt)$")

# Check if any files were found
if (length(tree_files) == 0) {
  stop("No tree files found in the specified directory.")
}

# Helper function to safely extract index values as numeric
safe_extract <- function(index) {
  if (is.null(index) || length(index) == 0) {
    return(NA)  # Return NA if index is NULL or empty
  } else {
    return(as.numeric(index))  # Return the index value as numeric
  }
}

# Initialize a global results dataframe
results_df <- data.frame(RootLabel = character(), LeafCount = numeric(), D0N = numeric(), 
                         D1N = numeric(), J1N = numeric(), D0S = numeric(), D1S = numeric(), 
                         J1S = numeric(), D0L = numeric(), D1L = numeric(), 
                         J1L = numeric(), TimeTaken = numeric(), stringsAsFactors = FALSE)

# Loop over all files and calculate indices
for (tree_file in tree_files) {
  cat("Processing file:", basename(tree_file), "\n")
  
  # Read the tree from the file
  tree <- tryCatch(read_convert(tree_file), 
                   error = function(e) {
                     warning(paste("Error reading file:", tree_file, "\n", e))
                     return(NULL)
                   })
  
  if (is.null(tree)) next
  
  # Extract subtrees
  all_subtrees <- subtrees(tree)
  
  # Ensure subtrees are valid phylo objects
  valid_subtrees <- lapply(all_subtrees, function(st) if (inherits(st, "phylo") && Ntip(st) <= 200) st else NULL)
  valid_subtrees <- valid_subtrees[!sapply(valid_subtrees, is.null)]
  
  # Sort valid subtrees by number of leaves
  valid_subtrees <- valid_subtrees[order(sapply(valid_subtrees, Ntip))]
  
  # Create a temporary dataframe for the results of this specific tree
  tree_results_df <- data.frame(RootLabel = character(), LeafCount = numeric(), D0N = numeric(), 
                                D1N = numeric(), J1N = numeric(), D0S = numeric(), D1S = numeric(), 
                                J1S = numeric(), D0L = numeric(), D1L = numeric(), 
                                J1L = numeric(), TimeTaken = numeric(), stringsAsFactors = FALSE)
  
  # Process each valid subtree
  for (sub_tree in valid_subtrees) {
    
    # Extract the root label of the subtree
    root_label <- if (!is.null(sub_tree$node.label) && length(sub_tree$node.label) > 0) {
      as.character(sub_tree$node.label[1])  # Ensure RootLabel is always a character
    } else {
      NA  # If no label is available, use NA
    }
    
    # Calculate indices for the subtree
    indices_time <- system.time({
      indices <- tryCatch(all_indices(sub_tree), 
                          error = function(e) {
                            warning(paste("Error calculating indices for a subtree:", e))
                            return(NULL)
                          })
    })
    
    if (is.null(indices)) next
    
    # Record time taken for the index calculation
    index_calculation_time <- indices_time["elapsed"]

    # Extract index values safely
    index_values <- map_dbl(c("D0N", "D1N", "J1N", "D0S", "D1S", "J1S", "D0L", "D1L", "J1L"),
                            ~ safe_extract(indices[[.]]))
    
    # Append the results to the temporary dataframe
    temp_df <- tibble(
      RootLabel = root_label,  # Root label of the subtree
      LeafCount = Ntip(sub_tree), # Number of tips in the subtree
      D0N = index_values[1],
      D1N = index_values[2],
      J1N = index_values[3],
      D0S = index_values[4],
      D1S = index_values[5],
      J1S = index_values[6],
      D0L = index_values[7],
      D1L = index_values[8],
      J1L = index_values[9],
      TimeTaken = index_calculation_time  # Elapsed time
    )
    tree_results_df <- bind_rows(tree_results_df, temp_df)
  }
  
  # Save results for this tree
  output_file <- paste0("results for ", gsub("\\.(nexus|nwk|tree|tre|txt)$", "", basename(tree_file)), " 100.csv")
  output_file_path <- file.path(dirname(tree_file), output_file)
  
  if (nrow(tree_results_df) > 0) {
    write.csv(tree_results_df, file = output_file_path, row.names = FALSE)
    print(paste("Results for", basename(tree_file), "exported successfully."))
  } else {
    warning(paste("No results to export for", basename(tree_file)))
  }
  
  # Append to global results
  tree_results_df$RootLabel <- as.character(tree_results_df$RootLabel)  # Ensure consistent type
  results_df <- bind_rows(results_df, tree_results_df)
}

print("Indices calculated and exported for all trees.")