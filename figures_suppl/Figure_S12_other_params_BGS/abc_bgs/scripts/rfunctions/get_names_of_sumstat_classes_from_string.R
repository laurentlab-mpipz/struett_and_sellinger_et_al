
get_names_of_sumstat_classes_from_string = function(my_string, which_arg) {
  # my_string seperates lists (using space as sep) which are separated by "NEXT"
  
  a <- strsplit(my_string, split = "NEXT")[[1]][which_arg]
  a <- strsplit(a, split = " ")[[1]]
  
  return(a)
}
