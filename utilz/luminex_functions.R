## HELPER FUNCTIONS FOR LUMINEX DATA PROCESSING + ANALYSIS ------

library(tidyverse)

read_luminex_data <- function(filepath) {
  # reads in data from luminex experiment, extracting median values, net_mfi, and bead counts
  # in that order. returns as list of three dataframes

  # figure out where we need to start to get the juicy stuff
  all_lines <- readLines(filepath)
  median_start_index <- grep("\"DataType:\",\"Median\"", all_lines)
  net_mfi_start_index <- grep("\"DataType:\",\"Net MFI\"", all_lines)
  count_start_index <- grep("\"DataType:\",\"Count\"", all_lines)
  # not super incredibly sure why -3 is right here but seems legit
  number_lines_in_expt <- net_mfi_start_index - median_start_index - 3

  # read in the juicy stuff as separate data frames
  median <- read_csv(filepath, 
                   skip = median_start_index, n_max = number_lines_in_expt)
  net_mfi <- read_csv(filepath, 
                   skip = net_mfi_start_index, n_max = number_lines_in_expt)
  counts <- read_csv(filepath, 
                   skip = count_start_index, n_max = number_lines_in_expt)

  luminex_data <- list(median, net_mfi, counts)

  luminex_data
}


get_estimated_concentrations_one_analyte <- function (data, top_standard_concentration, sample_label = 'well_id_timepoint', verbose=TRUE) {
  # fits nplr curves 
  # input: luminex data with a column labeled `well_id_timepoint` that has standards labeled "n-standard"
  # and a column labeled 'mfi'

  # step one: get standards + fit standard curve
  standards_summary <- data |>
    filter(str_detect(get(sample_label), 'tandard')) |>
    filter(str_detect(get(sample_label), 'tandard')) |>  # 7 is blank
    mutate(mmfi = mfi)
    # uncomment to aggregate standards before fitting curve
    # should add as a function param. oh well
    #group_by(well_id_timepoint) |>
    #summarize(mmfi = mean(mfi))

  top_standard_mmfi <- max(standards_summary$mmfi)

  standards <- standards_summary |>
    mutate(prop_of_top = mmfi/top_standard_mmfi) |>
    # calc standard concs : top_standard/(3^(n-1)), where n is the dilution number
    # 7 is blank
    mutate(standard_conc = case_when(
      TRUE ~ (top_standard_concentration / (3^(as.double(str_extract(get(sample_label), '\\d')) - 1)))
    ))

  model <- nplr::nplr(x = standards$standard_conc, y = standards$prop_of_top, useLog = TRUE)
  
  if (verbose){
    plot(model)
    print(model)
  }

  bottom_asymptote <- model@pars$bottom
  top_asymptote <- model@pars$top


  # now let's flag all the values that are win 3 std dev of blank
  blank <- data |>
    filter(str_detect(get(sample_label), '7') & str_detect(get(sample_label), 'tandard')) |>
    summarize(mmedian = mean(mfi),
              stdev = sd(mfi))


  too_low_basically_blank <- blank$mmedian + (3 * blank$stdev)

  targets <- data |>
    mutate(prop_of_top = mfi / top_standard_mmfi,
           check_asymptotes = case_when(
            prop_of_top > top_asymptote ~ str_c("> ", as.character(top_standard_concentration)),
            prop_of_top < bottom_asymptote ~ str_c("< ", as.character(top_standard_concentration/(3^6))),
            TRUE ~ 'linear range'
           ))

  nplr::getEstimates(model, targets = targets$prop_of_top) |>
    as.tibble() |>
    select(y, x) |>
    rename('concentration' = 'x') |>
    right_join(y=targets, by=join_by('y' == 'prop_of_top'), relationship = 'many-to-many') |>
    unique() |>
    # NA out the too low ones
    mutate(concentration = case_when(
      mfi < too_low_basically_blank ~ NA,
      TRUE ~ concentration
    )) 
}
# only really need to return well, analyte, and est conc


estimate_concentration_by_analyte <- function(name, data, sample_label) {
  # used in lapply for `estimate_concentrations_of_analytes`
  dataa <- data |>
    filter(analyte == name[1])
  get_estimated_concentrations_one_analyte(dataa, as.numeric(name[2]), sample_label = sample_label)
}


estimate_concentrations_of_analytes <- function(luminex_data, list_of_analytes, sample_label = "well_id_timepoint") {
  # run `get_estimated_concentrations_one_analyte` on a larger dataset with some list of analytes
  output <- lapply(list_of_analytes, estimate_concentration_by_analyte, data = luminex_data, sample_label = sample_label)
  bind_rows(output) 
}
