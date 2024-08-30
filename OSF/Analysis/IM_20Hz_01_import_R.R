#####
# Author: Julius Kricheldorff
# Import of 20Hz data
# Includes task data the Flanker-, Stop-Signal-, Go-NoGo and Simple RT task
# Date 11.01.2023
####

# set wd
setwd('C:/Users/doex9445/Dateien/Julius/20Hz/Data/Original/Daten final/')

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)

# To determine which run was which stimulation condition, we load the UPDRS table

updrs <- read_sav("C:/Users/doex9445/Dateien/Julius/20Hz/Data/Original/Ergebnisse/UPDRS.sav")
updrs <- as_tibble(updrs)

# Participants who were used in the final analysis
part_nrs <- c("001", "002", "004", "005", "006", "007", "008", "009", "010",
              "011", "012", "014", "016", "017", "018", "019", "020")
getwd()

flanker_tb <- tibble(
  Trial_Nr = as.numeric(),
  Required_Response = as.character(),
  nrbeding = as.numeric(),
  Correct_Response = as.numeric(),
  RT = as.numeric(),
  Stimulus_Switch = as.numeric(),
  Response_Switch = as.numeric(),
  Post_Error = as.numeric(),
  eingfehler = as.numeric(),
  Part_nr = as.character(),
  run_nr = as.numeric(),
  UPDRS = as.numeric(),
  Stim = as.numeric(),
  Stim_verb = as.character(),
  Congruency = as.character()
)

GnG_tb <- tibble(
  Trial_Nr = as.numeric(),
  Required_Response = as.character(),
  nrbeding = as.numeric(),
  Correct_Response = as.numeric(),
  RT = as.numeric(),
  Stimulus_Switch = as.numeric(),
  Response_Switch = as.numeric(),
  Post_Error = as.numeric(),
  eingfehler = as.numeric(),
  Part_nr = as.character(),
  run_nr = as.numeric(),
  UPDRS = as.numeric(),
  Stim = as.numeric(),
  Stim_verb = as.character(),
  GoNoGo = as.character()
)

Stopp_tb <- tibble(
  Trial_Nr = as.numeric(),
  Required_Response = as.character(),
  nrbeding = as.numeric(),
  Correct_Response = as.numeric(),
  RT = as.numeric(),
  delay = as.numeric(),
  Stimulus_Switch = as.numeric(),
  #Post_StopTrl = as.numeric(),
  eingfehler = as.numeric(),
  Part_nr = as.character(),
  run_nr = as.numeric(),
  UPDRS = as.numeric(),
  Stim = as.numeric(),
  Stim_verb = as.character(),
  StopTrl = as.character()
)

srt_tb <- 
  tibble(
    Trial_Nr = as.numeric(),
    Required_Response = as.character(),
    nrbeding = as.numeric(),
    Correct_Response = as.numeric(),
    RT = as.numeric(),
    Stimulus_Switch = as.numeric(),
    Post_Error = as.numeric(),
    eingfehler = as.numeric(),
    Part_nr = as.character(),
    run_nr = as.numeric(),
    UPDRS = as.numeric(),
    Stim = as.numeric(),
    Stim_verb = as.character(),
  )

# function to read .dat extension
read_datfile <- function(nr_rows, file_n, skip_s){
  if (is.na(nr_rows)){
    temp_dat <- read.table(file_n, header = FALSE, skip = skip_s)
  } else if (!is.na(nr_rows)){
    temp_dat <- read.table(file_n, header = FALSE, skip = skip_s, nrows = nr_rows)
  }
  
  # load the data
  
  header <- read.table(file_n, nrows = 1, header = FALSE, sep = '', stringsAsFactors = FALSE)
  colnames(temp_dat) <- unlist(header)
  
  return(temp_dat)
}

add_updrs_info <- function(data_f, upd){
  # Also add the stimulation and UPDRS score condition for each run by looking through the UPDRS file
  temp_up <- updrs %>%
    filter(Proband == as.numeric(pt),
           Durchgang == run) %>%
    mutate(Stim_verb = case_when(
      Stim == 1 ~ "130Hz",
      Stim == 2 ~ "20Hz",
      Stim == 3 ~ "OFF"
    ))
  
  #convert to tibble
  temp_tb <- as_tibble(data_f) %>%
    mutate(Part_nr = pt,
           run_nr = run,
           UPDRS = temp_up$UPDRS,
           Stim = temp_up$Stim,
           Stim_verb = temp_up$Stim_verb) 
  return(temp_tb)
}

# the original files are littered, have to read the data in line by line
read_dat_mod <- function(file_name, updrs, max_trial){
  file_size <- readLines(con = file_name)
  hdr <- file_size[1] # get the header
  hdr_rep <- which(file_size == hdr) #number of times the header is found
  skip_var = rep(1, length(hdr_rep))
  # find out how many lines are contaminated by non essential information (sometimes it is one, sometimes two)
  # goal is to find out how many lines to skip
  ind <- 0
  for (breaks in hdr_rep){
    ind <- ind + 1
    # get the header
    header <- read.table(file_name, nrows = 1, skip = breaks - 1, header = FALSE, sep = '', stringsAsFactors = FALSE)
    next_line <- 0
    while (TRUE){
      next_line <- next_line + 1
      line_t <- read.table(file_name, nrows = 1, skip = next_line, header = FALSE, sep = '', stringsAsFactors = FALSE)
      if (length(header)> length(line_t)){
        skip_var[ind] = skip_var[ind] +1
      } else {
        break
      }
    }
  }
  
  
  # first find out if the data is split if not simply load the data
  if (length(hdr_rep) < 2){
    temp_dat <- read_datfile(nr_rows = NA, file_n = file_name, skip_s = skip_var) # load data
    temp_dat <- add_updrs_info(data_f = temp_dat, upd = updrs) # load updrs infos
    return(temp_dat)
  } # when the split data together is smaller than the maximum number of trials you would expect, merge
  else if (max_trial> (length(file_size) - sum(skip_var))){
    temp_one <- read_datfile(nr_rows = hdr_rep[2] - skip_var[1] - 1 , file_n= file_name, skip_s = skip_var[1])
    temp_two <- read_datfile(nr_rows = NA , file_n= file_name, skip_s = hdr_rep[2] + skip_var[2] - 1)
    temp_dat <- bind_rows(temp_one, temp_two)
    temp_dat <- add_updrs_info(data_f = temp_dat, upd = updrs) # load updrs infos
    return(temp_dat)
  } else if ((hdr_rep[2] - (hdr_rep[1] + skip_var[1])) > (length(file_size) - (hdr_rep[2] + 1))){
    # take the first halfs
    temp_dat <- read_datfile(nr_rows = hdr_rep[2] - skip_var[1] - 1 , file_n= file_name, skip_s = skip_var[1])
    temp_dat <- add_updrs_info(data_f = temp_dat, upd = updrs) # load updrs infos
    return(temp_dat)
  } else if ((hdr_rep[2] - (hdr_rep[1] + skip_var)) < (length(file_size) - (hdr_rep[2] + skip_var[2] - 1))){
    # when the split data is larger, take the larger data set  
    temp_dat <- read_datfile(nr_rows = NA , file_n= file_name, skip_s = hdr_rep[2] + skip_var[2] - 1)
    temp_dat <- add_updrs_info(data_f = temp_dat, upd = updrs) # load updrs infos
    return(temp_dat)
    
  }
  return(temp_dat)
}

for (pt in part_nrs){
  # create path string
  path_new <- file.path(getwd(), pt)
  
  # find folders for run one, two and three
  for (run in seq(1,3)){
    r_path <- file.path(path_new, as.character(run))
    
    # find folders for run one, two and three
    files = list.files(path = r_path, pattern = ".dat", full.names = TRUE, recursive = FALSE)
    file_list <- files[!grepl("_test|memory", files)] # sort out the test files and memory files
    print(paste("Loading data of participant", pt, "for run", run, "and found", as.character(length(file_list)), "files"))
    
    # loop through files to load the raw data
    for (fl in file_list){
      # add a participant number and a run number
      if(grepl("flanker", fl)){
        
        temp_tb  <- read_dat_mod(file_name = fl, updrs = updrs, max_trial = 100)
        
        temp_tb <- temp_tb %>%
          rename(RT = zeit,
                 Stimulus_Switch = wechsel,
                 Response_Switch = seitwechsel,
                 Post_Error = postfehler,
                 Correct_Response = richtig,
                 Required_Response = beding,
                 Trial_Nr = nr) %>% # first rename relevant columns to English
          mutate (Stimulus_Switch = recode(Stimulus_Switch, '1' = 0L, '2' = 1L, .default = NULL), # recode variables to more sensible values
                  Post_Error = recode(Post_Error, '1' = 0L, '2' = 1L, .default = NULL),
                  Response_Switch = recode(Response_Switch, '1' = 0L, '2' = 1L, .default = NULL),
                  RT = as.numeric(gsub(",", ".", RT)), # because commas are used instead of dots, we have to replace those first
                  Congruency = case_when(
                    Required_Response == 1 | Required_Response == 2 ~ "congruent",
                    Required_Response == 3 | Required_Response == 4 ~ "incongruent"
                    ),
                  Required_Response = case_when(
                    Required_Response == 1 | Required_Response == 3 ~ "right", 
                    Required_Response == 2 | Required_Response == 4 ~ 'left')
          )
        
        # bind datasets
        flanker_tb <- bind_rows(flanker_tb, temp_tb)
          
        
      } else if (grepl("gng", fl)){
        
        temp_tb  <- read_dat_mod(file_name = fl, updrs = updrs, max_trial = 80)
        
        temp_tb <- temp_tb %>%
          rename(RT = zeit,
                 Stimulus_Switch = wechsel,
                 Response_Switch = seitwechsel,
                 Post_Error = postfehler,
                 Correct_Response = richtig,
                 Required_Response = beding,
                 Trial_Nr = nr) %>% # first rename relevant columns to English
          mutate (Stimulus_Switch = recode(Stimulus_Switch, '1' = 0L, '2' = 1L, .default = NULL), # recode variables to more sensible values
                  Post_Error = recode(Post_Error, '1' = 0L, '2' = 1L, .default = NULL),
                  Response_Switch = recode(Response_Switch, '1' = 0L, '2' = 1L, .default = NULL),
                  RT = as.numeric(gsub(",", ".", RT)), # because commas are used instead of dots, we have to replace those first
                  GoNoGo = case_when(
                    Required_Response == 1 | Required_Response == 2 ~ "Go",
                    Required_Response == 3  ~ "NoGo - Go",
                    Required_Response == 4  ~ "NoGo - Stop"
                  ),
                  Required_Response = case_when(
                    Required_Response == 1 | Required_Response == 3 ~ "right", 
                    Required_Response == 2 | Required_Response == 4 ~ 'left')
          )
        # bind datasets
        GnG_tb <- bind_rows(GnG_tb, temp_tb)
        
        
      } else if (grepl("stopp", fl) & !grepl("1D", fl)){
        temp_tb  <- read_dat_mod(file_name = fl, updrs = updrs, max_trial = 120)
        
        temp_tb <- temp_tb %>%
          rename(RT = zeit,
                 Stimulus_Switch = wechsel,
                 Post_Error = postfehler,
                 Correct_Response = richtig,
                 Required_Response = beding,
                # Post_StopTrl = posthupe,
                 Trial_Nr = nr) %>% # first rename relevant columns to English
          mutate (Stimulus_Switch = recode(Stimulus_Switch, '1' = 0L, '2' = 1L, .default = NULL), # recode variables to more sensible values
                  Post_Error = recode(Post_Error, '1' = 0L, '2' = 1L, .default = NULL),
                  RT = as.numeric(gsub(",", ".", RT)), 
                  #Post_StopTrl = recode(Post_StopTrl, '1' = 0L, '2' = 1L, .default = NULL),
                  StopTrl = case_when(
                    Required_Response == 1 | Required_Response == 2 ~ "Go",
                    Required_Response == 3 | Required_Response == 4  ~ "Stop",
                  ),
                  Required_Response = case_when(
                    Required_Response == 1 | Required_Response == 3 ~ "right", 
                    Required_Response == 2 | Required_Response == 4 ~ 'left')
          )
        # bind datasets
        Stopp_tb <- bind_rows(Stopp_tb, temp_tb)
        
       
        
      } else if (grepl("simplereaktion", fl)){
        #load data
        temp_tb  <- read_dat_mod(file_name = fl, updrs = updrs, max_trial = 80)
        
        temp_tb <- temp_tb %>%
          rename(RT = zeit,
                 Stimulus_Switch = wechsel,
                 Post_Error = postfehler,
                 Correct_Response = richtig,
                 Required_Response = beding,
                 Trial_Nr = counter) %>% # first rename relevant columns to English
          mutate (Stimulus_Switch = recode(Stimulus_Switch, '1' = 0L, '2' = 1L, .default = NULL), # recode variables to more sensible values
                  Required_Response = recode(Required_Response, '1' = 'right', '2' = 'left', .default = NULL),
                  Post_Error = recode(Post_Error, '1' = 0L, '2' = 1L, .default = NULL),
                  RT = as.numeric(gsub(",", ".", RT)), )
        # bind datasets
        srt_tb <- bind_rows(srt_tb, temp_tb)
        
        
      }
      
    }
    
  }
  
}


# NOTES: Paarticipant 1, run 1, signal stopp task, experiment was started twice, only start reading from line 22
# run 1, part 004 flanker, 37 trials extra.
# run 1 part 004 stopp, 20 trials extra,
# run 1 part 004, simplrt 10 trials extra
# run 2 part 004 simplert, only the first 49 trials
# run 3 part 004 simplert, first 25 trials extra

# let us check if we have missing datafiles

flanker_trials <- flanker_tb %>%
  group_by(Part_nr, run_nr) %>%
  summarise(nr_run = n())

# flanker task looks good, only a few trials missing in some conditions

GnG_trials <- GnG_tb %>%
  group_by(Part_nr, run_nr) %>%
  summarise(nr_run = n())

# The first 5 participants have a few trials twice? no idea why that is - not just repeats though performance differs

Stopp_trials <- Stopp_tb %>%
  group_by(Part_nr, run_nr) %>%
  summarise(nr_run = n())

# stop task looks good, only a few trials missing in some conditions

srt_trials <- srt_tb %>%
  group_by(Part_nr, run_nr) %>%
  summarise(nr_run = n())

# srt looks really good, only one participant has missing trials

# Now that we have checked the data, we can save it
write_csv(flanker_tb, file = "C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted/flanker.csv")
write_csv(GnG_tb, 'C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted/GoNoGo.csv' )
write_csv(Stopp_tb, 'C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted/StoppSignal.csv' )
write_csv(srt_tb, 'C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted/SimpleRT.csv' )



