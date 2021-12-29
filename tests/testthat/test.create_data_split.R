context("Test that the create data split function has the expected behavior")

test_that("Generating a data split using the dplyr approach doesn't throw an error or generate duplicates between cal. and val. data", {
  plot<- rep(c("plot1", "plot2", "plot3"), each=42)
  season<- rep(1:6, 21)
  disease<- c(rep(0,84), rep(1,42))
  d<- seq(1:126)
  df <- data.frame(plot,season,disease,d)
  df <- df %>% dplyr::mutate(id=dplyr::row_number())
  
  split_data <- spectratrait::create_data_split(dataset=df, approach="dplyr", 
                                                split_seed=7529075, prop=0.8, 
                                                group_variables=c("plot", 
                                                                "season", 
                                                                "disease"))
  expect_false(sum(split_data$cal_data$id %in% split_data$val_data$id)>0)
})

test_that("Generating a data split using the base approach doesn't throw an error or generate duplicates between cal. and val. data", {
  plot<- rep(c("plot1", "plot2", "plot3"), each=42)
  season<- rep(1:6, 21)
  disease<- c(rep(0,84), rep(1,42))
  d<- seq(1:126)
  df <- data.frame(plot,season,disease,d)
  df <- df %>% dplyr::mutate(id=dplyr::row_number())
  
  split_data <- spectratrait::create_data_split(dataset=df, approach="base", 
                                                split_seed=7529075, prop=0.8, 
                                                group_variables=c("plot", 
                                                                  "season", 
                                                                  "disease"))
  expect_false(sum(split_data$cal_data$id %in% split_data$val_data$id)>0)
})