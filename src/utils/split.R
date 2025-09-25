split_data <- function(df, train_prop, validation_prop, test_prop){
  
  # Sanity check
  if (train_prop > 1 | validation_prop > 1 | test_prop > 1){
    train_prop <- train_prop / 100
    validation_prop <- validation_prop / 100
    test_prop <- test_prop / 100
  }
  if (train_prop + validation_prop + test_prop != 1){
    print("Proportion should sum to 1 or 100 if in % \n
          We do a classical split 80/10/10")
  }
  
  n <- nrow(df) # data frame size
  train_size <- floor(train_prop * n) 
  validation_size <- floor(validation_prop * n) 
  test_size <- n - validation_size - train_size
  
  # Index for split
  trainIndex <- 1:train_size
  validIndex <- (train_size + 1):(train_size + validation_size)
  testIndex <- (train_size + validation_size + 1):n
  
  # Split data set into training, validation and test
  trainData <- df[trainIndex, ]
  validData <- df[validIndex, ]
  testData <- df[testIndex, ]
  

  return(list(trainData = trainData, validData = validData, testData = testData))
}
