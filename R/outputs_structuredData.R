#' Outputs for structuredData class
#' 
#' @export

setClassUnion("listorNULL", c('list','NULL'))

#'bru_sdm_data class
#'
#' @export

setClass('structuredData', 
         slots = c(dataPO = 'listorNULL',
                   dataPA = 'listorNULL'))


#'Print methods for structuredData object
#'
#' @exportS3Method
#' 

setMethod('show', 'structuredData',
          function(object) {
            
  cat('Summary output for structuredData object:')
  cat('\n\n')
  
  if (!is.null(object@dataPO)) {
  
    cat('Summary of the presense only datasets:')
    cat('\n\n')
    
    for (dataset in 1:length(object@dataPO)) {
      cat('Summary of:', names(object@dataPO)[[dataset]])
      cat('\n')
      print(summary(object@dataPO[[dataset]]))
      cat('\n')
              
      }
              
    cat('\n')  
              
    }
            
  if (!is.null(object@dataPA)) {
    
    cat('Summary of the presence absence datasets:')
    cat('\n\n')
            
    for (dataset in 1:length(object@dataPA)) {
      cat('Summary of:', names(object@dataPA)[[dataset]])
      cat('\n')        
      print(summary(object@dataPA[[dataset]]))
      cat('\n')
              
     }
            
  cat('\n')
            
    }  
            
    })
