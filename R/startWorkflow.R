#' @description
#' Function to initialize the reproducible workflow using integrated species distribution models.
#' @param Countries A vector of country names to complete the analysis over. If missing, a boundary object may be added to the model using \code{.$addArea}.
#' @param Species A vector of Species names (scientific) to include in the analysis.
#' @param Save Logical argument indicating if the model objects and outputs should be saved. Defaults to \code{TRUE}.
#' @param saveOptions A list containing two items: \code{projectDirectory} indicating where the objects should be saved, and \code{projectName} which indicates the name for the folder in the relevant directory.
#' @param Quiet Logical argument indicating if the workflow should provide the user messages during the setup and estimation process. Defaults to \code{TRUE}.
#'
#' @returns An R6 object of class \code{species_model}. This object contains a collection of slot functions to assist the user in customizing their workflow.
#'
#'

startWorkflow <- function(Countries, Species,
                          Projection,
                          Save = TRUE,
                          saveOptions = list(
                          projectDirectory = getwd(),
                          projectName = NULL),
                          Quiet = TRUE) {

  if (Save) {

  if (!is.character(saveOptions$projectName)) stop('Please provide projectName in the saveOptions list.')

  }

  if (missing(Species)) stop('At least one species name needs to be provided.')

  if(!inherits(Projection, 'CRS')) stop('Projection needs to be a CRS object.')

  if (Quiet) {

  cat('Initializing workflow for integrated species distribution model:\n\n')

  cat('Studied species:', paste(Species, collapse = ', '), '\n')

  if (!missing(Countries)) {

  if (length(Countries) == 1) cat('Studied country:', Countries, '\n\n')
  else cat('Studied countries:', paste(Countries, collapse = ', '), '\n\n')


  } else message('Countries not specified. You can add a sampling region with the `.$addArea` function.')


  cat('This function creates an object with "slot" functions to help you customize the ISDM for your workflow. These may be accessed by using the "$" after the name of the object created with this function.\n\n')
  cat('The following slot functions are availble with this object: \n')

  descriptionSlots <- data.frame(Name = c('---------------','plot',
                                          'addStructured',
                                          'addMesh',
                                          'addGBIF',
                                          'addArea',
                                          'addCovariates',
                                          'crossValidation',
                                          'modelOptions',
                                          'specifySpatial',
                                          'biasFields',
                                          'workflowOutput'),
                                 Desription = c('---------------','-> Plot data',
                                                 '-> Add structured data',
                                                 '-> Create an inla.mesh object',
                                                 '-> Add data from GBIF',
                                                 '-> Specify sampling region',
                                                 '-> Add spatial covariates',
                                                 '-> Specify CV method',
                                                 '-> Add INLA model options',
                                                 '-> Specify spatial effects',
                                                 '-> Add bias field',
                                                 '-> Output of workflow'))
print.data.frame(descriptionSlots, right = FALSE, row.names = FALSE)

cat('\nThe workflow may then be estimated using the "sdmWorkflow" function.\n\n')

if (Save) {

dir.create(path = paste0(saveOptions$projectDirectory, '/', saveOptions$projectName))

cat('Directory for model outputs is:\n', paste0(saveOptions$projectDirectory, '/', saveOptions$projectName))

}

}

  modelSetup <- species_model$new(Countries = Countries,
                                  Species = Species,
                                  Projection,
                                  Quiet = Quiet,
                                  nameProject = saveOptions$projectName,
                                  Directory = saveOptions$projectDirectory)

}
