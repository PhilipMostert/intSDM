## *intSDM* version 2.0

------------------------------------------------------------------------

### Package updates

-   Major changes to the primary functions of the package: both `structuredData` and `species_model` have now been depreciated.

    -   The package has now shifted to two new functions: `startWorkflow` and `sdmWorkflow` to help structure the workflow.

-   The package has depreciated any *sp* and *raster* support, and now uses *sf* and *terra* for the spatial data.

-   Cross-validation methods have been added to the workflow as a possible output.

-   The workflow now provides more meaningful messages throughout the whole process.

-   The workflow gives the user the option to save the model outputs in *R* formats.

-   The function now has the ability to provide metadata on the occurrence data obtained.
