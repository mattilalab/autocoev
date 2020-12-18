######################################################################################
#How to use

##########################################################
##Pass input and app location as arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
appLocation<-args[2]

shiny::runApp(appLocation, launch.browser = TRUE)
