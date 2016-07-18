###################
# WRAPPER FUNCTIONS
###################
readFile  <-  function(filepath){
    read.csv(filepath, header=TRUE, stringsAsFactors=FALSE, na.strings=c('','NA'), strip.white=TRUE)  
}

cutOffSolution97  <-  function(denPar) {
    # calculated based on the formula dV/dx / dV/d(x-1) = 0.97
    # where dV/dx = (a*b)/(b+x)^2, and a and b are parameters from Michaelis-Menten formula, so
    # (b+(x-1))^2 / (b+x)^2 = 0.97
    # parameters above are not to be confounded with a-d below which belong to Bhaskara formula
    a  <-  0.03
    b  <-  (0.06*denPar - 2)
    c  <-  (0.03*denPar^2 - 2*denPar + 1)
    d  <-  b^2 - 4*a*c
    (-b+sqrt(d))/(2*a) # only the positive solution from Bhaskara is wanted
}

#############################
# DEALING WITH BIB REFERENCES
#############################
getIndividualBibs  <-  function(path='~/bibtex_library/') {
    listRefs     <-  dir(path)[grep('.bib', dir(path))]
    listRefs[listRefs != 'library.bib']
}

readBibRefs  <-  function(individualBibs, path='~/bibtex_library/') {
    lapply(individualBibs, function(x, path)read.bib(file.path(path, x)), path=path)
}

listBibs  <-  function() {
    listRefs        <-  getIndividualBibs()
    theRefs         <-  readBibRefs(listRefs)
    names(theRefs)  <-  gsub('.bib', '', listRefs)
    theRefs
}

exporBibs  <-  function(theRefs=listBibs(), libraryPath='../library.bib', erase=TRUE) {
    if(erase)
        system(paste('rm -r', libraryPath))
    if(!file.exists(libraryPath))
        l_ply(theRefs, function(x, libraryPath)write.bibtex(x, file=libraryPath, append=TRUE), libraryPath=libraryPath)
}

myCite  <-  function(citationsVec) {
    paste0('[@', paste0(citationsVec, collapse=';@'),']')
}
