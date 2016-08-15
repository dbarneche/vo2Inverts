###################
# WRAPPER FUNCTIONS
###################
readFile  <-  function(filepath){
    read.csv(filepath, header=TRUE, stringsAsFactors=FALSE, na.strings=c('','NA'), strip.white=TRUE)  
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
