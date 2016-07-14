source('R/functions-analyses.R')

# load volume data and air saturation data
vol    <-  readFile('data/vol_2015_12_17.csv')
ast    <-  readFile('data/airSat_2015_12_17.csv')
spid   <-  readFile('data/speciesID.csv')
o2tab  <-  data.frame()
for(i in 1:ncol(vol)) {
	column   <-  as.numeric(sub('VO.', '', names(vol)[i]))
	spMatch  <-  spid$column == column
	dat      <-  data.frame(species=spid$species[spMatch], dry_mass_g=spid$dry_mass_g[spMatch], status=tolower(spid$status[spMatch]), shape=tolower(spid$shape[spMatch]), column=column, o2volume=vol[!is.na(vol[,i]),i], o2sat=ast[!is.na(ast[,i]),i], 	stringsAsFactors=FALSE)
	dat$boundedO2Vol               <-  dat$o2volume / max(dat$o2volume)
 	o2tab    <-  rbind(o2tab, dat)
}; rm(i, dat, vol, ast, spid)

o2tab$sppNum  <-  as.numeric(as.factor(tolower(gsub(' ', '', o2tab$species))))
