source('R/functions-analyses.R')

# load volume data and air saturation data
vol    <-  readFile('data/vol_2015_12_17.csv')
ast    <-  readFile('data/airSat_2015_12_17.csv')
spid   <-  readFile('data/speciesID.csv')
o2tab  <-  data.frame()
for(i in 1:ncol(vol)) {
	dat  <-  data.frame(species=spid$species[spid$column==i], dry_mass_g=spid$dry_mass_g[spid$column==i], status=tolower(spid$status[spid$column==i]), shape=tolower(spid$shape[spid$column==i]), column=i, o2volume=vol[!is.na(vol[,i]),i], o2sat=ast[!is.na(ast[,i]),i], 	stringsAsFactors=FALSE)
	dat$boundedO2Vol               <-  dat$o2volume / max(dat$o2volume)
	dat$massSpecificO2Vol          <-  dat$o2volume / dat$dry_mass_g
	dat$boundedMassSpecificO2Vol   <-  dat$massSpecificO2Vol / max(dat$massSpecificO2Vol)
 	o2tab    <-  rbind(o2tab, dat)
}; rm(i, dat, vol, ast, spid)

o2tab$sppNum  <-  as.numeric(as.factor(tolower(gsub(' ', '', o2tab$species))))
