makeTable1  <-  function(dest, o2tab) {
	table1  <-  ddply(o2tab, .(species), function(x) {
					data.frame('Growth shape'   =  unique(x$shape),
						       'Status'         =  unique(x$status),
						       'DBM'            =  round(mean(tapply(x$dry_mass_g, x$column, unique)), 2),
						       'DBM_sd'         =  round(sd(tapply(x$dry_mass_g, x$column, unique)), 2),
						       'n'              =  length(unique(x$column)), stringsAsFactors = FALSE)
				})
	write.csv(table1, dest, row.names = FALSE)
}

makeTable2  <-  function(dest, fieldFlow) {
	table2  <-  ddply(fieldFlow, .(LocationNum), function(x) {
					data.frame(Site = unique(x$Location),
						       Mean = round(mean(x$X.AS, na.rm = TRUE), 2),
						       SD   = round(sd(x$X.AS, na.rm = TRUE), 2),
						       Min  = round(min(x$X.AS, na.rm = TRUE), 2),
						       Max  = round(max(x$X.AS, na.rm = TRUE), 2), stringsAsFactors = FALSE)
				})
	write.csv(table2, dest, row.names = FALSE)
}