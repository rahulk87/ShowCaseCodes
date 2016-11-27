make_scatter_plot <- function(
	x0,
	x1,
	gene="PARP1",
	filename="plot.pdf",
	x0name="library",
	x1name="sample",
	response_colname="pptm_psuedo",
	main=""
	){
	
	pdf(filename, width=5, height=5)

	plot(
		x0[,response_colname],
		x1[,response_colname],
		log="xy",
		xlab=x0name,
		ylab=x1name,
		main=main
		)

	rows_to_mark <- grep(gene, x0$shrna.id)

	for(row in rows_to_mark){
		points(
			x0[row,response_colname],
			x1[row,response_colname],
			pch=19,
			col="red"
			)
	}

	dev.off()
}
