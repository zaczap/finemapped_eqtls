core.populations = c('CEU','CHB','GIH','JPT','LWK','YRI')
core.colors = c('#0000FF','#ADCD00','#9400d3','#008b00','#CC9933','#FFB933') 

#' Cross-population eQTL and fine-mapping API
#' 
#' @export
#' @examples
#' eqtls = fetch_eqtl_data('ORMDL3')
#' eqtls = fetch_eqtl_data('BRCA1')
fetch_eqtl_data = function(gene) {
	require(data.table)
	multipopulation = fread(paste0('http://52.37.206.146:7999/download?switch=multi&gene=',gene))
	finemapped = fread(paste0('http://52.37.206.146:7999/download?switch=fine&gene=',gene))
	list(multipopulation=multipopulation, finemapped=finemapped)
}

#' Plot eQTL data for each population
#' 
#' @export
#' @examples
#' eqtls = fetch_eqtl_data('BRCA1')
#' plot_multipopulation_eqtls(eqtls, populations='CEU')
#' plot_multipopulation_eqtls(eqtls, populations=c('CEU','JPT'))
#' plot_multipopulation_eqtls(eqtls, start = 41000000, end = 41500000)
plot_multipopulation_eqtls = function(obj, populations=c('all'), start=-1, end=-1) {
	require(ggplot2)
	if('all' %in% populations) 
		include_populations = core.populations
	else 
		include_populations = core.populations[match(populations, core.populations)]
	if(length(include_populations) < 1) 
		stop(paste0('There are no valid populations in ', populations,'.'))
	include_populations = sort(include_populations)
	color_scale = scale_color_manual(name='Population', limits=include_populations, values=core.colors[match(include_populations, core.populations)])

	plot_df = subset(obj$multipopulation, population %in% include_populations)

	if(start > 0 && end > 0)
		plot_df = subset(plot_df, pos >= start & pos <= end)

  	ggplot(plot_df, aes(x=pos,y=pvalue,color=population)) + geom_point() + color_scale + 
  		ylab('-log10 p-value') + xlab(paste('Position on chromosome', plot_df$chrom[1])) + 
  		theme_bw() + theme(legend.key = element_blank()) + 
		guides(color = guide_legend(override.aes = list(size=3)))
}

#' Plot fine-mapping eQTL data for a specific discovery population
#' 
#' @export
#' @examples
#' eqtls = fetch_eqtl_data('ORMDL3')
#' plot_finemapped_eqtls(eqtls, 'JPT')
#' plot_finemapped_eqtls(eqtls, 'JPT', start = 38000000, 38050000)
plot_finemapped_eqtls = function(obj, discovery='CEU', start=-1, end=-1) {
	require(ggplot2)

	if(!(discovery %in% core.populations))
		stop(paste(discovery, "is not a valid population."))

	plot_df = subset(obj$finemapped, discovery == population)
	if(start > 0 && end > 0)
		plot_df = subset(plot_df, pos >= start & pos <= end)

	if(nrow(plot_df) == 0)
		warning(paste0('There is no finemapping eQTL data for population ',discovery,' in this region.'))

  	ggplot(plot_df, aes(x=pos,y=pvalue)) + geom_point(color='black') + ylab('-log10 p-value') + xlab(paste('Position on chromosome', plot_df$chrom[1])) + theme_bw() + theme(legend.position='none')
}

#' Plot both multipopulation eQTL data and fine-mapping eQTL data for a specific population.
#' 
#' @export
#' @examples
#' eqtls = fetch_eqtl_data('SP1')
#' plot_eqtl_locus(eqtls)
#'
#' eqtls = fetch_eqtl_data('ORMDL3')
#' plot_eqtl_locus(eqtls, populations=c('CHB','JPT'), discovery='JPT')
plot_eqtl_locus = function(obj, populations=c('all'), discovery='CEU', start=-1, end=-1) {
	require(ggplot2)
	require(cowplot)

	if(!(discovery %in% core.populations))
		stop(paste(discovery, "is not a valid population."))
	fig1 = plot_multipopulation_eqtls(obj, as.vector(populations), start, end) + labs(x='')
	xmin = layer_scales(fig1)$x$range$range[1]
  	xmax = layer_scales(fig1)$x$range$range[2]

	ticks = ggplot_build(fig1)$panel$ranges[[1]]$x.major_source
	minor_ticks = ggplot_build(fig1)$panel$ranges[[1]]$x.minor_source

	plot_df = subset(obj$finemapped, discovery == population)
	plot_df = subset(plot_df, pos >= xmin & pos <= xmax)

	if(nrow(plot_df) == 0) {
	    plot_df = rbind(plot_df, list('', xmin, '', 0, 'No data'))
	    plot_df = rbind(plot_df, list('', xmax, '', 4, 'No data'))
	    populations = c('No data')
	    colors = c('white')
	    fig2 = ggplot(plot_df, aes(x=pos,y=pvalue,color=population)) + geom_point() + scale_color_manual(name='', breaks=populations,values=colors) + ylab('-log10 pvalue') + xlab(paste('Position on chromosome',obj$multipopulation$chrom[1])) + theme_bw() + scale_x_continuous(limits=c(xmin,xmax),breaks=c(ticks),minor_breaks=minor_ticks) + guides(color = guide_legend(override.aes = list(size=3))) + theme(legend.key = element_blank())
	} else {
		fig2 = ggplot(plot_df, aes(x=pos,y=pvalue, color=population)) + geom_point() + scale_color_manual(name='', breaks=discovery,values='black') + ylab('-log10 pvalue') + xlab(paste('Position on chromosome',obj$multipopulation$chrom[1])) + theme_bw() + scale_x_continuous(limits=c(xmin,xmax),breaks=c(ticks),minor_breaks=minor_ticks) + guides(color = guide_legend(override.aes = list(size=3))) + theme(legend.key = element_blank())
	}
	plot_grid(fig1, fig2, labels=c('A','B'), nrow=2, align='v')
}
