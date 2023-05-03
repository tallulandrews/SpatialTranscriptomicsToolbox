create_coords_grid <- function(obj, nrow=50, ncol=50) {
	tissue <- rep(1, nrow*ncol);
	get_coords <- function(x, n) {seq(from=min(x), to=max(x), length=n)}	
	row_coords <- get_coords(obj@images[[1]]@coordinates[,2], n=nrow)
	col_coords <- get_coords(obj@images[[1]]@coordinates[,3], n=ncol)
	row_coords2 <- get_coords(obj@images[[1]]@coordinates[,4], n=nrow)
	col_coords2 <- get_coords(obj@images[[1]]@coordinates[,5], n=ncol)
	df <- cbind(tissue, rep(row_coords, times=ncol), rep(col_coords, each=nrow),
					rep(row_coords2, times=ncol), rep(col_coords2, each=nrow))
	rownames(df) <- paste("Spot", 1:nrow(df), sep="_")
	return(df)
}

convert_hex_coord_to_cartesian <- function(hex_coords) {
	convert_mat <- rbind(c(1/2, 0, 1), c(sqrt(3)/2, sqrt(3), 0))
	cart <- convert_mat %in% hex_coords
	return(cart)
}

create_coords_hexes_wrong <- function(obj, nrow=50, ncol=50) {
	get_coords <- function(x, n) {seq(from=min(x), to=max(x), length=n/2)}
	col_coords_odd <- get_coords(obj@images[[1]]@coordinates[,3], n=ncol)
	row_coords_odd <- get_coords(obj@images[[1]]@coordinates[,2], n=nrow)
	row_coords_even <- row_coords_odd[-1] - diff(row_coords_odd)/2
	col_coords_even <- col_coords_odd[-1] - diff(col_coords_odd)/2
	
	all_coords <- cbind(
				c(rep(row_coords_odd, times=length(col_coords_odd)), rep(row_coords_even, times=length(col_coords_even))),
				c(rep(col_coords_odd, each=length(row_coords_odd)), rep(col_coords_even, each=length(row_coords_even)))
			)
}

create_coords_VISIUM_hexes <- function(nrow=77, ncol=127, pixels_x=c(1980,24014), pixels_y=c(2780, 23665)) {
	row_coords1 <- seq(1, nrow, by =2)
	row_coords2 <- seq(2, nrow, by =2)
	col_coords1 <- seq(1, ncol, by =2)
	col_coords2 <- seq(2, ncol, by =2)
	all_coords <- cbind(
				c(rep(row_coords1, times=length(col_coords1)), rep(row_coords2, times=length(col_coords2))),
				c(rep(col_coords1, each=length(row_coords1)), rep(col_coords2, each=length(row_coords2)))
			)
	rownames(all_coords) <- paste("Spot", 1:nrow(all_coords), sep="")
	pixels_x <- round(seq(pixels_x[1], pixels_x[2], length = max(c(row_coords1, row_coords2))))
	pixels_y <- round(seq(pixels_y[1], pixels_y[2], length = max(c(col_coords1, col_coords2))))
	
	all_coords_pixels <- cbind(pixels_x[all_coords[,1]], pixels_y[all_coords[,2]])

	return(cbind(rep(1, nrow(all_coords)), all_coords, all_coords_pixels))
}


get_stripes <- function(coords, nstripes=3, vertical=TRUE, plot=FALSE) {
	if (vertical) {
		out <- split(rownames(coords), cut(coords[,2], nstripes));
	} else {
		out <- split(rownames(coords), cut(coords[,3], nstripes));
	}
	names(out) <- paste("stripe", if (vertical){"V"}else{"H"}, 1:nstripes, sep="_");
	coords <- data.frame(coords)
	for (i in 1:length(out)) {
			coords$stripe[ rownames(coords) %in% out[[i]] ] <- names(out)[i]
	}


	if (plot) {
		cols <- RColorBrewer::brewer.pal(nstripes, "Set3")
		plot(coords[,2], coords[,3], pch=16, col=cols[factor(coords$stripe)], xlab="", ylab="")
	}
	return(coords)
}

get_concentric_circles <- function(coords, ncircles=3, plot=FALSE) {
	max_diameter <- min(max(coords[,2]), max(coords[,3]))
	center <- apply(coords[,2:3], 2, median)
	dist_thresholds <- seq(from=1, to=max_diameter/2, length=ncircles+1)
	spot_dists <- sqrt((coords[,2]-center[1])^2 + (coords[,3]-center[2])^2)
	out <- split(rownames(coords), cut(spot_dists, breaks=dist_thresholds))	
	names(out) <- paste("circle", 1:ncircles, sep="_");
	for (i in 1:length(out)) {
			coords$circle[ rownames(coords) %in% out[[i]] ] <- names(out)[i]
	}
	if (plot) {
		cols <- RColorBrewer::brewer.pal(ncircles+1, "Set3")
		plot(coords[,2], coords[,3], pch=16, col=cols[factor(coords$circle)], xlab="", ylab="")
	}
	return(coords)
}

get_random_points <- function(coords, nspots=0.1*nrow(coords), nsets=1, plot=FALSE) {
	coords$points <- rep(NA, nrow(coords))
	for (i in 1:nsets) {
		spots <- sample(1:nrow(coords), nspots)
		coords$points[spots] <- rep(paste("points", i, sep=""), length(spots))
	}
	if (plot) {
		cols <- RColorBrewer::brewer.pal(nsets, "Set3")
		plot(coords[,2], coords[,3], pch=16, col=cols[factor(coords$points)], xlab="", ylab="")
	}
	return(coords)
}


add_layer <- function(counts, to_sample, spot_ids, rescale=NULL) {
	to_add <- sample(1:ncol(to_sample), length(spot_ids), replace=TRUE)
	to_add <- counts[,to_add]
	if (!is.null(rescale)) {
		to_add <- t( t(to_add)/colSums(to_add) * colSums(counts[,spot_ids])*rescale )
	}
	counts[,spot_ids] <- counts[,spot_ids]+to_add
	return(counts)
}

sc <- readRDS("C:/Users/tandrews/Documents/UHNSonya/snRNAseq/SingleNuc_vs_SingleCell/DE_DataShare/Integrated_with_Subannotation.rds")
spatial <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Spatial/PSC011_4_C1_forLattice.rds")

		
