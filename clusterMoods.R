###################
## clusterMoods.R
##
## clusters the mood data
##
##################

require('rdist') ## for calculating distances
require('data.table') ## for the rbindlist function
require('Matrix') ## for converting the transition matrix to an edgelist

# radialToEuclid <- function(point){
# 	alpha <- atan2(point$activation, point$affect)
# 	x <- sin(alpha)
# 	y <- cos(alpha)
	
## helper function which binds the cluster ID to the 
## (global) moods data
mkCluster <- function(index, membershipList){
	members <- membershipList[[index]]
	cluster <-  rep(index,length(members))
	cbind(moods[ members, ], cluster)
}


## step (climbers, reference)
## perform a hillcimbing step
## return the new positions of all points
step <- function(climbers, reference, bandwidth=10){
	dmat <- cdist(climbers, reference) / bandwidth
	dmat[is.nan(dmat)] <- 0 ## for safety; 
	kern <- exp(-0.5 * dmat * dmat) 
	if(sum(is.na(kern)) > 0){ cat("--- error in hillclimbing -----\n") }
	allinfl <- rowSums(kern)
	xstep <- rowSums(kern %*% climbers[,1] ) / allinfl
	ystep <- rowSums(kern %*% climbers[,2] ) / allinfl
	newpos <- cbind(xstep, ystep)
	stepsize <- diag(cdist(climbers,newpos)) ## TODO!! not very efficient
	stepsize[is.nan(stepsize)] <- 0
	deltaF <- sum(allinfl)
	return (list(newpos = newpos, stepsize = stepsize,
							 deltaF = deltaF))
}


findClusters <- function(moods){
	## extract the point locations
	points <- cbind(moods$activation, moods$affect)
	## TODO! convert points to euclidean space!

	## set parameters
	epsilon <- 1e-3
	score <- c(0)
	deltaScore <- 100
	ssize <- matrix(0,ncol=2,nrow=dim(points)[1])
	numNotUnique <- 1

	## initialize hillclimbing
	hillclimbOut <- step(points,points)
	# plot(hillclimbOut$newpos)
	while(numNotUnique > 0){
		## hillclimb to classify
		while(deltaScore > epsilon){
			hillclimbOut <- step(hillclimbOut$newpos, hillclimbOut$newpos)
			# points(hillclimbOut$newpos,pch=12,col=colors()[10*length(score)])
			## calculate stopping parameter
			deltaScore <- abs(hillclimbOut$deltaF - score[length(score)])/hillclimbOut$deltaF
			ssize[,1] <- ssize[,2]
			ssize[,2] <- hillclimbOut$stepsize
			## keep a running score (just for fun, since a scalar would
			## work as well to compute the delta
			## score should increase-- as points get closer to each other,
			## the density increases.
			score <- c(score, hillclimbOut$deltaF)
			if(length(score)>50) break
			# cat(hillclimbOut$deltaF, '\n')
		}
		## check if this results in unique classification	
		dmat <- pdist(hillclimbOut$newpos)
		dmat[is.nan(dmat)] <- 0
		balls <- rowSums(ssize)
		isClassified <- rep(FALSE, dim(points)[1])
		## create boolean matrix, with rows the points,
		## and "true" indicating that a member of the row
		## is in the same cluster as the point
		clustMat <- apply(cbind(balls,dmat),1,
											function(row){row[-1] <= row[1]})
		## thus unique(clustMat) gives the clusters, and the items in each
		## row of clustMat are the cluster members.
		numNotUnique <- sum(colSums(unique(clustMat)) != 1)
		## we'll reduce episilon here, though this is only needed if 
		## the outer while loop is triggered
		epsilon <- min(epsilon * 0.2, deltaScore * 0.8)
	}
	## assign cluster IDs to each entry
	membership <- apply(unique(clustMat),1,function(row){which(row)})
	clustered <- rbindlist(lapply(seq_along(membership), 
		function(index){
			members <- membership[[index]]
			clusterID <-  rep(index,length(members))
			cbind(moods[ members, ], clusterID)
		}
		))
	return(as.data.frame(clustered))
}


## sorts the nodes based on timestamp.
## uses this to create a sparse transition 
## matrix, which is returned as an edge list.
findEdges <- function(clusters){
	## order nodes by time
	clustNames <- as.factor(clusters[order(clusters$timestamp),]$cluster)

	## make transition matrix
	nclust <- length(unique(clustNames))
	tmat <- Matrix(0, nrow=nclust, ncol=nclust)
	## loopy loop
	for(i in 2:length(clustNames)){
		prior <- as.numeric(clustNames[i-1])
		current <- as.numeric(clustNames[i])
		tmat[ prior, current ] <- tmat[ prior, current ] + 1
	}	
	edges <- summary(tmat)
	names(edges) <- c("source", "target", "value")
	return(edges)
}

#* @filter cors
 # cors <- function(res) {
 #     res$setHeader("Access-Control-Allow-Origin", "*")
 #     plumber::forward()
 # }


#* @post /moodgraph
function(moods){
	## error handling 
	# if(!moods){ return(list(error="no mood data")) }
	if(dim(moods)[1]< 20){ 
		cat('error\n')
		msg <- paste("insufficient mood history, only ",
								 dim(moods)[1], "entries")
		return(list(error=msg))
	}

	nodes <- findClusters(moods)
	links <- findEdges(nodes)
	return(list(nodes=nodes,links=links))
}

#* @get /echo
#* @post /echo
function(msg=""){
  list(msg = paste0("The message is: '", msg, "'"))
}

#* @get /
function() { Sys.Date() }
