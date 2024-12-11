## This function works with both 2D and 3D; no need to specify which it is
# X is a configuration of one specimen in the format of k landmarks by m dimensions
# "ventral" and "dorsal" are vectors indicating the names (NOT indices)
# of the ventral and dorsal landmarks, respectively.
# comV and comD are vectors indicating the common landmarks' names (NOT indices) for the ventral and dorsal sides respectively (in the same order);
# the function attaches the ventral side to the dorsal side
# (or any other two aspects of the same structure) based on the common landmarks.
# It finds the rotation matrix that minimizes the sum of squared deviations
# between the ventral and the dorsal
# (Rohlf 1990. Proceedings of the Michigan Morphometrics Workshop),
# and calculates the errors as the vector length between the ventral and the dorsal
# equivalent landmarks after fitting.
# If either the dorsal or the ventral copy of a common landmark is missing,
# that landmark is ignored for the unification.
# Any other missing data does not affect the procedure at all and can be safely included.
# If average=TRUE then the dorsal and ventral copies of the common landmarks are averaged
# and listed only for the dorsal, so the output will have fewer rows than the input.

# NOTE: the resulting configuration will be re-organized to have all the dorsal
# landmarks first, then all the ventral ones.
# It is therefore HIGHLY RECOMMENDED that landmarks will be designated by names in X
# (i.e., rownames(X) is a charcter vector rather than NULL)

# The output includes the unified data ($unified) and the unification errors ($errors)

# an example for an input file and a protocol that implements this function can be downloaded from
# http://home.uchicago.edu/~annat/
# Please email me with comments and questions annat22@gmail.com
# last updated Feb 28 2011
#' @export
unifyViews <- function(X, V1, V2, comLands, average=FALSE) {
		V <- X[comLands,]
		D <- X[comLands,]
		Xv <- X[V1,]
		Xd <- X[V2,]
		V[which(is.na(D))] <- NA
		D[which(is.na(V))] <- NA # making sure both ventral and dorsal of the same LM are NA's whenever one of them is
		mV <- matrix(apply(V, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xv), nc=ncol(Xv))
		Xvc <- Xv-mV # translating all the ventral LM's based on the centroid of the common ones
		Vc <- scale(V, scale=F)
		mD <- matrix(apply(D, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xd), nc=ncol(Xd))
		Xdc <- Xd-mD # translating all the dorsal LM's based on the centroid of the common ones
		Dc <- scale(D, scale=F)
		M <- t(na.omit(Dc)) %*% na.omit(Vc) # arbitrarily choosing the dorsal as the reference
		SVD <- svd(M)
		L <- diag(SVD$d)
		S <- ifelse(L<0, -1, L)
		S <- ifelse(L>0, 1, L)
		RM <- SVD$v %*% S %*% t(SVD$u) # the rotation matrix
		Xvr <- Xvc %*% RM # rotate all the translated ventral LM's
		dv <- rbind(Xdc, Xvr)
		dv[comV[which(is.na(dv[comV,1]))],] <- dv[comD[which(is.na(dv[comV,1]))],]
		dv[comD[which(is.na(dv[comD,1]))],] <- dv[comV[which(is.na(dv[comD,1]))],]
		if (average==TRUE) {
			dv[comD,] <- (dv[comD,]+dv[comV,])/2
			dv <- dv[-match(comV, rownames(dv)),]}
		list(unified=dv, errors=sqrt(rowSums((Xvr[comV,]-Dc)^2)))
		}
