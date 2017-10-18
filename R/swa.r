#' Subsampling Winner Algorithm for Variable Selection in the Classification Setting
#' 
#' This function conducts variable selection in the classification setting. The algorithm repeatedly subsamples variables and runs linear discriminant analysis (LDA) via linear regression based on the subsampled variables. Variables are scored based on the AUC and the t-statistics in linear regression. Variables then enter a competition and the semi-finalists will be used in a final round of LDA (or linear regression.) The semi-finalists and the final linear model are returned.
#'
#' @param X The p by n design matrix where p is the dimension and n is the sample size.
#' @param y The binary class labels (0 or 1.)
#' @param S Size(s) of subsamplings. This could be a set of candidate sizes, or a single size. In the former case, all sizes will be tried. A vector of default values is c(3,5,7,10,15,20,30,40).
#' @param q Number of semi-finalists. The default is 20.
#' @param m Number of subsamplings. The default is 500.  One may increase the value of m to see if the outcomes have been stabilized.
#' @param MPplot Boolean input showing whether the multi-panel diagnostic plot is drawn. The default is FALSE.
#' @param screening Boolean input showing whether the variables are pre-screened based on marginal variance. The default is TRUE.
#'
#' @return A list with the indices of the semi-finalists, summary of the final regression model, and summary of the stepwise selection of the final model.
#'
#' @keywords classification, variable selection, subsampling.
#'
#' @export
#' @import ROCR ggplot2 reshape
#' 
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(20*100), ncol=20)
#' b = c(0.5,1,1.5,2,3); 
#' y <- x[,1:5]%*%b + rnorm(100)
#' y <- as.numeric(y<0)
#' X <- t(x)
#' \donttest{swa(X,y,c(3,5,10,15,20),5,500,MPplot = TRUE) ;}
#' ## The MP plots show that the upper arm sets, i.e., the points above the elbow point of each 
#' ## panel plot, started to be stabilized at s = 10. Thus, we fix s = 10 below.
#' sum_swa <- swa(X,y,10,5,500)
#' sum_swa 
#' \donttest{
#' ## The final model perfectly recovers the five true variables. The stepwise selection removes
#' ## X1, which has the weakest signal (hence is not expected to be recovered.) We next check the
#' ## results with a larger m = 1000.
#' sum_swa <- swa(X,y,10,5,1000)   
#' sum_swa
#' # It now captured X2,X3,X4,X5, and also a noise variable X6 that is not as significant as the 
#' ## former 4 covariates X2-X4. We next look at the estimates of all 20 coefficients including 
#' ## those not selected.
#' coef_swa = numeric(20)
#' column_index = as.numeric(substring(rownames(sum_swa$step_sum$coef[-1,]),2,5))
#' coef_swa[column_index] <- sum_swa$step_sum$coef[-1,1]
#' 
#' #### Compare the coefficients from the SWA on all covariates with the coefficients from the 
#' ## oracle LDA. We first calculate the coefficients from the oracle LDA.
#' n = 100
#' n0 = sum(y==0)
#' n1 = sum(y==1)
#' yy=y
#' yy[y==0] <- -n/n0 ;
#' yy[y==1] <- n/n1 ;
#' X = X[1:5,]
#' X = X - apply(X,1,median);
#' X = X/apply(X,1,mad) ;
#' designX = t(X) ;
#' colnames(designX) <- paste('X',1:5,sep='') ;
#' dataXY = data.frame(cbind(yy,designX)) ; dim(dataXY) ;
#' LDA <- lm(yy ~ ., data = dataXY)
#' sum_LDA <- summary(LDA)
#' coef_LDA = numeric(20)
#' coef_LDA[as.numeric(substring(rownames(sum_LDA$coef[-1,]),2,5))] <- sum_LDA$coef[-1,1]
#' ####  Compare the coefficients using correlation.
#' cor(coef_swa,coef_LDA); rbind(coef_swa, coef_LDA); 
#' ####  There is an excellent correlation between the coefficients recovered from the SWA on 
#' ## the full data and the coefficients from the benchmark, ie., the oracle LDA.
#' 
#' ####  Next we try on a set of smaller coefficients that should be harder to recover from 
#' ## the noisy data
#' x <- matrix(rnorm(20*100), ncol=20)
#' b = c(0.3,0.6,0.9,1.2,1.5);
#' y <- x[,1:5]%*%b + rnorm(100)
#' y <- as.numeric(y<0)
#' X <- t(x)
#' swa(X,y,10,5,500) ;
#' ####  SWA does a perfect job in recovering the true variables.}

swa <- function(X, y, S = c(3,5,7,10,15,20,30,40), q = 20, m = 500, MPplot = FALSE, screening = TRUE){

	# require(ROCR)
	
	# if(MPplot){
		# require(reshape)
		# require(ggplot2)
	# }
	
	d = dim(X)[1] ;
	n = dim(X)[2] ;
	# if(length(S)==1){q = S}
	
	Tmat = array(NA,dim=c(d,m,length(S))) ;
	AUCv = matrix(NA,m,length(S)) ;
	scores = matrix(NA,d,length(S)) ;
	finalindex = matrix(NA,d,length(S)) ;
		
	n0 = sum(y==0)
	n1 = sum(y==1)
	keep = 1:d
	
	if(screening){
		fullid = 1:d ;
		ttestXQ <- function(X,y){return(stats::t.test(X[y==0],X[y==1])$statistic)}
		keep = fullid[abs(apply(X,1,ttestXQ,y=y)) > 1e-3] ;
	}
	
	yy=y
	yy[y==0] <- -n/n0 ;
	yy[y==1] <- n/n1 ;
	
	X = X - apply(X,1,stats::median);
	X = X/apply(X,1,stats::mad) ;
	designX = t(X) ;
	
	if(is.null(rownames(X))){
		colnames(designX) <- paste('X',1:d,sep='') ;
	}else{
		colnames(designX) <- rownames(X) ;
	}
	
	for(is in 1:length(S)){
		s = S[is]
		counter = 0 ;
		repeat{
			if(counter == m){break ;}
			index = sample(keep,s) ;
			dsgnX = designX[,index] ;
			dataXY = data.frame(cbind(yy,dsgnX)) ;
			lmmodel <- stats::lm(yy ~ ., data = dataXY)
			counter = counter + 1 ;
			fits = as.numeric(lmmodel$fitted.values) ;
			
			pred <- ROCR::prediction(fits,y)
			AUCv[counter,is] = ROCR::performance(pred,'auc')@y.values[[1]] ;
			if(length(index)>0) {Tmat[index,counter,is] = summary(lmmodel)$coef[-1,"t value"]}
		}
		indexKept =	order(AUCv[,is],decreasing=TRUE)[!duplicated(sort(AUCv[,is],decreasing=TRUE))][1:q] ; # ignore the duplicated resamplings.
		reduceTmat = Tmat[,indexKept,is] ; # choose the best q models out of the m replications
		reduceAUC = AUCv[indexKept,is] ; # choose the best q models out of the m replications
		weightedTmat = abs(t(t(reduceTmat)*(reduceAUC^0.5))) ;
		score_temp = sort(rowMeans(weightedTmat,na.rm = TRUE),decreasing = TRUE)
		scores[1:length(score_temp),is] = score_temp ;
		finalindex_temp = order(rowMeans(weightedTmat,na.rm = TRUE),decreasing = TRUE) ; # choose the best q variables.
		finalindex[1:length(score_temp),is] = finalindex_temp[1:length(score_temp)] ;
	}
	
	if(MPplot){
		scores_df = reshape::melt(scores)
		scores_df = data.frame(scores_df)
		names(scores_df) = c('index','s','score')
		scores_df = scores_df[!is.na(scores_df$score),]
		scores_df$s = factor(S[scores_df$s])
		
		finalindex_df = reshape::melt(finalindex)
		finalindex_df = data.frame(finalindex_df)
		names(finalindex_df) = c('index','s','variable')
		finalindex_df = finalindex_df[!is.na(finalindex_df$variable),]
		finalindex_df$s = factor(S[finalindex_df$s])
		
		combine = merge(scores_df,finalindex_df)

		score<-variable<-NULL
		ggplot2::ggplot(subset(combine,index<=25), ggplot2::aes(x=index, y=score, color=s, label = variable))	+ ggplot2::geom_text() + ggplot2::facet_grid(s ~ .)
		
	}else if(length(S)==1){
		finalindex1 = finalindex[1:q,1]
		dsgnX = designX[,finalindex1] ;
		dataXY = data.frame(cbind(yy,dsgnX)) ;
		final.lmmodel <- stats::lm(yy ~ ., data = dataXY)
        final.sum = summary(final.lmmodel)
		step_sum = summary(stats::step(final.lmmodel,trace=0))
        return(list(index = finalindex1, summary = final.sum, step_sum = step_sum))
	}
}