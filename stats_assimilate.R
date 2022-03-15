#Assimilating mean of all indices over all dates
library(gplots)
library(lattice)
#collecting all stats file and reading them into a list
files <- list.files("./", recursive = T, pattern = "norm_Indices_stats.csv")
In_stats <- list()
for(i in seq_along(files))
{
In_stats[[i]] <- read.csv(files[i])
}
#calculating basic image statistics for all indexes through the acquisition dates
MN <- lapply(In_stats, function(x) x$Mean)
MD <- lapply(In_stats, function(x) x$MD)
SD <- lapply(In_stats, function(x) x$SD)
CV <- lapply(In_stats, function(x) x$CV)
SK <- lapply(In_stats, function(x) x$SK)
Mean_stats <- do.call(cbind, MN)
Median_stats <- do.call(cbind, MD)
Stdv_stats <- do.call(cbind, SD)
CVar_stats <- do.call(cbind, CV)
Skew_stats <- do.call(cbind, SK)
L <- list(Mean_stats,Median_stats,Stdv_stats,CVar_stats,Skew_stats)
names(L)<-c("Mean_stats","Median_stats","Stdv_stats","CVar_stats","Skew_stats")
L <- lapply(L, function(x){colnames(x) <- c("92DAS","105DAS","118DAS","133DAS")
 							x
				  }
		)
L <- lapply(L, function(x){rownames(x) <- In_names
							x
				  }
		)
lapply(1:length(L), function(i) write.csv(L[[i]],file = paste0(names(L[i]), ".csv"))

}