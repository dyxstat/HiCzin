#Normalize metagenomic Hi-C data and detect spurious contacts using zero-inflated Negative Binominal regression frameworks
#Auther and maintainer: Yuxuan Du <yuxuandu@usc.edu>
#HiCzin R script depends on 'glmmTMB' package and 'optparse' package

library("optparse")
library("glmmTMB")
options(scipen = 999)

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA,
              help = "input raw metagenomic Hi-C contacts"),
  make_option(c("-s", "--sample"), type = "character", default = NA,
              help = "input samples of the intra-species contacts"),
  make_option(c("-o", "--output"), type = "character", default = NA,
              help = "output normalized Hi-C contacts by HiCzin"),
  make_option(c("-d", "--discard"), type = "character", default = NA,
              help = "output normalized Hi-C contacts that are regarded as spurious contacts and discarded"), 
  make_option(c("-f", "--info"), type = "character", default = NA,
              help = "contig information file"),
  make_option(c("-c", "--coverage"), type = "character", default = NA,
              help = "contigs' converage file"),
  make_option(c("-t", "--threshold"), type = "character", default = 0.1,
              help = "preselected threshold for discarding spurious contacts"),
  make_option(c("-u", "--unlabeled"), action = "store_true", default = FALSE,
              help = "use ulabeled mode of HiCzin"),
  make_option(c("-n", "--nb"), action = "store_true", default = FALSE,
              help = "use negative binomial regression")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(is.na(opt$input) ) {
  stop("Missing input option of raw metagenomic Hi-C contacts.")
}

if(is.na(opt$output) ) {
  stop("Missing output option of normalized metagenomic Hi-C contacts.")
}

if(is.na(opt$info) ) {
  stop("Missing contig information file.")
}

if(is.na(opt$coverage) ) {
  stop("Missing contigs' coverage file.")
}

if(!opt$unlabeled & is.na(opt$sample)){
  stop("Missing sample data to fit the zero-inflated negative binomial model")
}



allcontact_file_name = opt$input
sample_file_name =opt$sample
output_file_name = opt$output
discarded_contact_name = opt$discard
contig_info_name = opt$info
coverage_name = opt$coverage
thres = as.numeric(opt$threshold)

if ( !file.exists(allcontact_file_name) ) {
  stop("Cannot open the input file of raw metagenomic Hi-C contacts..")
}
if ( !file.exists(sample_file_name) ) {
  stop("Cannot open sample data of the intra-species contacts.")
}
if ( !file.exists(contig_info_name) ) {
  stop("Cannot open the file of contig information")
}
if ( !file.exists(coverage_name) ) {
  stop("Cannot open the coverage file")
}

all_contacts = read.csv(allcontact_file_name , header = F , sep = ',' )
all_contacts = as.data.frame(all_contacts)
colnames(all_contacts) = c('index1' , 'index2' , 'contacts')

sample_data = read.csv(sample_file_name , header = F , sep = ',' )
sample_data = as.data.frame(sample_data)
colnames(sample_data) = c('index1' , 'index2' , 'contacts')

contig_info = read.csv(contig_info_name , header = F , sep = ',' )
contig_info = as.data.frame(contig_info)

if(ncol(contig_info) == 3){
  colnames(contig_info) = c('contig_name' , 'site' , 'length')
  contig_info[contig_info$site==0 , 2] = 1
  
  coverage = read.table(coverage_name , header = F , sep = '' , skip= 1)[ , c(1 , 6)]
  if(sum(contig_info[ , 1]!=coverage[ , 1])>0)
  {
    stop('Order mistakes exist in converage file')
  }
  contig_info$coverage = coverage[ , 2]
  
  
  sample_len = rep(0 , nrow(sample_data))
  sample_site = rep(0 , nrow(sample_data))
  sample_cov = rep(0 , nrow(sample_data))
  
  for(i in 1:nrow(sample_data))
  {
    sample_site[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 2]) * 
                           as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 2]))
    
    sample_len[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 3]) * 
                          as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 3]))
    
    sample_cov[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 4]) * 
                          as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 4]))
  }
  
  sampleCon = as.numeric(sample_data[ , 3])
  
  mean_site = mean(sample_site)
  sd_site = sd(sample_site)
  mean_len = mean(sample_len)
  sd_len = sd(sample_len)
  mean_cov = mean(sample_cov)
  sd_cov = sd(sample_cov)
  
  sample_site = (sample_site-mean_site)/sd_site
  sample_len = (sample_len-mean_len)/sd_len
  sample_cov = (sample_cov-mean_cov)/sd_cov
  
  data_sample = cbind(sample_site , sample_len , sample_cov , sampleCon)
  data_sample = as.data.frame(data_sample)
  colnames(data_sample) = c('sample_site' , 'sample_len' , 'sample_cov' , 'sampleCon')
  
  
  all_len = rep(0 , nrow(all_contacts))
  all_site = rep(0 , nrow(all_contacts))
  all_cov = rep(0 , nrow(all_contacts))
  
  for(i in 1:nrow(all_contacts))
  {
    all_site[i] = log(as.numeric(contig_info[as.numeric(all_contacts[i , 1]) , 2]) * 
                        as.numeric(contig_info[as.numeric(all_contacts[i , 2]) , 2]))
    
    all_len[i] = log(as.numeric(contig_info[as.numeric(all_contacts[i , 1]) , 3]) * 
                       as.numeric(contig_info[as.numeric(all_contacts[i , 2]) , 3]))
    
    all_cov[i] = log(as.numeric(contig_info[as.numeric(all_contacts[i , 1]) , 4]) * 
                       as.numeric(contig_info[as.numeric(all_contacts[i , 2]) , 4]))
  }
  
  allCon = as.numeric(all_contacts[ , 3])
  all_site = (all_site-mean_site)/sd_site
  all_len = (all_len-mean_len)/sd_len
  all_cov = (all_cov-mean_cov)/sd_cov
  
  
  tryCatch(
    {
      message(paste("normalizing",sep=" "))
      
      if(opt$unlabeled){
        fit1 = glmmTMB(sampleCon~sample_site+sample_len+sample_cov, data = data_sample,
                       ziformula=~1,family=nbinom2)
        
        
      }else if(opt$nb){
        fit1 = glmmTMB(sampleCon~sample_site+sample_len+sample_cov, data = data_sample,
                       ziformula=~0,family=nbinom2)
        
      }else{
        fit1 = glmmTMB(sampleCon~sample_site+sample_len+sample_cov, data = data_sample,
                       ziformula=~sample_site+sample_len+sample_cov , family=nbinom2)
      }
    },
    error = function(e){
      message(e)
      message(paste("\nskip",  sep=" "))
    },
    warning = function(w){
      message(w)
      message(paste("\nskip",  sep=" "))
    }
  )
  
  
  coeff = as.numeric(fit1$fit$par)
  res_sample = sampleCon/exp(coeff[1] + coeff[2]*sample_site + coeff[3]*sample_len+ coeff[4]*sample_cov)
  mu_sample = exp(coeff[1] + coeff[2]*sample_site + coeff[3]*sample_len+ coeff[4]*sample_cov)
  index_nonzero = (res_sample > 0)
  res_sample_nonzero = res_sample[index_nonzero]
  mu_sample_nonzero = mu_sample[index_nonzero]
  
  res_all =  allCon/exp(coeff[1] + coeff[2]*all_site + coeff[3]*all_len+ coeff[4]*all_cov)
  mu_all = exp(coeff[1] + coeff[2]*all_site + coeff[3]*all_len+ coeff[4]*all_cov)
}else{
  colnames(contig_info) = c('contig_name' , 'length')
  coverage = read.table(coverage_name , header = F , sep = '' , skip= 1)[ , c(1 , 6)]
  if(sum(contig_info[ , 1]!=coverage[ , 1])>0)
  {
    stop('Order mistakes exist in converage file')
  }
  contig_info$coverage = coverage[ , 2]
  
  sample_len = rep(0 , nrow(sample_data))
  sample_cov = rep(0 , nrow(sample_data))
  
  for(i in 1:nrow(sample_data))
  {
    
    sample_len[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 2]) * 
                          as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 2]))
    
    sample_cov[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 3]) * 
                          as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 3]))
  }
  
  sampleCon = as.numeric(sample_data[ , 3])
  
  mean_len = mean(sample_len)
  sd_len = sd(sample_len)
  mean_cov = mean(sample_cov)
  sd_cov = sd(sample_cov)
  
  sample_len = (sample_len-mean_len)/sd_len
  sample_cov = (sample_cov-mean_cov)/sd_cov
  
  data_sample = cbind(sample_len , sample_cov , sampleCon)
  data_sample = as.data.frame(data_sample)
  colnames(data_sample) = c('sample_len' , 'sample_cov' , 'sampleCon')
  
  
  all_len = rep(0 , nrow(all_contacts))
  all_cov = rep(0 , nrow(all_contacts))
  
  for(i in 1:nrow(all_contacts))
  {
    all_len[i] = log(as.numeric(contig_info[as.numeric(all_contacts[i , 1]) , 2]) * 
                       as.numeric(contig_info[as.numeric(all_contacts[i , 2]) , 2]))
    
    all_cov[i] = log(as.numeric(contig_info[as.numeric(all_contacts[i , 1]) , 3]) * 
                       as.numeric(contig_info[as.numeric(all_contacts[i , 2]) , 3]))
  }
  
  allCon = as.numeric(all_contacts[ , 3])
  all_len = (all_len-mean_len)/sd_len
  all_cov = (all_cov-mean_cov)/sd_cov
  
  
  tryCatch(
    {
      message(paste("normalizing",sep=" "))
      
      if(opt$unlabeled){
        fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                       ziformula=~1,family=nbinom2)
        
        
      }else if(opt$nb){
        fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                       ziformula=~0,family=nbinom2)
        
      }else{
        fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                       ziformula=~sample_len+sample_cov , family=nbinom2)
      }
    },
    error = function(e){
      message(e)
      message(paste("\nskip",  sep=" "))
    },
    warning = function(w){
      message(w)
      message(paste("\nskip",  sep=" "))
    }
  )
  
  
  coeff = as.numeric(fit1$fit$par)
  res_sample = sampleCon/exp(coeff[1] + coeff[2]*sample_len + coeff[3]*sample_cov)
  mu_sample = exp(coeff[1] + coeff[2]*sample_len + coeff[3]*sample_cov)
  index_nonzero = (res_sample > 0)
  res_sample_nonzero = res_sample[index_nonzero]
  mu_sample_nonzero = mu_sample[index_nonzero]
  
  res_all =  allCon/exp(coeff[1] + coeff[2]*all_len + coeff[3]*all_cov)
  mu_all = exp(coeff[1] + coeff[2]*all_len + coeff[3]*all_cov)
}


###########detect spurious contacts#################
sigma = summary(fit1)$sigma

pvalue_all = pnbinom(allCon , size = sigma  ,mu = mu_all)
pvalue_sample = pnbinom(sampleCon[index_nonzero] , size =  sigma  , mu = mu_sample_nonzero)

index_spur = (pvalue_all<quantile(pvalue_sample , thres) | res_all < quantile(res_sample_nonzero, thres))

all_contacts$contacts = res_all

res_all_valid = all_contacts[!index_spur , ]
res_all_spur = all_contacts[index_spur , ]


write.table(res_all_valid , file=output_file_name, row.names=F , col.names = F , sep = ',')

if(!is.na(discarded_contact_name)){
  write.table(res_all_spur , file=discarded_contact_name, row.names=F , col.names = F , sep = ',')
}




