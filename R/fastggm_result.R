library(FastGGM)

x = read.csv(
  "data/Human__TCGA_KIRP__JHU_USC__Methylation__Meth450__01_28_2016__BI__Gene__Firehose_Methylation_Prepocessor.cct.gz",
  sep = "\t"
) %>%
  .[colMeans(apply(., 1, is.na)) < 0.3,]

kirp.gene = read.csv("data/gene_list.csv")

x = x %>%
  filter(attrib_name %in% kirp.gene) %>%
  select(-attrib_name) %>%
  as.matrix(.) %>% 
  impute::impute.knn(.) %>% 
  .$data %>% 
  t() %>%
  as.matrix()

colnames(x) <- kirp.gene

tuning = 
  foreach(lambda = exp(seq(log(0.01),log(length(kirp.gene)),len = 20)),
          .combine = rbind,
          .packages = c("FastGGM","tidyverse")) %do% {
            
            Omega = list()
            
            for(i in 1:20){
            y = x[sample(1:nrow(x),0.5*nrow(x)),]
            Omega[[i]] = FastGGM_Parallel(y,lambda = lambda)$precision
            }
            
            Stab = Reduce("+",lapply(Omega, function(x) x>0))/20
            Stab = sum(Stab*(1-Stab))/choose(ncol(Stab),2)
            
            tibble(lambda = lambda,
                   stab = Stab)
          }