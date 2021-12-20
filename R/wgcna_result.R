library(WGCNA)

x = read.csv(
  "data/Human__TCGA_KIRP__JHU_USC__Methylation__Meth450__01_28_2016__BI__Gene__Firehose_Methylation_Prepocessor.cct.gz",
  sep = "\t"
) %>%
  .[colMeans(apply(., 1, is.na)) < 0.3,]

kirp.gene = x$attrib_name[1:100]

x = x %>%
  filter(attrib_name %in% kirp.gene) %>%
  select(-attrib_name) %>%
  as.matrix(.) %>% 
  impute::impute.knn(.) %>% 
  .$data %>% 
  t() %>%
  as.matrix()

colnames(x) <- kirp.gene

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  x,
  dataIsExpr = TRUE,
  powerVector = powers,
  corFnc = cor,
  corOptions = list(use = 'p'),
  networkType = "unsigned"
)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


softPower = 3;

#calclute the adjacency matrix
adj= adjacency(x,type = "unsigned", power = softPower);
TOM=TOMsimilarityFromExpr(x,networkType = "unsigned", TOMType = "unsigned", power = softPower);

