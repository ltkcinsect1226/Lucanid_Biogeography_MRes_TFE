data<-read.csv("lucanid gene and biogeography data2.csv", head=T)
head(data)

##Sort out species with both 16S and 28S available
selectd<-read.csv("species selection.csv", head=T)
head(selectd)
selectd1<-selectd[-(262),]
selectd1[is.na(selectd1)] <- 0

a<-c()
for (i in 1:nrow(selectd1)) {
  if(selectd1$X16S[i]=="1" & selectd1$X28S[i]=="1"){
    a<-c(a, print(selectd1$Species[i]))
  }
}

a

## Steps in phylogenetic tree building Were done outside R. 

## DEC using Lagrange
# Install optimx
install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")
library(optimx)

# Also get snow (for parallel processing)
install.packages("snow")
library(snow)

# Install phylobase
install.packages("phylobase", dependencies=TRUE, repos="http://cran.rstudio.com")

# From October 2018 onwards, install BioGeoBEARS from GitHub:
# https://github.com/nmatzke/BioGeoBEARS
install.packages("devtools")
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")


# Load the package (after installation, see above).
install.packages("GenSA")
library(GenSA)    # GenSA is better than optimx (although somewhat slower)
install.packages("FD")
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
install.packages("parallel")
library(parallel)
installed.packages("rexpokit")
library(rexpokit)
installed.packages("CladoRcpp")
library(cladoRcpp)
library(BioGeoBEARS)

#set the wd to the folder with the newick trees and range txt

trfn =   "MCC wo outgroup.newick"     

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)


geogfn = "species dis3.txt"  
# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

#Find example data directory, then manually copy all files to your wd
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir # Now you can open the example dataset to know the correct formatting and file type etc.

myextdata_dir<- "C:/Users/Calvin T K Leung/Desktop/real transfer/calvin's stuff/Imperial/dissertation/lucanid phylogeny/sequences/Lagrange"

tree_file_name = np(paste(addslash(myextdata_dir), "MCC edited.newick", sep=""))
tr = read.tree(tree_file_name)

# plot tree
plot(tr)
axisPhylo()

#Remove non binary root because of outgroup deletion
library(ape)
multi2di(tr)
di2multi(tr, tol=0.00000001)
is.binary(tr)
write.tree(tr, file = "C:/Users/Calvin T K Leung/Desktop/real transfer/calvin's stuff/Imperial/dissertation/lucanid phylogeny/sequences/Lagrange/DECTr.newick", append = FALSE, digits = 17)

geo_file_name = np(paste(addslash(myextdata_dir), "species dis final.data", sep=""))
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges

#Setup DEC model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = tree_file_name # feed in tree
BioGeoBEARS_run_object$geogfn = geo_file_name # feed in range data
BioGeoBEARS_run_object$max_range_size = 6 #number of biogeographic realms
# The rest as default
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
#Run DEC!
results_DEC = bears_optim_run(BioGeoBEARS_run_object)
results_DEC


# Setup DEC+J model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name
BioGeoBEARS_run_object$max_range_size = 6
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# use d and e from DEC
dstart = results_DEC$outputs@params_table["d","est"]
estart = results_DEC$outputs@params_table["e","est"]
# set starting values
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
# Set J as free value
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
results_DECJ = bears_optim_run(BioGeoBEARS_run_object)

# set a pdf for DEC and DEC J comparison
pdffn = "lucanid_DEC_vs_DEC+J.pdf"
pdf(pdffn, width=12, height=28)

analysis_titletxt = "DEC Lucanidae"
results_object = results_DEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

analysis_titletxt ="DEC+J Lucanidae"
results_object = results_DECJ
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


#Calculate AICs using LnL
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(results_DEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(results_DECJ)

numparams1 = 3
numparams2 = 2

stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats$AIC1 #205.6598 
stats$AIC2 #279.9003
stats$pval #2.51142304575804e-18
## Therefore DEC+J is better

## Read relative probability of each state at each node
results_DECJ$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
write.csv(results_DECJ$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS, "C:/Users/Calvin T K Leung/Desktop/real transfer/calvin's stuff/Imperial/dissertation/lucanid phylogeny/rpbranchattop.csv")

## Now find ML of each node
results_DECJ$ML_marginal_prob_each_state_at_branch_top_AT_node
write.csv(results_DECJ$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS, "C:/Users/Calvin T K Leung/Desktop/real transfer/calvin's stuff/Imperial/dissertation/lucanid phylogeny/marginalML.csv")

## Testing age difference between tropical and temperate lineages
agedata<-read.csv("lineage age.csv",header = T)

#See how the age distribution wihin groups looks like
library(ggplot2)
ggplot(agedata, aes(x = age..Myr.)) +
  geom_histogram(aes(color = tropicortemperate), fill = "white",
                 position = "identity", bins = 15) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))

dev.off()

## Very discrete discrete distribution in temperate age... parametric test is almost impossible
## So use 2-sample Wilcoxon test
mod<-wilcox.test(age..Myr.~tropicortemperate, data = agedata) 
mod #W = 37, p-value = 0.1081, insignificant.

#Lucaninae Wilcoxon
luage<-subset(agedata, agedata$Lucanine.=="Y")
luage
mod1<-wilcox.test(age..Myr.~tropicortemperate, data = luage)
mod1 #W = 18, p-value = 0.149. Insignificant

#plot boxplot for age comparison
bx<-boxplot(age..Myr.~ tropicortemperate, data=agedata, ylab = "age (Myr)", xlab= "lineage type")
bx

ggbx <- ggplot(agedata, aes(x=tropicortemperate, y=age..Myr.)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  labs(title = "A", x="lineage type", y = "age (Mya)")+theme(
    plot.title = element_text(size=14, face="bold"),
  axis.title.x = element_text(size=10, face="bold"),
axis.title.y = element_text(size=10, face="bold"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggbx

## Lucaninae only boxplot

lggbx <- ggplot(luage, aes(x=tropicortemperate, y=age..Myr.)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  labs(title = "B", x="lineage type", y = "age (Mya)")+theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=10, face="bold"),
    axis.title.y = element_text(size=10, face="bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

lggbx

library(gridExtra)
library(grid)
grid.arrange(ggbx, lggbx, nrow = 1)

## Net direction of transition
td<-read.csv("transitions.csv", header=T)
td
x2<-table(td$original.realm, td$type.transition)
x2
chisq.test(x2) # X-squared = 1.0651, df = 1, p-value = 0.3021
## therefore no difference
## warning said 'Chi-squared approximation may be incorrect' because of small sample size. 
## Let's try with eaxt binomial test

binom.test(c(3,10),1/2) # p-value = 0.09229, 95CI 0.05038107 0.53813154
binom.test(c(4,3), 1/2) # p-value = 1, 95ci 0.1840516 0.9010117

##Lucaninae only
ltd<-td[-(1:5),] # first 5 rows are not lucanines
ltd
bn<-table(ltd$original.realm, ltd$type.transition)
bn
binom.test(c(3,9),1/2) # p-value = 0.146, 95CI 0.05486064 0.57185846
binom.test(c(2,1), 1/2) #p = 1, 95CI 0.09429932 0.99159624
