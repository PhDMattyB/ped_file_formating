setwd('~/PhD/SNP Demographic Modelling/Working Directory')
library(LEA)
library(tidyverse)

sNMF_all_pops_large_run = snmf('icelandic_pops_sNMF_poly_largerun.geno',
                     K = 1:15, 
                     entropy = TRUE,
                     repetitions = 50,
                     project = 'new')

##To load the sNMF progject use:
#sNMF_all_pops = load.snmfProject('icelandic_pops_sNMF_poly.snmfProject')

#ex = export.snmfProject('icelandic_pops_sNMF_poly_largerun.snmfProject')

sNMF_all_pops_large_run = load.snmfProject('icelanidc_pops_sNMF_poly_largerun.snmfProject')

sNMF_all_pops_large_run = import.snmfProject('icelandic_pops_sNMF_poly_largerun.snmf')

sNMF_all_pops_large_run = import.snmfProject("icelandic_pops_sNMF_poly_largerun.zip", 'test')
##Find the K value with the lowest cross.entropy
show(sNMF_all_pops_large_run)
summary(sNMF_all_pops_large_run)

##use cross.entropy() function to get the CE for a K value

##Plot the CE values to determine K
plot(sNMF_all_pops_large_run)
##minimum cross entropy (CE) needed to find the best K
##update the large run K value
ce = cross.entropy(sNMF_all_pops_large_run, K = 10)
best = which.min(ce)

##BARPLOT SHOWING ADMIXTURE BETWEEN INDIVIDUALS AND POPULATIONS K = v??
#bp = barchart(sNMF_all_pops, K = 10, run = best, border = NA, space = 0,
#             col = my_pallet, xlab = 'Indviduals', ylab = 'Ancestry proportions', 
#            horiz = F, names.arg = snmf_names, las = 1, cex.names = 0.35, cex = 0.7)

#axis(1, at = 1:length(bp$order), labels = snmf_names, las = 1, cex.axis = 0.5)

##Names for the different plots
names = read_csv('Name_order_ped_file.csv')
#snmf_names = names$pop
snmf_lamorph = names$lamorph
#snmf_indiv = names$indiv

##Colour pallet
my_pallet = c('cadetblue3', 'chocolate2', 'brown3', 'burlywood3',
              'darkorchid3', 'seagreen3', 'skyblue4', 'springgreen4',
              'sienna3', 'plum')
#BP_pallet = c('cadetblue3', 'aquamarine4', 'lightpink2', 'lightsalmon3', 'indianred3',
              #'seagreen3', 'palegreen', 'burlywood3', 'brown3', 'plum')
# pwdth = 10
# pheight = 3
# plot.window(c(0, pwdth), c(0, pheight))


##Nicer barchart
bp = barchart(sNMF_all_pops_large_run, K = 10, run = best, border = NA, space = 0,
              col = my_pallet, xlab = 'Indviduals', ylab = 'Ancestry proportions',
              horiz = F, names.arg = snmf_lamorph, las = 1, cex.names = 0.35,
              cex = 0.7)

axis(1, at = 1:length(bp$order), labels = bp$order, las = 2, cex.axis = 0.4)
#    cex.axis = 0.5, srt = 45)
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##NEED TO MAKE A NICE COLOUR PALLET SO THE SNMF BARPLOT DOESNT SUCK!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


##Order of individuals corresponding the Qmatrix plot!
snmf_order_large = bp$order
write.csv(snmf_order_large, 'Order_sNMF_AllIcelanidc_large_polymorphic.csv')
#write.tsv(snmf_order, 'Order_sNMF_AllIcelandic_pops_polymorphic.txt')

ped = read_tsv('Sept192018_plink_input_icelandic_pops3.ped',col_names = F)
head(ped)
ped = ped[-1,]

ped_name_column = ped[,1:2]
write.csv(ped_name_column, 'Name_order_ped_file.csv')
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##Population differentiation test using sNMF !!!!!!!!!!!!!!!!!!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Find pavalues for SNP differentiation
snmf_pval = snmf.pvalues(sNMF_all_pops, entropy = T, ploidy = 2, 
                         K = 10)

pval = snmf_pval$pvalues
#par(mfrow = c(1,1))
hist(pval, col = 'seagreen3')
#axis(2, at = seq(0, 3000, 50))
plot(-log10(pval), pch = 19, col = 'black', cex = 0.4)
