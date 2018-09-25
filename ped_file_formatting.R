setwd('~/PhD/SNP Demographic modelling/Working Directory')
library(tidyverse)
library(magrittr)

##Read in tsv files

##Read in the map file for the snps in the ped file

map = read_tsv('fljotaa.map', col_names = T, skip = 3)
##Order the map file by the marker ID
map2 = map[order(map$`Marker ID`),]
##Rename the marker Id's to match the summary file so we can filter out monomorphic snps
colnames(map2) = c('Chromosome', 'probeset_id', 'Genetic distance', 'Physical distance')


##Read in the affymetrix id file with all of the markers
affy_id = read_tsv('array_markers.tsv', col_names = TRUE)
##Order the affymetrix marker file by the marker names
affy_id2 = affy_id[order(affy_id$Name),]

##Read in the chromosome index file that Cam gave me. 
##This indexes chromosome ids with what chromosome they belong on as well as 
##unplaced contigs
chr_index = read_tsv('chr_index.txt', skip = 30, col_names = T)
#chr_index$`RefSeq-Accn`

##Read in the summary file from affymetrix and filter out the polymorphic snps
poly = read_tsv('fljotaa_summary.txt', col_names = T) %>%
  filter(ConversionType == 'PolyHighResolution')
##Order the polymorphic snps by the marker names, so everything is ordered right
poly = poly[order(poly$probeset_id),]

##join the map and the polymorphic snp summary file by the matching name so that
##The map file only retains polymorphic snps
poly_snps_only = left_join(poly, map2, by = 'probeset_id') %>% 
  select(Chromosome, probeset_id, `Genetic distance`, `Physical distance`)
##Rename the columns to match downstreams
colnames(poly_snps_only) = c('Chromosome', 'Name', 'Genetic distance', 'Physical position')

##Join the polymorphic snp map file with the affymetrix id's
update_names = left_join(poly_snps_only, affy_id2, by = 'Name')%>%
  select(Chromosome.x, `Affx ID`, Name, `Genetic distance`, 
         `Physical position`)
##Rename the ID column to match the vcf file
colnames(update_names)[2] = 'ID'

sort_snps = update_names[order(update_names$ID),]

##read in vcf file with new positions
vcf = read_tsv('array_markers.NCBI.Salp.V2.vcf', skip = 3, col_names = T)
#vcf$`#CHROM`
sort_vcf = vcf[order(vcf$ID),]
#Merge the two files based on the common id column
merge_update_vcf = left_join(sort_snps, sort_vcf, by = 'ID')
#merge_update_vcf$`#CHROM`

##Select what we need
update_pos = merge_update_vcf %>% select(`#CHROM`, ID, `Genetic distance`,
                                         POS)

#update_pos = merge_update_vcf %>% select(`#CHROM`, ID, POS, `Physical position`)
colnames(update_pos) = '#Chromosome'
colnames(update_pos)[2:4] = c('Marker ID', 'Genetic distance', 'Physical position')
update_pos = na.omit(update_pos)

# update_pos$`#Chromosome`
# 
# View(update_pos)


##The chromosome dictionary this has each chromosome id matched to a chromosome 
##number
chrom_num = read_tsv('chromeosome_dict.txt', col_names = T)

update_num = left_join(update_pos, chrom_num, by = '#Chromosome') %>%
  select(ChromNum, `Marker ID`, `Genetic distance`,`Physical position`)
colnames(update_num) = c('#Chromosome', 'Marker ID','Genetic distance', 'Physical position')


write_tsv(update_num, 'updated_fljotaa.map', na = 'NA')


##.Ped file format #####

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP!  
##BEFORE THIS STEP YOU NEED TO RUN CAMS PSEUDO CHROMOSOME MAKER PYTHON SCRIPT!!!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##Read in the updated .map file
pseduo_chrome = read_tsv('updated_fljotaa_pseudochrome.map', col_names = TRUE) 
  


##Need to use this to assign the marker names to the ped file. 
header = c('Sample', pseduo_chrome$`Marker ID`)
##This is the correct ped file
ped = read_tsv('fljotaa.ped', col_names = header)
#ped = read_tsv('all_icelandic_snps.ped', col_names = F)
ped = ped[-1,] #%<>%
##Sort by sample column
sorted_ped = ped[order(ped$Sample),]


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP!!
##WE LEFT OFF CODING THE SCRIPT BETTER FOR AN EFFICIENT WORK FLOW!!!!!!!!!!!!!!!!!!!
##WE NEED TO MAKE THE FIRST 6 COLUMNS OF THE PED FILE AND THEN LOAD THIS MOTHER IN!!
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##The grouping for the morphs from the ped file with all icelandic populations:
##  1. T.LGB
##  2. V.BR
##  3. V.SIL
##  4. S.PL
##  5. S.PI
##  6. T.PL
##  7. T. SB
##  8. S.LGB
##  9. G.SB
##  10. G.PI
##  11. fljotaa
##  12. Fljotaa

names_of_the_ped = read_tsv('pedfile_sample_names.txt', col_names = TRUE) %>%
  group_by(FamilyID) %>%
  filter(FamilyID == '12')
# samples = ped %>% select(`#Sample Filename`)
# 
# write_tsv(samples, 'pedfile_sample_names.txt')


##Shortend names file
#names = read_csv('galtabol_ped_sample_names.csv', col_names = T)
colnames(names_of_the_ped) = c('#FamilyID', 'IndividualID', 'PaternalID', 'MaternalID',
                    'Sex', 'Phenotype', 'Sample', 'ShortName')
##Need to sort by sample
sorted_names = names_of_the_ped[order(names_of_the_ped$Sample),]
head(sorted_names$Sample)
colnames(sorted_names) = c('#FamilyID', 'IndividualID', 'PaternalID', 'MaternalID',
                    'Sex', 'Phenotype', 'Sample', 'ShortName')
head(sorted_names$ShortName)

short_ped = left_join(sorted_names, sorted_ped, by = 'Sample') %>%
  select(-Sample, -ShortName)
  
#short_ped = short_ped %>% select(-Sample, -ShortName)


write_tsv(short_ped, 'Sept242018_plink_input_fljotaa.ped')



#complete_ped = read_tsv('Sept192018_plink_input_icelandic_pops3.ped', col_names = T)



