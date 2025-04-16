library(dplyr)

#load tsv copied from https://drive.google.com/drive/folders/1NCqSf9C5ka964DEtmaYRemwaUONn4iUu by Patel et al.
#Or, you can find copied files here: https://drive.google.com/open?id=1yMw1pYloLfD4NqEFovFaMUYHLb1hFnX7&usp=drive_fs
df_1 <- read.delim("group1.tsv", header=TRUE, sep="\t")
df_2 <- read.delim("group2.tsv", header=TRUE, sep="\t")
df_3 <- read.delim("group3.tsv", header=TRUE, sep="\t")
df_meta <- read.delim("india_metadata.tsv", header=TRUE, sep="\t")

sample_size <- c(length(unique(df_1$SAMPLE_ID)), length(unique(df_2$SAMPLE_ID)), length(unique(df_3$SAMPLE_ID)))

sample_size

# merge df_1/2/3 and df_meta by SAMPLE_ID==Accession.ID to check age
df_1meta <- left_join(df_1, df_meta, by=c("SAMPLE_ID"="Accession.ID"))
unique(df_1meta$Patient.age)
df_2meta <- left_join(df_2, df_meta, by=c("SAMPLE_ID"="Accession.ID"))
unique(df_2meta$Patient.age)
df_3meta <- left_join(df_3, df_meta, by=c("SAMPLE_ID"="Accession.ID"))
unique(df_3meta$Patient.age)


#combine
df_all <- rbind(df_1, df_2, df_3)
df_all$group <- c(rep(1, nrow(df_1)), rep(2, nrow(df_2)), rep(3, nrow(df_3)))

df_all$mut <- paste0(df_all$REF, df_all$POS, df_all$ALT)

#list of all sub
unique_sub <- unique(df_all$mut)

#count from each group
count1 <- c()
count2 <- c()
count3 <- c()

for (i in unique_sub) {
  # count number of rows in df_all. group ==1, mut==i
  temp1 <- nrow(df_all %>% filter(group == 1, mut == i))
  count1 <- c(count1, temp1)
  temp2 <- nrow(df_all %>% filter(group == 2, mut == i))
  count2 <- c(count2, temp2)
  temp3 <- nrow(df_all %>% filter(group == 3, mut == i))
  count3 <- c(count3, temp3)
}



#count-table
df_count <- data.frame(mut = unique_sub, count1 = count1, count2 = count2, count3 = count3)

# keep only >=5 mutations
df_count_all <- df_count
df_count <- df_count %>% filter(count1 + count2 + count3 >= 5)

# count all subs in each group
nrow(df_count)
nrow(df_count %>% filter(count1 >0))
nrow(df_count %>% filter(count2 >0))
nrow(df_count %>% filter(count3 >0))


#extract group-unique seq: e.g. count1>=1, count2==0, count3==0
df_unique1 <- df_count %>% filter(count1 >= 1, count2 == 0, count3 == 0)
df_unique2 <- df_count %>% filter(count1 == 0, count2 >= 1, count3 == 0)
df_unique3 <- df_count %>% filter(count1 == 0, count2 == 0, count3 >= 1)

#combine
df_unique <- rbind(df_unique1, df_unique2, df_unique3)










#print unique subs and counts
#group1
df_unique1$mut
nrow(df_unique1)
sum(df_unique1$count1)

#group2
df_unique2$mut
nrow(df_unique2)
sum(df_unique2$count2)

#group3
df_unique3$mut
nrow(df_unique3)
sum(df_unique3$count3)










#subsampling probability to miss 5 sequences with a certain sub
# group 2 missing prop
sample_prop <- 23325 / 171363
sample_prop
original_sub_seq <- 1

allsample_size <- 171363
subsample_size <- round(allsample_size * sample_prop)
phyper(0, original_sub_seq, (allsample_size - original_sub_seq), subsample_size)






  




# test prevalence for all subs
pvalue_seq <- c()
for (i in 1:nrow(df_count)) {
  temp <- c(df_count[i, 2:4], (sample_size - df_count[i, 2:4]))
  temp <- as.numeric(temp)
  data <- matrix(temp, 
                 nrow=3, byrow=FALSE,
                 dimnames = list(Group=c("G1", "G2", "G3"),
                                 Outcome=c("sub", "NoSub")))
  chiseq <- chisq.test(data)
  pvalue_seq <- c(pvalue_seq, chiseq$p.value)
}

df_count$pvalue <- pvalue_seq
df_count$bpvalue <- ifelse(df_count$pvalue * nrow(df_count) < 1, df_count$pvalue * nrow(df_count), 1)




# combine to df_unique
df_unique <- left_join(df_unique, df_count, by="mut")

# print when pvalue < 0.05
df_unique_p5 <- df_unique %>% filter(pvalue < 0.05) %>% select(mut, count1.x, count2.x, count3.x, pvalue, bpvalue)
df_unique_p5

nrow(df_unique_p5 %>% filter(count1.x > 0))
nrow(df_unique_p5 %>% filter(count2.x > 0))
nrow(df_unique_p5 %>% filter(count3.x > 0))



#count total seqs with unique subs
# sum(df_unique_p5$count1.x)
# sum(df_unique_p5$count2.x)
# sum(df_unique_p5$count3.x)



# print when bpvalue < 0.05
df_unique_bp5 <- df_unique %>% filter(bpvalue < 0.05) %>% select(mut, count1.x, count2.x, count3.x, pvalue, bpvalue)
df_unique_bp5

nrow(df_unique_bp5 %>% filter(count1.x > 0))
nrow(df_unique_bp5 %>% filter(count2.x > 0))
nrow(df_unique_bp5 %>% filter(count3.x > 0))


#count total seqs with unique subs
# sum(df_unique_bp5$count1.x)
# sum(df_unique_bp5$count2.x)
# sum(df_unique_bp5$count3.x)










# unique subs not necessarily 0 in 2 groups
df_count_sig <- df_count %>% filter(bpvalue < 0.05) %>% select(mut, count1, count2, count3, pvalue, bpvalue)
df_count_sig
nrow(df_count_sig)
