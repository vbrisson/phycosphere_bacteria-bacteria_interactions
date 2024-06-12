# import data
size_bact <- read.table("data/biovolume.txt", header = TRUE, sep = ',')
size_pt <- read.table('data/size_alga.csv', header = TRUE, sep = ',')

# calculate biovolume from Pt
size_pt$Biovolume <- pi/12 * size_pt$Length * (size_pt$Width)^2

# calculate C mass from biovolume
size_bact$mass <- (88.6 * size_bact$Biovolume^0.59) * 0.86  # [fg] Mayali et al 2023 Nat Comm
size_bact$mass <- exp(4.28) * (size_bact$Biovolume^1.12)  # [fg] Fagerbakke et al 1996 Aquat Microb Ecol
size_pt$mass <- 288 * size_pt$Biovolume^0.811  # [fg] Menden-Deuer et al 2000 Limnol Oceangr

# export data
df <- size_bact
df[nrow(df)+1, 1] <- 'Alga'
df[4, 2] <- median(size_pt$Biovolume)
df[4, 3] <- median(size_pt$mass)
write.table(df, 'data/mass.csv', sep=',', row.names=FALSE)