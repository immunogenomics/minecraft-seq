####################################################################################################
## Script purpose: Process PAX5 amplicon data with allele tables
## Author: Zepeng Mu
## Date: Thu Jan 16 13:24:23 2025
####################################################################################################
library(data.table)
library(tidyverse)
library(pbmcapply)
library(cowplot)
library(ggridges)
library(PrettyCols)
library(RColorBrewer)
library(Phoenix)
library(ComplexHeatmap)
"%&%" <- function(a, b) paste0(a, b)

# Define functions ====
calcDeriv <- function(x) {
  len <- length(x)
  
  if (len == 1) return(0)
  
  x <- sort(x, decreasing = T)
  outX <- rep(0, len)
  outX[1:(len - 1)] <- x[2:len] - x[1:(len - 1)]
  return(outX)
}

idxSort <- fread("../idxSort.txt")

# Region A ====
allCellsA <- fread("../allCells_allele_A.txt")

targetA <- "GGACATGGAGGAGTGAATCAGCTTGGGGGGGTTTTTGTGAATG"

allCellsASumm <- allCellsA %>% 
  mutate(Target_Sequence = str_extract(Aligned_Sequence, str_replace_all(targetA, fixed("T"), "[TC]")),
         Reference_Sequence = str_extract(Reference_Sequence, targetA)) %>% 
  filter(!is.na(Target_Sequence) & !is.na(Reference_Sequence)) %>% 
  dplyr::group_by(smp, Target_Sequence) %>% 
  summarise(`#Reads` = sum(`#Reads`), `%Reads` = sum(`%Reads`), Reference_Sequence = unique(Reference_Sequence)) %>% 
  arrange(desc(`#Reads`)) %>% 
  mutate(rowNumber = row_number()) %>% 
  ungroup()

ggplot(allCellsASumm, aes(rowNumber, `#Reads`)) +
  theme_zm() +
  geom_point() +
  facet_wrap(~smp, scales = "free")

allCellsAFlt <- allCellsA %>% 
  mutate(
    Target_Sequence = str_extract(Aligned_Sequence, str_replace_all(targetA, fixed("T"), "[TC]")),
    Reference_Sequence = str_extract(Reference_Sequence, targetA)
  ) %>% 
  filter(!is.na(Target_Sequence) & !is.na(Reference_Sequence)) %>% 
  dplyr::group_by(smp, bc, Target_Sequence) %>% 
  summarise(`#Reads` = sum(`#Reads`), Reference_Sequence = unique(Reference_Sequence)) %>% 
  mutate(A.T_84 = str_sub(Target_Sequence, 14, 14),
         A.T_105 = str_sub(Target_Sequence, 35, 35),
         A.T_106 = str_sub(Target_Sequence, 36, 36),
         A.T_108 = str_sub(Target_Sequence, 38, 38),
         haplotype = str_glue("{A.T_84}{A.T_105}{A.T_106}{A.T_108}")) %>% 
  ungroup()

allCellsAFltHaplo %>% 
  left_join(idxSort, by = c("smp", "bc")) %>% 
  filter(cellType == "NTC") %>% 
  dplyr::group_by(smp, haplotype) %>% 
  tally() %>% 
  ggplot(aes(x = haplotype, y = n, fill = smp)) +
  theme_zm() +
  geom_col(position = position_dodge())

allCellsAFltGeno <- allCellsAFltHaplo %>% 
  dplyr::group_by(smp, bc) %>% 
  slice_max(nReads, n = 2) %>% 
  mutate(hapName = "hap"%&%row_number()) %>% 
  pivot_wider(id_cols = c(smp, bc, totalReads), names_from = hapName, values_from = haplotype) %>% 
  mutate(hap2 = case_when(is.na(hap2) ~ hap1, T ~ hap2),
         A.T_84 = (str_sub(hap1, 1, 1) == "C") + (str_sub(hap2, 1, 1) == "C"),
         A.T_105 = (str_sub(hap1, 2, 2) == "C") + (str_sub(hap2, 2, 2) == "C"),
         A.T_106 = (str_sub(hap1, 3, 3) == "C") + (str_sub(hap2, 3, 3) == "C"),
         A.T_108 = (str_sub(hap1, 4, 4) == "C") + (str_sub(hap2, 4, 4) == "C"))

fwrite(allCellsAFltGeno, col.names = T, row.names = F, sep = "\t", quote = F,
       file = "../allCells_haplo_geno_A.txt.gz")

allCellsABMGeno <- fread("../allCells_summ_genoNew.txt.gz")

allCellsAFltGeno <- allCellsAFltGeno %>% 
  left_join(allCellsABMGeno %>% dplyr::select(smp, bc, A.T_84.geno:A.T_108.geno), by = c("smp", "bc"))

allCellsAFltGeno %>% 
  filter(totalReads >= 100) %>% 
  dplyr::group_by(smp, A.T_84, A.T_84.geno) %>% 
  tally() %>% 
  ggplot(aes(factor(A.T_84), factor(A.T_84.geno), size = n)) +
  theme_zm() +
  geom_point() +
  scale_size_continuous(range = c(2, 12)) +
  facet_wrap(~ smp) +
  theme(aspect.ratio = 1)

allCellsAFltHaplo1 <- allCellsAFlt %>% 
  dplyr::select(-Target_Sequence, -Reference_Sequence, -(A.T_84:A.T_108)) %>% 
  dplyr::group_by(smp, bc, haplotype) %>% 
  summarise(nReads = sum(`#Reads`)) %>% 
  dplyr::group_by(smp, bc) %>% 
  mutate(totalReads = sum(nReads), propReads = nReads / totalReads) %>% 
  dplyr::arrange(desc(propReads), .by_group = T) %>% 
  mutate(cumPropReads = cumsum(propReads),
         cumCut = cumsum(cumPropReads >= 0.9)) %>% 
  filter(cumCut <= 1)

allCellsAFltHaplo1 %>% 
  left_join(idxSort, by = c("smp", "bc")) %>% 
  dplyr::group_by(smp, bc, cellType) %>% 
  tally() %>% 
  ggplot(aes(x = factor(n), fill = smp)) +
  theme_zm() +
  geom_bar(position = position_dodge()) +
  facet_wrap(~ cellType, scales = "free")

allCellsAFltHaplo1 %>% 
  left_join(idxSort, by = c("smp", "bc")) %>% 
  filter(cellType == "NTC") %>% 
  dplyr::group_by(smp, haplotype) %>% 
  tally() %>% 
  ggplot(aes(x = smp, y = n, fill = haplotype)) +
  theme_zm() +
  geom_col() +
  scale_y_continuous(expand = c(0, 0))

allCellsAFltHaploNTCAUC <- lapply(seq(0.05, 1, 0.02), function(s) {
  tmp <- allCellsAFlt %>% 
    dplyr::select(-Target_Sequence, -Reference_Sequence, -(A.T_84:A.T_108)) %>% 
    dplyr::group_by(smp, bc, haplotype) %>% 
    dplyr::summarise(nReads = sum(`#Reads`), .groups = "keep") %>% 
    dplyr::group_by(smp, bc) %>% 
    mutate(totalReads = sum(nReads), propReads = nReads / totalReads) %>% 
    dplyr::arrange(desc(propReads), .by_group = T) %>% 
    mutate(cumPropReads = cumsum(propReads),
           cumCut = cumsum(cumPropReads >= s)) %>% 
    filter(cumCut <= 1) %>% 
    left_join(idxSort, by = c("smp", "bc")) %>% 
    filter(totalReads >= 100) %>% 
    mutate(callNTC = length(unique(haplotype)) == 1 & unique(haplotype) == "TTTT") %>% 
    distinct(smp, bc, callNTC, cellType)
  
  tp <- sum(tmp$cellType == "NTC" & tmp$callNTC)
  tn <- sum(tmp$cellType != "NTC" & !tmp$callNTC)
  fp <- sum(tmp$cellType != "NTC" & tmp$callNTC)
  fn <- sum(tmp$cellType == "NTC" & !tmp$callNTC)
  
  recall <- tp / (tp + fn)
  precision <- tp / (tp + fp)
  
  return(data.frame(s, precision, recall, fp, tp))
}) %>% base::Reduce(rbind, .)

ggplot(allCellsAFltHaploNTCAUC, aes(precision, recall)) +
  geom_point() +
  geom_line()

ggplot(allCellsAFltHaploNTCAUC, aes(fp, tp)) +
  geom_point() +
  geom_line()

targetB <- "ATGACACCGTGCCTAGCGTCAGTTCCATCAACA"

# Region B ====
allCellsB <- fread("../allCells_allele_B.txt")

allCellsBFlt <- allCellsB %>% 
  mutate(
    Target_Sequence = str_extract(Aligned_Sequence, str_replace_all(targetB, fixed("T"), "[TC]")),
    Reference_Sequence = str_extract(Reference_Sequence, targetB)
  ) %>% 
  filter(!is.na(Target_Sequence) & !is.na(Reference_Sequence)) %>% 
  dplyr::group_by(smp, bc, Target_Sequence) %>% 
  summarise(`#Reads` = sum(`#Reads`), Reference_Sequence = unique(Reference_Sequence)) %>% 
  mutate(B.T_64 = str_sub(Target_Sequence, 14, 14),
         B.T_73 = str_sub(Target_Sequence, 23, 23),
         B.T_74 = str_sub(Target_Sequence, 24, 24),
         B.T_78 = str_sub(Target_Sequence, 28, 28),
         haplotype = str_glue("{B.T_64}{B.T_73}{B.T_74}{B.T_78}")) %>% 
  ungroup()

allCellsBFltHaplo <- allCellsBFlt %>% 
  dplyr::select(-Target_Sequence, -Reference_Sequence, -(B.T_64:B.T_78)) %>% 
  dplyr::group_by(smp, bc, haplotype) %>% 
  summarise(nReads = sum(`#Reads`)) %>% 
  ungroup() %>% 
  filter(nReads >= 10) %>% 
  dplyr::group_by(smp, bc) %>% 
  dplyr::arrange(desc(nReads), .by_group = T) %>% 
  mutate(totalReads = sum(nReads), propReads = nReads / totalReads,
         readDeriv = calcDeriv(log10(nReads)),
         isMin = cumsum(cumsum(readDeriv == min(readDeriv)))) %>% 
  filter(isMin <= 1)

allCellsBFltHaplo %>% 
  left_join(idxSort, by = c("smp", "bc")) %>% 
  dplyr::group_by(smp, bc, cellType) %>% 
  tally() %>% 
  ggplot(aes(x = factor(n), fill = smp)) +
  theme_zm() +
  geom_bar(position = position_dodge()) +
  facet_wrap(~ cellType, scales = "free")

allCellsBFltHaplo %>% 
  left_join(idxSort, by = c("smp", "bc")) %>% 
  filter(cellType == "NTC") %>% 
  dplyr::group_by(smp, haplotype) %>% 
  tally() %>% 
  ggplot(aes(x = haplotype, y = n, fill = smp)) +
  theme_zm() +
  geom_col(position = position_dodge())

allCellsBFltGeno <- allCellsBFltHaplo %>% 
  dplyr::group_by(smp, bc) %>% 
  slice_max(nReads, n = 2) %>% 
  mutate(hapName = "hap"%&%row_number()) %>% 
  pivot_wider(id_cols = c(smp, bc, totalReads), names_from = hapName, values_from = haplotype) %>% 
  mutate(hap2 = case_when(is.na(hap2) ~ hap1, T ~ hap2),
         B.T_64 = (str_sub(hap1, 1, 1) == "C") + (str_sub(hap2, 1, 1) == "C"),
         B.T_73 = (str_sub(hap1, 2, 2) == "C") + (str_sub(hap2, 2, 2) == "C"),
         B.T_74 = (str_sub(hap1, 3, 3) == "C") + (str_sub(hap2, 3, 3) == "C"),
         B.T_78 = (str_sub(hap1, 4, 4) == "C") + (str_sub(hap2, 4, 4) == "C"))

fwrite(allCellsBFltGeno, col.names = T, row.names = F, sep = "\t", quote = F,
       file = "../allCells_haplo_geno_B.txt.gz")

allCellsBFltGeno <- allCellsBFltGeno %>% 
  left_join(allCellsABMGeno %>% dplyr::select(smp, bc, B.T_64.geno:B.T_64.geno), by = c("smp", "bc"))

allCellsBFltGeno %>% 
  filter(totalReads >= 100) %>% 
  dplyr::group_by(smp, B.T_64, B.T_64.geno) %>% 
  tally() %>% 
  ggplot(aes(factor(B.T_64), factor(B.T_64.geno), size = n)) +
  theme_zm() +
  geom_point() +
  scale_size_continuous(range = c(2, 12)) +
  facet_wrap(~ smp) +
  theme(aspect.ratio = 1)

allCellsBFltHaploNTCAUC <- lapply(seq(0.05, 1, 0.02), function(s) {
  tmp <- allCellsBFlt %>% 
    dplyr::select(-Target_Sequence, -Reference_Sequence, -(B.T_64:B.T_78)) %>% 
    dplyr::group_by(smp, bc, haplotype) %>% 
    dplyr::summarise(nReads = sum(`#Reads`), .groups = "keep") %>% 
    dplyr::group_by(smp, bc) %>% 
    mutate(totalReads = sum(nReads), propReads = nReads / totalReads) %>% 
    dplyr::arrange(desc(propReads), .by_group = T) %>% 
    mutate(cumPropReads = cumsum(propReads),
           cumCut = cumsum(cumPropReads >= s)) %>% 
    filter(cumCut <= 1) %>% 
    left_join(idxSort, by = c("smp", "bc")) %>% 
    filter(totalReads >= 100) %>% 
    mutate(callNTC = length(unique(haplotype)) == 1 & unique(haplotype) == "TTTT") %>% 
    distinct(smp, bc, callNTC, cellType)
  
  tp <- sum(tmp$cellType == "NTC" & tmp$callNTC)
  tn <- sum(tmp$cellType != "NTC" & !tmp$callNTC)
  fp <- sum(tmp$cellType != "NTC" & tmp$callNTC)
  fn <- sum(tmp$cellType == "NTC" & !tmp$callNTC)
  
  recall <- tp / (tp + fn)
  precision <- tp / (tp + fp)
  
  return(data.frame(s, precision, recall, fp, tp))
}) %>% base::Reduce(rbind, .)

ggplot(allCellsBFltHaploNTCAUC, aes(fp / max(fp), tp / max(tp))) +
  theme_zm() +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_point(data = allCellsAFltHaploNTCAUC, aes(fp / max(fp), tp / max(tp))) +
  geom_line(data = allCellsAFltHaploNTCAUC, aes(fp / max(fp), tp / max(tp))) +
  theme(aspect.ratio = 1)

# Combine A and B ====
allA <- fread("../allCells_haplo_geno_A.txt.gz")
allB <- fread("../allCells_haplo_geno_B.txt.gz")

cmbHaploGeno <- allA %>% 
  dplyr::rename(A.totalReads = totalReads, A.hap1 = hap1, A.hap2 = hap2) %>% 
  dplyr::full_join(
    allB %>% dplyr::rename(B.totalReads = totalReads, B.hap1 = hap1, B.hap2 = hap2) %>% 
      dplyr::select(-hap3),
    by = c("bc", "smp")) %>% 
  dplyr::left_join(idxSort, by = c("bc", "smp"))

fwrite(cmbHaploGeno, col.names = T, row.names = F, sep = "\t", quote = F,
       file = "../allCells_haplo_geno_AB.txt.gz")
