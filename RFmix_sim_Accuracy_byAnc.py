#get_ipython().system(u'pwd')


import numpy as np
import string
import re


outFile = open('RFmix_simRep_Accuracy_byAnc.Bantu.txt', 'w') # to output the correctly called Bantu sites
outFile1 = open('RFmix_simRep_Accuracy_byAnc.French.txt', 'w') # to output the correctly called Bantu sites
outFile2 = open('RFmix_simRep_Accuracy_byAnc.KhoeSan.txt', 'w') # to output the correctly called Bantu sites


#RFmix haplotypes file:
#RFwins = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_khomani/SimSan_Rfmix_refAdmixed_autoAll.Viterbi1.bed", "r")
RFwins = open("SimSan_Rfmix_replicate_autoAll.Viterbi1.bed", "r")
wins = np.genfromtxt(RFwins, dtype=None)  #import the text of the windows file as a numpy array

#simulation bp pos true anc calls
#AncSites = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_khomani/MEGA_SimKhm1_chr1.AncCalls", "r")
AncSitesSim = open("MEGA_SimKhm1_allchr.AncCalls", "r")
AncSites = np.genfromtxt(AncSitesSim, dtype=None)


for line in range(0, len(wins)):
#look up SNPs in the sim file in terms of the RFmix windows.
# specifically, log or the start and stop of each haplotype window as well as each individuals ancestry in there
    
    Hapstart = wins[line][1]
    Hapstop = wins[line][2]
    chrom = wins[line][0]
    SimSan1_rf = wins[line][7]
    SimSan2_rf = wins[line][8]
    SimSan3_rf = wins[line][9]
    SimSan4_rf = wins[line][10]
    SimSan5_rf = wins[line][11]
    SimSan6_rf = wins[line][12]
    SimSan7_rf = wins[line][13]
    SimSan8_rf = wins[line][14]
    SimSan9_rf = wins[line][15]
    SimSan10_rf = wins[line][16]

    
    count = 0  #count is the number of SNPs in this window
    countTrue_1 = 0  #the number of SNPs in this window that have correct anc for the ind
    countTrue_2 = 0
    countTrue_3 = 0
    countTrue_4 = 0
    countTrue_5 = 0
    countTrue_6 = 0
    countTrue_7 = 0
    countTrue_8 = 0
    countTrue_9 = 0
    countTrue_10 = 0
    countKhoeSan_1 = 0 
    countKhoeSan_1 = 0
    countKhoeSan_2 = 0
    countKhoeSan_3 = 0
    countKhoeSan_4 = 0
    countKhoeSan_5 = 0
    countKhoeSan_6 = 0
    countKhoeSan_7 = 0
    countKhoeSan_8 = 0
    countKhoeSan_9 = 0
    countKhoeSan_10 = 0
    countBantu_1 = 0 
    countBantu_1 = 0
    countBantu_2 = 0
    countBantu_3 = 0
    countBantu_4 = 0
    countBantu_5 = 0
    countBantu_6 = 0
    countBantu_7 = 0
    countBantu_8 = 0
    countBantu_9 = 0
    countBantu_10 = 0
    countFrench_1 = 0 
    countFrench_1 = 0
    countFrench_2 = 0
    countFrench_3 = 0
    countFrench_4 = 0
    countFrench_5 = 0
    countFrench_6 = 0
    countFrench_7 = 0
    countFrench_8 = 0
    countFrench_9 = 0
    countFrench_10 = 0
    countBantu_snps = 0
    countFrench_snps = 0
    countKhoeSan_snps = 0
    countBantu_snps1 = 0
    countBantu_snps2 = 0
    countBantu_snps3 = 0
    countBantu_snps4 = 0
    countBantu_snps5 = 0
    countBantu_snps6 = 0
    countBantu_snps7 = 0
    countBantu_snps8 = 0
    countBantu_snps9 = 0
    countBantu_snps10 = 0
    countBantu_snps11 = 0
    countBantu_snps12 = 0
    countBantu_snps13 = 0
    countBantu_snps14 = 0
    countBantu_snps15 = 0
    countBantu_snps16 = 0
    countBantu_snps17 = 0
    countBantu_snps18 = 0
    countBantu_snps19 = 0
    countBantu_snps20 = 0
    countBantu_snps21 = 0
    countBantu_snps22 = 0
    countFrench_snps1 = 0
    countFrench_snps2 = 0
    countFrench_snps3 = 0
    countFrench_snps4 = 0
    countFrench_snps5 = 0
    countFrench_snps6 = 0
    countFrench_snps7 = 0
    countFrench_snps8 = 0
    countFrench_snps9 = 0
    countFrench_snps10 = 0
    countFrench_snps11 = 0
    countFrench_snps12 = 0
    countFrench_snps13 = 0
    countFrench_snps14 = 0
    countFrench_snps15 = 0
    countFrench_snps16 = 0
    countFrench_snps17 = 0
    countFrench_snps18 = 0
    countFrench_snps19 = 0
    countFrench_snps20 = 0
    countFrench_snps21 = 0
    countFrench_snps22 = 0
    countKhoeSan_snps1 = 0
    countKhoeSan_snps2 = 0
    countKhoeSan_snps3 = 0
    countKhoeSan_snps4 = 0
    countKhoeSan_snps5 = 0
    countKhoeSan_snps6 = 0
    countKhoeSan_snps7 = 0
    countKhoeSan_snps8 = 0
    countKhoeSan_snps9 = 0
    countKhoeSan_snps10 = 0
    countKhoeSan_snps11 = 0
    countKhoeSan_snps12 = 0
    countKhoeSan_snps13 = 0
    countKhoeSan_snps14 = 0
    countKhoeSan_snps15 = 0
    countKhoeSan_snps16 = 0
    countKhoeSan_snps17 = 0
    countKhoeSan_snps18 = 0
    countKhoeSan_snps19 = 0
    countKhoeSan_snps20 = 0
    countKhoeSan_snps21 = 0
    countKhoeSan_snps22 = 0    


    
    for line in range(0, len(AncSites)):
        siteChrom = AncSites[line][0]
        sitePos = AncSites[line][1]
        SimSan1_sim = AncSites[line][3] # since they're only 10 simulated individuals just load them in individually
        SimSan2_sim = AncSites[line][4]
        SimSan3_sim = AncSites[line][5]
        SimSan4_sim = AncSites[line][6]
        SimSan5_sim = AncSites[line][7]
        SimSan6_sim = AncSites[line][8]
        SimSan7_sim = AncSites[line][9]
        SimSan8_sim = AncSites[line][10]
        SimSan9_sim = AncSites[line][11]
        SimSan10_sim = AncSites[line][12]
        
        if chrom == siteChrom:   #screen to match chr
            if Hapstart <= sitePos < Hapstop:
                count = count + 1   #tally up the SNPs in the window
                if SimSan1_sim == 0:
                    countBantu_snps1 = countBantu_snps1 + 1 # Count how many snps in the window were bantu
                if SimSan1_sim == 1:
                    countFrench_snps1 = countFrench_snps1 + 1
                if SimSan1_sim == 2:
                    countKhoeSan_snps1 = countKhoeSan_snps1 + 1
                if SimSan1_rf == SimSan1_sim:
                    countTrue_1 = countTrue_1 + 1
                    if SimSan1_rf == 0:
                        countBantu_1 = countBantu_1 + 1 # then can get how many snps were correctly called for each ancestry
                    if SimSan1_rf == 1:
                        countFrench_1 = countFrench_1 + 1
                    if SimSan1_rf == 2:
                        countKhoeSan_1 = countKhoeSan_1 + 1
                
                if SimSan2_sim == 0:
                    countBantu_snps2 = countBantu_snps2 + 1 
                if SimSan2_sim == 1:
                    countFrench_snps2 = countFrench_snps2 + 1
                if SimSan2_sim == 2:
                    countKhoeSan_snps2 = countKhoeSan_snps2 + 1                
                if SimSan2_rf == SimSan2_sim:
                    countTrue_2 = countTrue_2 + 1
                    if SimSan2_rf == 0:
                        countBantu_2 = countBantu_2 + 1
                    if SimSan2_rf == 1:
                        countFrench_2 = countFrench_2 + 1
                    if SimSan2_rf == 2:
                        countKhoeSan_2 = countKhoeSan_2 + 1
                
                if SimSan3_sim == 0:
                    countBantu_snps3 = countBantu_snps3 + 1 
                if SimSan3_sim == 1:
                    countFrench_snps3 = countFrench_snps3 + 1
                if SimSan3_sim == 2:
                    countKhoeSan_snps3 = countKhoeSan_snps3 + 1                
                if SimSan3_rf == SimSan3_sim:
                    countTrue_3 = countTrue_3 + 1
                    if SimSan3_rf == 0:
                        countBantu_3 = countBantu_3 + 1
                    if SimSan3_rf == 1:
                        countFrench_3 = countFrench_3 + 1
                    if SimSan3_rf == 2:
                        countKhoeSan_3 = countKhoeSan_3 + 1
                
                if SimSan4_sim == 0:
                    countBantu_snps4 = countBantu_snps4 + 1 
                if SimSan4_sim == 1:
                    countFrench_snps4 = countFrench_snps4 + 1
                if SimSan4_sim == 2:
                    countKhoeSan_snps4 = countKhoeSan_snps4 + 1                
                if SimSan4_rf == SimSan4_sim:
                    countTrue_4 = countTrue_4 + 1
                    if SimSan4_rf == 0:
                        countBantu_4 = countBantu_4 + 1
                    if SimSan4_rf == 1:
                        countFrench_4 = countFrench_4 + 1
                    if SimSan4_rf == 2:
                        countKhoeSan_4 = countKhoeSan_4 + 1
                
                if SimSan5_sim == 0:
                    countBantu_snps5 = countBantu_snps5 + 1 
                if SimSan5_sim == 1:
                    countFrench_snps5 = countFrench_snps5 + 1
                if SimSan5_sim == 2:
                    countKhoeSan_snps5 = countKhoeSan_snps5 + 1                
                if SimSan5_rf == SimSan5_sim:
                    countTrue_5 = countTrue_5 + 1
                    if SimSan5_rf == 0:
                        countBantu_5 = countBantu_5 + 1
                    if SimSan5_rf == 1:
                        countFrench_5 = countFrench_5 + 1
                    if SimSan5_rf == 2:
                        countKhoeSan_5 = countKhoeSan_5 + 1
                
                if SimSan6_sim == 0:
                    countBantu_snps6 = countBantu_snps6 + 1 
                if SimSan6_sim == 1:
                    countFrench_snps6 = countFrench_snps6 + 1
                if SimSan6_sim == 2:
                    countKhoeSan_snps6 = countKhoeSan_snps6 + 1                
                if SimSan6_rf == SimSan6_sim:
                    countTrue_6 = countTrue_6 + 1
                    if SimSan6_rf == 0:
                        countBantu_6 = countBantu_6 + 1
                    if SimSan6_rf == 1:
                        countFrench_6 = countFrench_6 + 1
                    if SimSan6_rf == 2:
                        countKhoeSan_6 = countKhoeSan_6 + 1
                
                if SimSan7_sim == 0:
                    countBantu_snps7 = countBantu_snps7 + 1 
                if SimSan7_sim == 1:
                    countFrench_snps7 = countFrench_snps7 + 1
                if SimSan7_sim == 2:
                    countKhoeSan_snps7 = countKhoeSan_snps7 + 1                   
                if SimSan7_rf == SimSan7_sim:
                    countTrue_7 = countTrue_7 + 1
                    if SimSan7_rf == 0:
                        countBantu_7 = countBantu_7 + 1
                    if SimSan7_rf == 1:
                        countFrench_7 = countFrench_7 + 1
                    if SimSan7_rf == 2:
                        countKhoeSan_7 = countKhoeSan_7 + 1
             
                
                if SimSan8_sim == 0:
                    countBantu_snps8 = countBantu_snps8 + 1 
                if SimSan8_sim == 1:
                    countFrench_snps8 = countFrench_snps8 + 1
                if SimSan8_sim == 2:
                    countKhoeSan_snps8 = countKhoeSan_snps8 + 1                
                if SimSan8_rf == SimSan8_sim:
                    countTrue_8 = countTrue_8 + 1
                    if SimSan8_rf == 0:
                        countBantu_8 = countBantu_8 + 1
                    if SimSan8_rf == 1:
                        countFrench_8 = countFrench_8 + 1
                    if SimSan8_rf == 2:
                        countKhoeSan_8 = countKhoeSan_8 + 1
                
                if SimSan9_sim == 0:
                    countBantu_snps9 = countBantu_snps9 + 1 
                if SimSan9_sim == 1:
                    countFrench_snps9 = countFrench_snps9 + 1
                if SimSan9_sim == 2:
                    countKhoeSan_snps9 = countKhoeSan_snps9 + 1                
                if SimSan9_rf == SimSan9_sim:
                    countTrue_9 = countTrue_9 + 1
                    if SimSan9_rf == 0:
                        countBantu_9 = countBantu_9 + 1
                    if SimSan9_rf == 1:
                        countFrench_9 = countFrench_9 + 1
                    if SimSan9_rf == 2:
                        countKhoeSan_9 = countKhoeSan_9 + 1
                
                if SimSan10_sim == 0:
                    countBantu_snps10 = countBantu_snps10 + 1 
                if SimSan10_sim == 1:
                    countFrench_snps10 = countFrench_snps10 + 1
                if SimSan10_sim == 2:
                    countKhoeSan_snps10 = countKhoeSan_snps10 + 1                
                if SimSan10_rf == SimSan10_sim:
                    countTrue_10 = countTrue_10 + 1
                    if SimSan10_rf == 0:
                        countBantu_10 = countBantu_10 + 1
                    if SimSan10_rf == 1:
                        countFrench_10 = countFrench_10 + 1
                    if SimSan10_rf == 2:
                        countKhoeSan_10 = countKhoeSan_10 + 1



    #output_global = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countTrue_1),str(countTrue_2),str(countTrue_3),str(countTrue_4),str(countTrue_5),str(countTrue_6),str(countTrue_7),str(countTrue_8),str(countTrue_9),str(countTrue_10))
        #outputs the chromosome, the start and stop sites of the particular window, and true-called snp counts for all 10 simulated San
    output_Bantu = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countBantu_snps1), str(countBantu_1),str(countBantu_snps2), str(countBantu_2), str(countBantu_snps3), str(countBantu_3), str(countBantu_snps4), str(countBantu_4), str(countBantu_snps5), str(countBantu_5), str(countBantu_snps6), str(countBantu_6),str(countBantu_snps7), str(countBantu_7), str(countBantu_snps8), str(countBantu_8),str(countBantu_snps9), str(countBantu_9), str(countBantu_snps10), str(countBantu_10))    
    output_French = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countFrench_snps1), str(countFrench_1),str(countFrench_snps2), str(countFrench_2), str(countFrench_snps3), str(countFrench_3), str(countFrench_snps4), str(countFrench_4), str(countFrench_snps5), str(countFrench_5), str(countFrench_snps6), str(countFrench_6),str(countFrench_snps7), str(countFrench_7), str(countFrench_snps8), str(countFrench_8),str(countFrench_snps9), str(countFrench_9), str(countFrench_snps10), str(countFrench_10))    
    output_KhoeSan = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countKhoeSan_snps1), str(countKhoeSan_1),str(countKhoeSan_snps2), str(countKhoeSan_2), str(countKhoeSan_snps3), str(countKhoeSan_3), str(countKhoeSan_snps4), str(countKhoeSan_4), str(countKhoeSan_snps5), str(countKhoeSan_5), str(countKhoeSan_snps6), str(countKhoeSan_6),str(countKhoeSan_snps7), str(countKhoeSan_7), str(countKhoeSan_snps8), str(countKhoeSan_8),str(countKhoeSan_snps9), str(countKhoeSan_9), str(countKhoeSan_snps10), str(countKhoeSan_10))    
    outFile.write(output_Bantu)
    outFile1.write(output_French)
    outFile2.write(output_KhoeSan)
    
outFile.close()
outFile1.close()
outFile2.close()
  

