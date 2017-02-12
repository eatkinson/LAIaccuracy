#get_ipython().system(u'pwd')


import numpy as np
import string
import re


outFile = open('RFmix_sim_globalAccuracy.txt', 'w')


#RFmix haplotypes file:
#RFwins = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_khomani/SimSan_Rfmix_refAdmixed_autoAll.Viterbi.bed", "r")
RFwins = open("/Users/eatkinson/Dropbox/Projects/AGR/Simulation vs RFMix/KhoeSan_rep_Elizabeth/Accuracy/SimSan_Rfmix_refAdmixed_autoAll.Viterbi1.bed", "r")
wins = np.genfromtxt(RFwins, dtype=None)  #import the text of the windows file as a numpy array

#simulation bp pos true anc calls
#AncSites = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_khomani/MEGA_SimKhm1_chr1.AncCalls", "r")
AncSitesSim = open("/Users/eatkinson/Dropbox/Projects/AGR/Simulation vs RFMix/KhoeSan_rep_Elizabeth/Accuracy/MEGA_SimKhm1_allchr.AncCalls", "r")
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
                if SimSan1_rf == SimSan1_sim:
                    countTrue_1 = countTrue_1 + 1
                if SimSan2_rf == SimSan2_sim:
                    countTrue_2 = countTrue_2 + 1
                if SimSan3_rf == SimSan3_sim:
                    countTrue_3 = countTrue_3 + 1
                if SimSan4_rf == SimSan4_sim:
                    countTrue_4 = countTrue_4 + 1
                if SimSan5_rf == SimSan5_sim:
                    countTrue_5 = countTrue_5 + 1
                if SimSan6_rf == SimSan6_sim:
                    countTrue_6 = countTrue_6 + 1
                if SimSan7_rf == SimSan7_sim:
                    countTrue_7 = countTrue_7 + 1
                if SimSan8_rf == SimSan8_sim:
                    countTrue_8 = countTrue_8 + 1
                if SimSan9_rf == SimSan9_sim:
                    countTrue_9 = countTrue_9 + 1
                if SimSan10_rf == SimSan10_sim:
                    countTrue_10 = countTrue_10 + 1



    output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countTrue_1),str(countTrue_2),str(countTrue_3),str(countTrue_4),str(countTrue_5),str(countTrue_6),str(countTrue_7),str(countTrue_8),str(countTrue_9),str(countTrue_10))
        #outputs the chromosome, the start and stop sites of the particular window, and true-called snp counts for all 10 simulated San
        
    outFile.write(output)
outFile.close()
  

