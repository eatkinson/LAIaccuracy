
# coding: utf-8

# In[ ]:

#get_ipython().system(u'pwd')
#####calculating global accuracy comparing simulations to RFmix ancestry calls

##formatted for East African 2k ancestry calls, 2 ancestries and 25 individuals
###0 call = Bantu, 1 call = Nilotic

import numpy as np
import string
import re


outFile = open('EastAfrican_2kRFmix_globalAccuracy1.txt', 'w')


#RFmix haplotypes file:
RFwins = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_AGR/EastAfrican/EastAfSim_2k_Rfmix_chr_autoAll.Viterbi.bed", "r")
wins = np.genfromtxt(RFwins, dtype=None)  #import the text of the windows file as a numpy array

#simulation bp pos true anc calls
AncSitesSim = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_AGR/EastAfrican/EastAfSim_2k_Rfmix_allchr.AncCalls", "r")
#AncSitesSim = open("/Users/eatkinson/Dropbox/Projects/AGR/Simulation vs RFMix/KhoeSan_rep_Elizabeth/Accuracy/MEGA_SimKhm1_allchr.AncCalls", "r")
AncSites = np.genfromtxt(AncSitesSim, dtype=None)


for line in range(0, len(wins)):
#look up SNPs in the sim file in terms of the RFmix windows.
# specifically, log  the start and stop of each haplotype window as well as each individuals ancestry in there
    
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
    SimSan11_rf = wins[line][17]
    SimSan12_rf = wins[line][18]
    SimSan13_rf = wins[line][19]
    SimSan14_rf = wins[line][20]
    SimSan15_rf = wins[line][21]
    SimSan16_rf = wins[line][22]
    SimSan17_rf = wins[line][23]
    SimSan18_rf = wins[line][24]
    SimSan19_rf = wins[line][25]
    SimSan20_rf = wins[line][26]
    SimSan21_rf = wins[line][27]
    SimSan22_rf = wins[line][28]
    SimSan23_rf = wins[line][29]
    SimSan24_rf = wins[line][30]
    SimSan25_rf = wins[line][31]

    
    count = 0  #count is the number of SNPs in this window. 
    #also tally up the number of SNPs in this window that have correct anc for the ind
    countTrue_1 = 0  
    countTrue_2 = 0
    countTrue_3 = 0
    countTrue_4 = 0
    countTrue_5 = 0
    countTrue_6 = 0
    countTrue_7 = 0
    countTrue_8 = 0
    countTrue_9 = 0
    countTrue_10 = 0
    countTrue_11 = 0  
    countTrue_12 = 0
    countTrue_13 = 0
    countTrue_14 = 0
    countTrue_15 = 0
    countTrue_16 = 0
    countTrue_17 = 0
    countTrue_18 = 0
    countTrue_19 = 0
    countTrue_20 = 0
    countTrue_21 = 0  
    countTrue_22 = 0
    countTrue_23 = 0
    countTrue_24 = 0
    countTrue_25 = 0

    
    for line in range(0, len(AncSites)):
        siteChrom = AncSites[line][0]
        sitePos = AncSites[line][1]
        #  load them in the individuals individually... should automate this maybe
        SimSan1_sim = AncSites[line][3] 
        SimSan2_sim = AncSites[line][4]
        SimSan3_sim = AncSites[line][5]
        SimSan4_sim = AncSites[line][6]
        SimSan5_sim = AncSites[line][7]
        SimSan6_sim = AncSites[line][8]
        SimSan7_sim = AncSites[line][9]
        SimSan8_sim = AncSites[line][10]
        SimSan9_sim = AncSites[line][11]
        SimSan10_sim = AncSites[line][12]
        SimSan11_sim = AncSites[line][13] 
        SimSan12_sim = AncSites[line][14]
        SimSan13_sim = AncSites[line][15]
        SimSan14_sim = AncSites[line][16]
        SimSan15_sim = AncSites[line][17]
        SimSan16_sim = AncSites[line][18]
        SimSan17_sim = AncSites[line][19]
        SimSan18_sim = AncSites[line][20]
        SimSan19_sim = AncSites[line][21]
        SimSan20_sim = AncSites[line][22]
        SimSan21_sim = AncSites[line][23] 
        SimSan22_sim = AncSites[line][24]
        SimSan23_sim = AncSites[line][25]
        SimSan24_sim = AncSites[line][26]
        SimSan25_sim = AncSites[line][27]
        
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
                if SimSan11_rf == SimSan11_sim:
                    countTrue_11 = countTrue_11 + 1
                if SimSan12_rf == SimSan12_sim:
                    countTrue_12 = countTrue_12 + 1
                if SimSan13_rf == SimSan13_sim:
                    countTrue_13 = countTrue_13 + 1
                if SimSan14_rf == SimSan14_sim:
                    countTrue_14 = countTrue_14 + 1
                if SimSan15_rf == SimSan15_sim:
                    countTrue_15 = countTrue_15 + 1
                if SimSan16_rf == SimSan16_sim:
                    countTrue_16 = countTrue_16 + 1
                if SimSan17_rf == SimSan17_sim:
                    countTrue_17 = countTrue_17 + 1
                if SimSan18_rf == SimSan18_sim:
                    countTrue_18 = countTrue_18 + 1
                if SimSan19_rf == SimSan19_sim:
                    countTrue_19 = countTrue_19 + 1
                if SimSan20_rf == SimSan20_sim:
                    countTrue_20 = countTrue_20 + 1
                if SimSan21_rf == SimSan21_sim:
                    countTrue_21 = countTrue_21 + 1
                if SimSan22_rf == SimSan22_sim:
                    countTrue_22 = countTrue_22 + 1
                if SimSan23_rf == SimSan23_sim:
                    countTrue_23 = countTrue_23 + 1
                if SimSan24_rf == SimSan24_sim:
                    countTrue_24 = countTrue_24 + 1
                if SimSan25_rf == SimSan25_sim:
                    countTrue_25 = countTrue_25 + 1



    output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countTrue_1),str(countTrue_2),str(countTrue_3),str(countTrue_4),str(countTrue_5),str(countTrue_6),str(countTrue_7),str(countTrue_8),str(countTrue_9),str(countTrue_10),str(countTrue_11),str(countTrue_12),str(countTrue_13),str(countTrue_14),str(countTrue_15),str(countTrue_16),str(countTrue_17),str(countTrue_18),str(countTrue_19),str(countTrue_20),str(countTrue_21),str(countTrue_22),str(countTrue_23),str(countTrue_24),str(countTrue_25))
        #outputs the chromosome, the start and stop sites of the particular window, and true-called snp counts for all 10 simulated San
        
    outFile.write(output)
outFile.close()
  



# In[ ]:



