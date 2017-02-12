
# coding: utf-8

# In[ ]:

#get_ipython().system(u'pwd')


import numpy as np
import string
import re


outFile = open('EastAf_2k_RFmix_Acc_Bantu.txt', 'w') # to output the correctly called Bantu sites
outFile1 = open('EastAf_2k_RFmix_Acc_Nilotic.txt', 'w') # to output the correctly called Bantu sites


#RFmix haplotypes file:
RFwins = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_AGR/EastAfrican/EastAfSim_2k_Rfmix_chr_autoAll.Viterbi.bed", "r")
wins = np.genfromtxt(RFwins, dtype=None)  #import the text of the windows file as a numpy array

#simulation bp pos true anc calls
AncSitesSim = open("/vault/henn/people/elizabeth/ADRP/RFmix/simulations_AGR/EastAfrican/EastAfSim_2k_Rfmix_allchr.AncCalls", "r")
AncSites = np.genfromtxt(AncSitesSim, dtype=None)

for line in range(0, len(wins)):
#look up SNPs in the sim file in terms of the RFmix windows.
# specifically, log  the start and stop of each haplotype window as well as each individuals ancestry in there
    
    Hapstart = wins[line][1]
    Hapstop = wins[line][2]
    chrom = wins[line][0]
    SimInd1_rf = wins[line][7]
    SimInd2_rf = wins[line][8]
    SimInd3_rf = wins[line][9]
    SimInd4_rf = wins[line][10]
    SimInd5_rf = wins[line][11]
    SimInd6_rf = wins[line][12]
    SimInd7_rf = wins[line][13]
    SimInd8_rf = wins[line][14]
    SimInd9_rf = wins[line][15]
    SimInd10_rf = wins[line][16]
    SimInd11_rf = wins[line][17]
    SimInd12_rf = wins[line][18]
    SimInd13_rf = wins[line][19]
    SimInd14_rf = wins[line][20]
    SimInd15_rf = wins[line][21]
    SimInd16_rf = wins[line][22]
    SimInd17_rf = wins[line][23]
    SimInd18_rf = wins[line][24]
    SimInd19_rf = wins[line][25]
    SimInd20_rf = wins[line][26]
    SimInd21_rf = wins[line][27]
    SimInd22_rf = wins[line][28]
    SimInd23_rf = wins[line][29]
    SimInd24_rf = wins[line][30]
    SimInd25_rf = wins[line][31]

    
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
    countBantu_11 = 0  
    countBantu_12 = 0
    countBantu_13 = 0
    countBantu_14 = 0
    countBantu_15 = 0
    countBantu_16 = 0
    countBantu_17 = 0
    countBantu_18 = 0
    countBantu_19 = 0
    countBantu_20 = 0
    countBantu_21 = 0  
    countBantu_22 = 0
    countBantu_23 = 0
    countBantu_24 = 0
    countBantu_25 = 0
    
    countNilotic_1 = 0 
    countNilotic_1 = 0
    countNilotic_2 = 0
    countNilotic_3 = 0
    countNilotic_4 = 0
    countNilotic_5 = 0
    countNilotic_6 = 0
    countNilotic_7 = 0
    countNilotic_8 = 0
    countNilotic_9 = 0
    countNilotic_10 = 0
    countNilotic_11 = 0  
    countNilotic_12 = 0
    countNilotic_13 = 0
    countNilotic_14 = 0
    countNilotic_15 = 0
    countNilotic_16 = 0
    countNilotic_17 = 0
    countNilotic_18 = 0
    countNilotic_19 = 0
    countNilotic_20 = 0
    countNilotic_21 = 0  
    countNilotic_22 = 0
    countNilotic_23 = 0
    countNilotic_24 = 0
    countNilotic_25 = 0
    
    countBantu_snps = 0
    countNilotic_snps = 0
    
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
    countBantu_snps23 = 0
    countBantu_snps24 = 0
    countBantu_snps25 = 0

    
    countNilotic_snps1 = 0
    countNilotic_snps2 = 0
    countNilotic_snps3 = 0
    countNilotic_snps4 = 0
    countNilotic_snps5 = 0
    countNilotic_snps6 = 0
    countNilotic_snps7 = 0
    countNilotic_snps8 = 0
    countNilotic_snps9 = 0
    countNilotic_snps10 = 0
    countNilotic_snps11 = 0
    countNilotic_snps12 = 0
    countNilotic_snps13 = 0
    countNilotic_snps14 = 0
    countNilotic_snps15 = 0
    countNilotic_snps16 = 0
    countNilotic_snps17 = 0
    countNilotic_snps18 = 0
    countNilotic_snps19 = 0
    countNilotic_snps20 = 0
    countNilotic_snps21 = 0
    countNilotic_snps22 = 0
    countNilotic_snps23 = 0
    countNilotic_snps24 = 0
    countNilotic_snps25 = 0


    
    for line in range(0, len(AncSites)):
        siteChrom = AncSites[line][0]
        sitePos = AncSites[line][1]
        SimInd1_sim = AncSites[line][3] # since they're only 10 simulated individuals just load them in individually
        SimInd2_sim = AncSites[line][4]
        SimInd3_sim = AncSites[line][5]
        SimInd4_sim = AncSites[line][6]
        SimInd5_sim = AncSites[line][7]
        SimInd6_sim = AncSites[line][8]
        SimInd7_sim = AncSites[line][9]
        SimInd8_sim = AncSites[line][10]
        SimInd9_sim = AncSites[line][11]
        SimInd10_sim = AncSites[line][12]
        SimInd11_sim = AncSites[line][3] # since they're only 10 simulated individuals just load them in individually
        SimInd12_sim = AncSites[line][4]
        SimInd13_sim = AncSites[line][5]
        SimInd14_sim = AncSites[line][6]
        SimInd15_sim = AncSites[line][7]
        SimInd16_sim = AncSites[line][8]
        SimInd17_sim = AncSites[line][9]
        SimInd18_sim = AncSites[line][10]
        SimInd19_sim = AncSites[line][11]
        SimInd20_sim = AncSites[line][12] # since they're only 10 simulated individuals just load them in individually
        SimInd21_sim = AncSites[line][3] # since they're only 10 simulated individuals just load them in individually
        SimInd22_sim = AncSites[line][4]
        SimInd23_sim = AncSites[line][5]
        SimInd24_sim = AncSites[line][6]
        SimInd25_sim = AncSites[line][7]

        
        if chrom == siteChrom:   #screen to match chr
            if Hapstart <= sitePos < Hapstop:
                count = count + 1   #tally up the SNPs in the window
                if SimInd1_sim == 0:
                    countBantu_snps1 = countBantu_snps1 + 1 # Count how many snps in the window were bantu
                if SimInd1_sim == 1:
                    countNilotic_snps1 = countNilotic_snps1 + 1
                if SimInd1_rf == SimInd1_sim:
                    countTrue_1 = countTrue_1 + 1
                    if SimInd1_rf == 0:
                        countBantu_1 = countBantu_1 + 1 # then can get how many snps were correctly called for each ancestry
                    if SimInd1_rf == 1:
                        countNilotic_1 = countNilotic_1 + 1
                
                if SimInd2_sim == 0:
                    countBantu_snps2 = countBantu_snps2 + 1 
                if SimInd2_sim == 1:
                    countNilotic_snps2 = countNilotic_snps2 + 1      
                if SimInd2_rf == SimInd2_sim:
                    countTrue_2 = countTrue_2 + 1
                    if SimInd2_rf == 0:
                        countBantu_2 = countBantu_2 + 1
                    if SimInd2_rf == 1:
                        countNilotic_2 = countNilotic_2 + 1
                
                if SimInd3_sim == 0:
                    countBantu_snps3 = countBantu_snps3 + 1 
                if SimInd3_sim == 1:
                    countNilotic_snps3 = countNilotic_snps3 + 1         
                if SimInd3_rf == SimInd3_sim:
                    countTrue_3 = countTrue_3 + 1
                    if SimInd3_rf == 0:
                        countBantu_3 = countBantu_3 + 1
                    if SimInd3_rf == 1:
                        countNilotic_3 = countNilotic_3 + 1
                
                if SimInd4_sim == 0:
                    countBantu_snps4 = countBantu_snps4 + 1 
                if SimInd4_sim == 1:
                    countNilotic_snps4 = countNilotic_snps4 + 1       
                if SimInd4_rf == SimInd4_sim:
                    countTrue_4 = countTrue_4 + 1
                    if SimInd4_rf == 0:
                        countBantu_4 = countBantu_4 + 1
                    if SimInd4_rf == 1:
                        countNilotic_4 = countNilotic_4 + 1
                
                if SimInd5_sim == 0:
                    countBantu_snps5 = countBantu_snps5 + 1 
                if SimInd5_sim == 1:
                    countNilotic_snps5 = countNilotic_snps5 + 1            
                if SimInd5_rf == SimInd5_sim:
                    countTrue_5 = countTrue_5 + 1
                    if SimInd5_rf == 0:
                        countBantu_5 = countBantu_5 + 1
                    if SimInd5_rf == 1:
                        countNilotic_5 = countNilotic_5 + 1
                
                if SimInd6_sim == 0:
                    countBantu_snps6 = countBantu_snps6 + 1 
                if SimInd6_sim == 1:
                    countNilotic_snps6 = countNilotic_snps6 + 1
                if SimInd6_rf == SimInd6_sim:
                    countTrue_6 = countTrue_6 + 1
                    if SimInd6_rf == 0:
                        countBantu_6 = countBantu_6 + 1
                    if SimInd6_rf == 1:
                        countNilotic_6 = countNilotic_6 + 1
                
                if SimInd7_sim == 0:
                    countBantu_snps7 = countBantu_snps7 + 1 
                if SimInd7_sim == 1:
                    countNilotic_snps7 = countNilotic_snps7 + 1         
                if SimInd7_rf == SimInd7_sim:
                    countTrue_7 = countTrue_7 + 1
                    if SimInd7_rf == 0:
                        countBantu_7 = countBantu_7 + 1
                    if SimInd7_rf == 1:
                        countNilotic_7 = countNilotic_7 + 1
                
                if SimInd8_sim == 0:
                    countBantu_snps8 = countBantu_snps8 + 1 
                if SimInd8_sim == 1:
                    countNilotic_snps8 = countNilotic_snps8 + 1        
                if SimInd8_rf == SimInd8_sim:
                    countTrue_8 = countTrue_8 + 1
                    if SimInd8_rf == 0:
                        countBantu_8 = countBantu_8 + 1
                    if SimInd8_rf == 1:
                        countNilotic_8 = countNilotic_8 + 1
                
                if SimInd9_sim == 0:
                    countBantu_snps9 = countBantu_snps9 + 1 
                if SimInd9_sim == 1:
                    countNilotic_snps9 = countNilotic_snps9 + 1            
                if SimInd9_rf == SimInd9_sim:
                    countTrue_9 = countTrue_9 + 1
                    if SimInd9_rf == 0:
                        countBantu_9 = countBantu_9 + 1
                    if SimInd9_rf == 1:
                        countNilotic_9 = countNilotic_9 + 1
                
                if SimInd10_sim == 0:
                    countBantu_snps10 = countBantu_snps10 + 1 
                if SimInd10_sim == 1:
                    countNilotic_snps10 = countNilotic_snps10 + 1              
                if SimInd10_rf == SimInd10_sim:
                    countTrue_10 = countTrue_10 + 1
                    if SimInd10_rf == 0:
                        countBantu_10 = countBantu_10 + 1
                    if SimInd10_rf == 1:
                        countNilotic_10 = countNilotic_10 + 1

                if SimInd11_sim == 0:
                    countBantu_snps11 = countBantu_snps11 + 1 
                if SimInd11_sim == 1:
                    countNilotic_snps11 = countNilotic_snps11 + 1              
                if SimInd11_rf == SimInd11_sim:
                    countTrue_11 = countTrue_11 + 1
                    if SimInd11_rf == 0:
                        countBantu_11 = countBantu_11 + 1
                    if SimInd11_rf == 1:
                        countNilotic_11 = countNilotic_11 + 1

                if SimInd12_sim == 0:
                    countBantu_snps12 = countBantu_snps12 + 1 
                if SimInd12_sim == 1:
                    countNilotic_snps12 = countNilotic_snps12 + 1              
                if SimInd12_rf == SimInd12_sim:
                    countTrue_12 = countTrue_12 + 1
                    if SimInd12_rf == 0:
                        countBantu_12 = countBantu_12 + 1
                    if SimInd12_rf == 1:
                        countNilotic_12 = countNilotic_12 + 1

                if SimInd13_sim == 0:
                    countBantu_snps13 = countBantu_snps13 + 1 
                if SimInd13_sim == 1:
                    countNilotic_snps13 = countNilotic_snps13 + 1              
                if SimInd13_rf == SimInd13_sim:
                    countTrue_13 = countTrue_13 + 1
                    if SimInd13_rf == 0:
                        countBantu_13 = countBantu_13 + 1
                    if SimInd13_rf == 1:
                        countNilotic_13 = countNilotic_13 + 1

                if SimInd14_sim == 0:
                    countBantu_snps14 = countBantu_snps14 + 1 
                if SimInd14_sim == 1:
                    countNilotic_snps14 = countNilotic_snps14 + 1              
                if SimInd14_rf == SimInd14_sim:
                    countTrue_14 = countTrue_14 + 1
                    if SimInd14_rf == 0:
                        countBantu_14 = countBantu_14 + 1
                    if SimInd14_rf == 1:
                        countNilotic_14 = countNilotic_14 + 1
                        
                if SimInd15_sim == 0:
                    countBantu_snps15 = countBantu_snps15 + 1 
                if SimInd15_sim == 1:
                    countNilotic_snps15 = countNilotic_snps15 + 1              
                if SimInd15_rf == SimInd15_sim:
                    countTrue_15 = countTrue_15 + 1
                    if SimInd15_rf == 0:
                        countBantu_15 = countBantu_15 + 1
                    if SimInd15_rf == 1:
                        countNilotic_15 = countNilotic_15 + 1

                if SimInd16_sim == 0:
                    countBantu_snps16 = countBantu_snps16 + 1 
                if SimInd16_sim == 1:
                    countNilotic_snps16 = countNilotic_snps16 + 1              
                if SimInd16_rf == SimInd16_sim:
                    countTrue_16 = countTrue_16 + 1
                    if SimInd16_rf == 0:
                        countBantu_16 = countBantu_16 + 1
                    if SimInd16_rf == 1:
                        countNilotic_16 = countNilotic_16 + 1

                if SimInd17_sim == 0:
                    countBantu_snps17 = countBantu_snps17 + 1 
                if SimInd17_sim == 1:
                    countNilotic_snps17 = countNilotic_snps17 + 1              
                if SimInd17_rf == SimInd17_sim:
                    countTrue_17 = countTrue_17 + 1
                    if SimInd17_rf == 0:
                        countBantu_17 = countBantu_17 + 1
                    if SimInd17_rf == 1:
                        countNilotic_17 = countNilotic_17 + 1
                        
                if SimInd18_sim == 0:
                    countBantu_snps18 = countBantu_snps18 + 1 
                if SimInd18_sim == 1:
                    countNilotic_snps18 = countNilotic_snps18 + 1              
                if SimInd18_rf == SimInd18_sim:
                    countTrue_18 = countTrue_18 + 1
                    if SimInd18_rf == 0:
                        countBantu_18 = countBantu_18 + 1
                    if SimInd18_rf == 1:
                        countNilotic_18 = countNilotic_18 + 1
                        
                if SimInd19_sim == 0:
                    countBantu_snps19 = countBantu_snps19 + 1 
                if SimInd19_sim == 1:
                    countNilotic_snps19 = countNilotic_snps19 + 1              
                if SimInd19_rf == SimInd19_sim:
                    countTrue_19 = countTrue_19 + 1
                    if SimInd19_rf == 0:
                        countBantu_19 = countBantu_19 + 1
                    if SimInd19_rf == 1:
                        countNilotic_19 = countNilotic_19 + 1
                                               
                if SimInd20_sim == 0:
                    countBantu_snps20 = countBantu_snps20 + 1 
                if SimInd20_sim == 1:
                    countNilotic_snps20 = countNilotic_snps20 + 1              
                if SimInd20_rf == SimInd20_sim:
                    countTrue_20 = countTrue_20 + 1
                    if SimInd20_rf == 0:
                        countBantu_20 = countBantu_20 + 1
                    if SimInd20_rf == 1:
                        countNilotic_20 = countNilotic_20 + 1

                if SimInd21_sim == 0:
                    countBantu_snps21 = countBantu_snps21 + 1 
                if SimInd21_sim == 1:
                    countNilotic_snps21 = countNilotic_snps21 + 1              
                if SimInd21_rf == SimInd21_sim:
                    countTrue_21 = countTrue_21 + 1
                    if SimInd21_rf == 0:
                        countBantu_21 = countBantu_21 + 1
                    if SimInd21_rf == 1:
                        countNilotic_21 = countNilotic_21 + 1

                if SimInd22_sim == 0:
                    countBantu_snps22 = countBantu_snps22 + 1 
                if SimInd22_sim == 1:
                    countNilotic_snps22 = countNilotic_snps22 + 1              
                if SimInd22_rf == SimInd22_sim:
                    countTrue_22 = countTrue_22 + 1
                    if SimInd22_rf == 0:
                        countBantu_22 = countBantu_22 + 1
                    if SimInd22_rf == 1:
                        countNilotic_22 = countNilotic_22 + 1

                if SimInd23_sim == 0:
                    countBantu_snps23 = countBantu_snps23 + 1 
                if SimInd23_sim == 1:
                    countNilotic_snps23 = countNilotic_snps23 + 1              
                if SimInd23_rf == SimInd23_sim:
                    countTrue_23 = countTrue_23 + 1
                    if SimInd23_rf == 0:
                        countBantu_23 = countBantu_23 + 1
                    if SimInd23_rf == 1:
                        countNilotic_23 = countNilotic_23 + 1

                if SimInd24_sim == 0:
                    countBantu_snps24 = countBantu_snps24 + 1 
                if SimInd24_sim == 1:
                    countNilotic_snps24 = countNilotic_snps24 + 1              
                if SimInd24_rf == SimInd24_sim:
                    countTrue_24 = countTrue_24 + 1
                    if SimInd24_rf == 0:
                        countBantu_24 = countBantu_24 + 1
                    if SimInd24_rf == 1:
                        countNilotic_24 = countNilotic_24 + 1
                        
                        
                if SimInd25_sim == 0:
                    countBantu_snps25 = countBantu_snps25 + 1 
                if SimInd25_sim == 1:
                    countNilotic_snps25 = countNilotic_snps25 + 1              
                if SimInd25_rf == SimInd25_sim:
                    countTrue_25 = countTrue_25 + 1
                    if SimInd25_rf == 0:
                        countBantu_25 = countBantu_25 + 1
                    if SimInd25_rf == 1:
                        countNilotic_25 = countNilotic_25 + 1
                                                
                        
    #output_Bantu = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countBantu_snps1), str(countBantu_1),str(countBantu_snps2), str(countBantu_2), str(countBantu_snps3), str(countBantu_3), str(countBantu_snps4), str(countBantu_4), str(countBantu_snps5), str(countBantu_5), str(countBantu_snps6), str(countBantu_6),str(countBantu_snps7), str(countBantu_7), str(countBantu_snps8), str(countBantu_8),str(countBantu_snps9), str(countBantu_9), str(countBantu_snps10), str(countBantu_10)), str(countBantu_snps11), str(countBantu_11),str(countBantu_snps12), str(countBantu_12), str(countBantu_snps13), str(countBantu_13), str(countBantu_snps14), str(countBantu_14), str(countBantu_snps15), str(countBantu_15), str(countBantu_snps16), str(countBantu_16),str(countBantu_snps17), str(countBantu_17), str(countBantu_snps18), str(countBantu_18),str(countBantu_snps19), str(countBantu_19), str(countBantu_snps20), str(countBantu_20), str(countBantu_snps21), str(countBantu_21),str(countBantu_snps22), str(countBantu_22), str(countBantu_snps23), str(countBantu_23), str(countBantu_snps24), str(countBantu_24), str(countBantu_snps25), str(countBantu_25))
    #output_Nilotic = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countNilotic_snps1), str(countNilotic_1),str(countNilotic_snps2), str(countNilotic_2), str(countNilotic_snps3), str(countNilotic_3), str(countNilotic_snps4), str(countNilotic_4), str(countNilotic_snps5), str(countNilotic_5), str(countNilotic_snps6), str(countNilotic_6),str(countNilotic_snps7), str(countNilotic_7), str(countNilotic_snps8), str(countNilotic_8),str(countNilotic_snps9), str(countNilotic_9), str(countNilotic_snps10), str(countNilotic_10)), str(countNilotic_snps11), str(countNilotic_11),str(countNilotic_snps12), str(countNilotic_12), str(countNilotic_snps13), str(countNilotic_13), str(countNilotic_snps14), str(countNilotic_14), str(countNilotic_snps15), str(countNilotic_15), str(countNilotic_snps16), str(countNilotic_16),str(countNilotic_snps17), str(countNilotic_17), str(countNilotic_snps18), str(countNilotic_18),str(countNilotic_snps19), str(countNilotic_19), str(countNilotic_snps20), str(countNilotic_20), str(countNilotic_snps21), str(countNilotic_21),str(countNilotic_snps22), str(countNilotic_22), str(countNilotic_snps23), str(countNilotic_23), str(countNilotic_snps24), str(countNilotic_24), str(countNilotic_snps25), str(countNilotic_25))
    #output_Bantu = "'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}','{19}','{20}','{21}','{22}','{23}','{24}','{25}','{26}','{27}','{28}','{29}','{30}','{31}','{32}','{33}','{34}','{35}','{36}','{37}','{38}','{39}','{40}','{41}','{42}','{43}','{44}','{45}','{46}','{47}','{48}','{49}','{50}','{51}','{52}','{53}'" .format(str(chrom), str(Hapstart), str(Hapstop), str(count), str(countBantu_snps1), str(countBantu_1),str(countBantu_snps2), str(countBantu_2), str(countBantu_snps3), str(countBantu_3), str(countBantu_snps4), str(countBantu_4), str(countBantu_snps5), str(countBantu_5), str(countBantu_snps6), str(countBantu_6),str(countBantu_snps7), str(countBantu_7), str(countBantu_snps8), str(countBantu_8),str(countBantu_snps9), str(countBantu_9), str(countBantu_snps10), str(countBantu_10)), str(countBantu_snps11), str(countBantu_11),str(countBantu_snps12), str(countBantu_12), str(countBantu_snps13), str(countBantu_13), str(countBantu_snps14), str(countBantu_14), str(countBantu_snps15), str(countBantu_15), str(countBantu_snps16), str(countBantu_16),str(countBantu_snps17), str(countBantu_17), str(countBantu_snps18), str(countBantu_18),str(countBantu_snps19), str(countBantu_19), str(countBantu_snps20), str(countBantu_20), str(countBantu_snps21), str(countBantu_21),str(countBantu_snps22), str(countBantu_22), str(countBantu_snps23), str(countBantu_23), str(countBantu_snps24), str(countBantu_24), str(countBantu_snps25), str(countBantu_25))
    #output_Bantu = "'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}','{19}','{20}','{21}','{22}','{23}','{24}','{25}','{26}','{27}','{28}','{29}','{30}','{31}','{32}','{33}','{34}','{35}','{36}','{37}','{38}','{39}','{40}','{41}','{42}','{43}','{44}','{45}','{46}','{47}','{48}','{49}','{50}','{51}','{52}','{53}'" .format(chrom, Hapstart, Hapstop, count, countBantu_snps1, countBantu_1,countBantu_snps2, countBantu_2, countBantu_snps3, countBantu_3, countBantu_snps4, countBantu_4, countBantu_snps5, countBantu_5, countBantu_snps6, countBantu_6,countBantu_snps7, countBantu_7, countBantu_snps8, countBantu_8,countBantu_snps9, countBantu_9, countBantu_snps10, countBantu_10, countBantu_snps11, countBantu_11,countBantu_snps12, countBantu_12, countBantu_snps13, countBantu_13, countBantu_snps14, countBantu_14, countBantu_snps15, countBantu_15, countBantu_snps16, countBantu_16,countBantu_snps17, countBantu_17, countBantu_snps18, countBantu_18,countBantu_snps19, countBantu_19, countBantu_snps20, countBantu_20, countBantu_snps21, countBantu_21,countBantu_snps22, countBantu_22, countBantu_snps23, countBantu_23, countBantu_snps24, countBantu_24, countBantu_snps25, countBantu_25)
    #output_Nilotic = "'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}','{19}','{20}','{21}','{22}','{23}','{24}','{25}','{26}','{27}','{28}','{29}','{30}','{31}','{32}','{33}','{34}','{35}','{36}','{37}','{38}','{39}','{40}','{41}','{42}','{43}','{44}','{45}','{46}','{47}','{48}','{49}','{50}','{51}','{52}','{53}'" .format(chrom, Hapstart, Hapstop, count, countNilotic_snps1, countNilotic_1,countNilotic_snps2, countNilotic_2, countNilotic_snps3, countNilotic_3, countNilotic_snps4, countNilotic_4, countNilotic_snps5, countNilotic_5, countNilotic_snps6, countNilotic_6,countNilotic_snps7, countNilotic_7, countNilotic_snps8, countNilotic_8,countNilotic_snps9, countNilotic_9, countNilotic_snps10, countNilotic_10, countNilotic_snps11, countNilotic_11,countNilotic_snps12, countNilotic_12, countNilotic_snps13, countNilotic_13, countNilotic_snps14, countNilotic_14, countNilotic_snps15, countNilotic_15, countNilotic_snps16, countNilotic_16,countNilotic_snps17, countNilotic_17, countNilotic_snps18, countNilotic_18,countNilotic_snps19, countNilotic_19, countNilotic_snps20, countNilotic_20, countNilotic_snps21, countNilotic_21,countNilotic_snps22, countNilotic_22, countNilotic_snps23, countNilotic_23, countNilotic_snps24, countNilotic_24, countNilotic_snps25, countNilotic_25)
    
    #output as a tuple, a bit cleaner to see and keep track of how many fields I am outputting. This way will output things as comma delimited 
    #output_Bantu = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53} \n" .format(chrom, Hapstart, Hapstop, count, countBantu_snps1, countBantu_1,countBantu_snps2, countBantu_2, countBantu_snps3, countBantu_3, countBantu_snps4, countBantu_4, countBantu_snps5, countBantu_5, countBantu_snps6, countBantu_6,countBantu_snps7, countBantu_7, countBantu_snps8, countBantu_8,countBantu_snps9, countBantu_9, countBantu_snps10, countBantu_10, countBantu_snps11, countBantu_11,countBantu_snps12, countBantu_12, countBantu_snps13, countBantu_13, countBantu_snps14, countBantu_14, countBantu_snps15, countBantu_15, countBantu_snps16, countBantu_16,countBantu_snps17, countBantu_17, countBantu_snps18, countBantu_18,countBantu_snps19, countBantu_19, countBantu_snps20, countBantu_20, countBantu_snps21, countBantu_21,countBantu_snps22, countBantu_22, countBantu_snps23, countBantu_23, countBantu_snps24, countBantu_24, countBantu_snps25, countBantu_25)
    #output_Nilotic = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53} \n" .format(chrom, Hapstart, Hapstop, count, countNilotic_snps1, countNilotic_1,countNilotic_snps2, countNilotic_2, countNilotic_snps3, countNilotic_3, countNilotic_snps4, countNilotic_4, countNilotic_snps5, countNilotic_5, countNilotic_snps6, countNilotic_6,countNilotic_snps7, countNilotic_7, countNilotic_snps8, countNilotic_8,countNilotic_snps9, countNilotic_9, countNilotic_snps10, countNilotic_10, countNilotic_snps11, countNilotic_11,countNilotic_snps12, countNilotic_12, countNilotic_snps13, countNilotic_13, countNilotic_snps14, countNilotic_14, countNilotic_snps15, countNilotic_15, countNilotic_snps16, countNilotic_16,countNilotic_snps17, countNilotic_17, countNilotic_snps18, countNilotic_18,countNilotic_snps19, countNilotic_19, countNilotic_snps20, countNilotic_20, countNilotic_snps21, countNilotic_21,countNilotic_snps22, countNilotic_22, countNilotic_snps23, countNilotic_23, countNilotic_snps24, countNilotic_24, countNilotic_snps25, countNilotic_25) 
    output_Bantu = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\t{35}\t{36}\t{37}\t{38}\t{39}\t{40}\t{41}\t{42}\t{43}\t{44}\t{45}\t{46}\t{47}\t{48}\t{49}\t{50}\t{51}\t{52}\t{53} \n" .format(chrom, Hapstart, Hapstop, count, countBantu_snps1, countBantu_1,countBantu_snps2, countBantu_2, countBantu_snps3, countBantu_3, countBantu_snps4, countBantu_4, countBantu_snps5, countBantu_5, countBantu_snps6, countBantu_6,countBantu_snps7, countBantu_7, countBantu_snps8, countBantu_8,countBantu_snps9, countBantu_9, countBantu_snps10, countBantu_10, countBantu_snps11, countBantu_11,countBantu_snps12, countBantu_12, countBantu_snps13, countBantu_13, countBantu_snps14, countBantu_14, countBantu_snps15, countBantu_15, countBantu_snps16, countBantu_16,countBantu_snps17, countBantu_17, countBantu_snps18, countBantu_18,countBantu_snps19, countBantu_19, countBantu_snps20, countBantu_20, countBantu_snps21, countBantu_21,countBantu_snps22, countBantu_22, countBantu_snps23, countBantu_23, countBantu_snps24, countBantu_24, countBantu_snps25, countBantu_25)
    output_Nilotic = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\t{35}\t{36}\t{37}\t{38}\t{39}\t{40}\t{41}\t{42}\t{43}\t{44}\t{45}\t{46}\t{47}\t{48}\t{49}\t{50}\t{51}\t{52}\t{53} \n" .format(chrom, Hapstart, Hapstop, count, countNilotic_snps1, countNilotic_1,countNilotic_snps2, countNilotic_2, countNilotic_snps3, countNilotic_3, countNilotic_snps4, countNilotic_4, countNilotic_snps5, countNilotic_5, countNilotic_snps6, countNilotic_6,countNilotic_snps7, countNilotic_7, countNilotic_snps8, countNilotic_8,countNilotic_snps9, countNilotic_9, countNilotic_snps10, countNilotic_10, countNilotic_snps11, countNilotic_11,countNilotic_snps12, countNilotic_12, countNilotic_snps13, countNilotic_13, countNilotic_snps14, countNilotic_14, countNilotic_snps15, countNilotic_15, countNilotic_snps16, countNilotic_16,countNilotic_snps17, countNilotic_17, countNilotic_snps18, countNilotic_18,countNilotic_snps19, countNilotic_19, countNilotic_snps20, countNilotic_20, countNilotic_snps21, countNilotic_21,countNilotic_snps22, countNilotic_22, countNilotic_snps23, countNilotic_23, countNilotic_snps24, countNilotic_24, countNilotic_snps25, countNilotic_25) 
 

    outFile.write(output_Bantu)
    outFile1.write(output_Nilotic)
    
outFile.close()
outFile1.close()
  



# In[ ]:



