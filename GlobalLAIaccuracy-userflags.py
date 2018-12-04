#get_ipython().system(u'pwd')
#####calculating global accuracy comparing simulations to RFmix ancestry calls
__author__ = 'egatkinson'

##formatted for Grady AA ancestry calls, 2 ancestries and 25 individuals
import argparse
import numpy as np
import string
import re
import pprint


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--Anc', help='path stem to simulated truth ancestry calls', required=True) 
  parser.add_argument('--msp', help='path stem to RFmix ancestry calls msp file', required=True)
  args = parser.parse_args()
  return(args)    

args = parse_args()

#true anc calls from simulations
AncSitesSim = open(args.Anc)
AncSites = np.genfromtxt(AncSitesSim, dtype=None)

#RFmix haplotypes file:
RFwins = open(args.msp)
wins = np.genfromtxt(RFwins, dtype=None)  #import the text of the windows file as a numpy array

outFile = open(args.Anc + '.globalAccuracy.txt', 'w')

for line in range(0, len(wins)):
#look up SNPs in the sim file in terms of the RFmix windows.
# specifically, log  the start and stop of each haplotype window as well as each individuals ancestry in there
    
    Hapstart = wins[line][1]
    Hapstop = wins[line][2]
    chrom = wins[line][0]
    SimGrady1_rf = wins[line][ 6 ]
    SimGrady2_rf = wins[line][ 7 ]
    SimGrady3_rf = wins[line][ 8 ]
    SimGrady4_rf = wins[line][ 9 ]
    SimGrady5_rf = wins[line][ 10 ]
    SimGrady6_rf = wins[line][ 11 ]
    SimGrady7_rf = wins[line][ 12 ]
    SimGrady8_rf = wins[line][ 13 ]
    SimGrady9_rf = wins[line][ 14 ]
    SimGrady10_rf = wins[line][ 15 ]
    SimGrady11_rf = wins[line][ 16 ]
    SimGrady12_rf = wins[line][ 17 ]
    SimGrady13_rf = wins[line][ 18 ]
    SimGrady14_rf = wins[line][ 19 ]
    SimGrady15_rf = wins[line][ 20 ]
    SimGrady16_rf = wins[line][ 21 ]
    SimGrady17_rf = wins[line][ 22 ]
    SimGrady18_rf = wins[line][ 23 ]
    SimGrady19_rf = wins[line][ 24 ]
    SimGrady20_rf = wins[line][ 25 ]
    SimGrady21_rf = wins[line][ 26 ]
    SimGrady22_rf = wins[line][ 27 ]
    SimGrady23_rf = wins[line][ 28 ]
    SimGrady24_rf = wins[line][ 29 ]
    SimGrady25_rf = wins[line][ 30 ]
    
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
        SimGrady1_sim = AncSites[line][2] #individuals start on col 3, after chr and SNP pos
        SimGrady2_sim = AncSites[line][3] 
        SimGrady3_sim = AncSites[line][4]
        SimGrady4_sim = AncSites[line][5]
        SimGrady5_sim = AncSites[line][6]
        SimGrady6_sim = AncSites[line][7]
        SimGrady7_sim = AncSites[line][8]
        SimGrady8_sim = AncSites[line][9]
        SimGrady9_sim = AncSites[line][10]
        SimGrady10_sim = AncSites[line][11]
        SimGrady11_sim = AncSites[line][12]
        SimGrady12_sim = AncSites[line][13] 
        SimGrady13_sim = AncSites[line][14]
        SimGrady14_sim = AncSites[line][15]
        SimGrady15_sim = AncSites[line][16]
        SimGrady16_sim = AncSites[line][17]
        SimGrady17_sim = AncSites[line][18]
        SimGrady18_sim = AncSites[line][19]
        SimGrady19_sim = AncSites[line][20]
        SimGrady20_sim = AncSites[line][21]
        SimGrady21_sim = AncSites[line][22] 
        SimGrady22_sim = AncSites[line][23] 
        SimGrady23_sim = AncSites[line][24]
        SimGrady24_sim = AncSites[line][25]
        SimGrady25_sim = AncSites[line][26]
        
        if chrom == siteChrom:   #screen to match chr
            if Hapstart <= sitePos < Hapstop:
                count = count + 1   #tally up the SNPs in the window
                if SimGrady1_rf == SimGrady1_sim:
                    countTrue_1 = countTrue_1 + 1
                if SimGrady2_rf == SimGrady2_sim:
                    countTrue_2 = countTrue_2 + 1
                if SimGrady3_rf == SimGrady3_sim:
                    countTrue_3 = countTrue_3 + 1
                if SimGrady4_rf == SimGrady4_sim:
                    countTrue_4 = countTrue_4 + 1
                if SimGrady5_rf == SimGrady5_sim:
                    countTrue_5 = countTrue_5 + 1
                if SimGrady6_rf == SimGrady6_sim:
                    countTrue_6 = countTrue_6 + 1
                if SimGrady7_rf == SimGrady7_sim:
                    countTrue_7 = countTrue_7 + 1
                if SimGrady8_rf == SimGrady8_sim:
                    countTrue_8 = countTrue_8 + 1
                if SimGrady9_rf == SimGrady9_sim:
                    countTrue_9 = countTrue_9 + 1
                if SimGrady10_rf == SimGrady10_sim:
                    countTrue_10 = countTrue_10 + 1
                if SimGrady11_rf == SimGrady11_sim:
                    countTrue_11 = countTrue_11 + 1
                if SimGrady12_rf == SimGrady12_sim:
                    countTrue_12 = countTrue_12 + 1
                if SimGrady13_rf == SimGrady13_sim:
                    countTrue_13 = countTrue_13 + 1
                if SimGrady14_rf == SimGrady14_sim:
                    countTrue_14 = countTrue_14 + 1
                if SimGrady15_rf == SimGrady15_sim:
                    countTrue_15 = countTrue_15 + 1
                if SimGrady16_rf == SimGrady16_sim:
                    countTrue_16 = countTrue_16 + 1
                if SimGrady17_rf == SimGrady17_sim:
                    countTrue_17 = countTrue_17 + 1
                if SimGrady18_rf == SimGrady18_sim:
                    countTrue_18 = countTrue_18 + 1
                if SimGrady19_rf == SimGrady19_sim:
                    countTrue_19 = countTrue_19 + 1
                if SimGrady20_rf == SimGrady20_sim:
                    countTrue_20 = countTrue_20 + 1
                if SimGrady21_rf == SimGrady21_sim:
                    countTrue_21 = countTrue_21 + 1
                if SimGrady22_rf == SimGrady22_sim:
                    countTrue_22 = countTrue_22 + 1
                if SimGrady23_rf == SimGrady23_sim:
                    countTrue_23 = countTrue_23 + 1
                if SimGrady24_rf == SimGrady24_sim:
                    countTrue_24 = countTrue_24 + 1
                if SimGrady25_rf == SimGrady25_sim:
                    countTrue_25 = countTrue_25 + 1



    output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(Hapstart), str(Hapstop), str(count), str(countTrue_1),str(countTrue_2),str(countTrue_3),str(countTrue_4),str(countTrue_5),str(countTrue_6),str(countTrue_7),str(countTrue_8),str(countTrue_9),str(countTrue_10),str(countTrue_11),str(countTrue_12),str(countTrue_13),str(countTrue_14),str(countTrue_15),str(countTrue_16),str(countTrue_17),str(countTrue_18),str(countTrue_19),str(countTrue_20),str(countTrue_21),str(countTrue_22),str(countTrue_23),str(countTrue_24),str(countTrue_25))
        #outputs the chromosome, the start and stop sites of the particular window, and true-called snp counts 
        
    outFile.write(output)
outFile.close()
  