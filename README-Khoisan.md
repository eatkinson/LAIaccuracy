# LAIaccuracy
Calculations of accuracy comparing simulations to RFmix runs


Ancestry call accuracy from rfmix v2 was compared to known simulated individuals, using Amy Williams's script to simulate KhoeSan individuals. Then ran RFmix Version 2 on the simulated files. The included python scripts tabulate how often sites were truly called as the ancestry they were over an individual globally, and how often sites for a particular ancestry were correctly called.

We estimated the accuracy of RFmix_v2 local ancestry calls on SNP array data for 3 admixed population models corresponding to the Zulu, Nama, and East African Bantu-speaking groups. The goal of these simulations was to assess how well LAI works for 1) ancient gene flow events, 2) between minimally-divergent ancestral groups, and 3) for datasets where we have few proxy individuals for a given component ancestry. These scenarios are not yet well-characterized in published LAI papers, and a lack of information is especially apparent for African populations. In general, we found that the largest issue plaguing LAI is access to large, appropriate proxy ancestral populations. LAI suffered most greatly from failure to sufficiently capture representative haplotypes from contributing groups, as most prominently evidenced for the KhoeSan populations, for whom there are few existing reference genomes available. 
