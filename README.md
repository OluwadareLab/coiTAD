# coiTAD: Topologically Associating Domains Detection based on Clustering of Circular Influence Features from Hi-C data
To use coiTAD:
1. Change filepath variables to match your system.
2. Ensure Chr_Data refers to the dataset that you are working on in both the Feature Generation and Main code files
3. Run the Program

Finding Results:
1. coiTAD will print the best radius selected to the console. Remember this radius.
2. Open the Data_Results file and then TADs. Locate text domain file associated with your best radius.
(This file will follow the structure HDBSCAN_Radius#_domain.txt)
3. This file contains the TAD domains. It can be converted to a .bed file and used for comparison. 
