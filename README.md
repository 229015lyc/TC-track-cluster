**Two main coodes:**

# 1. runCluster.py

   This code is modified from CCToolbox(http://www.datalab.uci.edu/resources/CCT/), which is a Matlab code.
   
   Input data include TC track data (latitude, lontitude) and the index number of TC.
   
   You can set the number of clusters (variable "c" in this code) to try how many clusters are acceptable. 
# 2. calculation_VI.py
   
   This code is an example to show you how I calculate VI. Before running it, PI and wind shear should be calculated and saved as files. It is recommended that computing PI via pyPI v1.3 (Gilford, D. M.,2020). Many variables (vmax,diseq,sst,q,t) follow the output format of it.  
