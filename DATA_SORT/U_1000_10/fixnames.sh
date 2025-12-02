#!/bin/bash
cd /glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/U_1000_10/L32
for iyear in `seq 1971 2020`; do
   mv U_1000_10_L83_init11_$iyear'.nc' U_1000_10_L32_init11_$iyear'.nc'
done
