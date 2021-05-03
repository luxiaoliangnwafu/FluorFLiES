#!/bin/csh

 g77 -funroll-all-loops -O3 scatvec.f scatangletbl.f  transform.f  frndi.f frnd.f  imath.f math.f igtbl.f
