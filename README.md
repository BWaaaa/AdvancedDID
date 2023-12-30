# AdvancedDID
Replication and Extension on DID methods
Date: Dec 30 2023
Author: Brian Yisong Wang

This is a STATA replication AND extension package for part of the "The air quality and well-being effects of low emission zones" by Sarmeinto et.al.(2023).
In particular, this package focuses around Table 6 and Table 8, but 
1. Extends the analysis to include other Advanced DID methods
2. Codes in STATA while the original paper uses R

Here is a short description of the files

:: pol.dta    :    a balanced panel data of pollutants, used for 2-by-2 DID and Bacon Decomposition
:: allpol.dta :    the raw data (unbalanced) of pollutants, used for main regression results

:: Advanced_DIDweapon.do       :     Reader-friendly version of the codes
:: Advanced_DIDweapon_noc.do   :     Concisely written codes

:: Graphs/   :    directory for 17 replication AND extesion figures for Stata-lab purpose.

Here is a quick view:
![fivec](https://github.com/BWaaaa/AdvancedDID/blob/1333f5b57966537a8092154887cd067a393af60f/Graphs/figure17_fivecombine.png)
