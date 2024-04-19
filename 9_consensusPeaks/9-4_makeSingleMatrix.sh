#Ellen Ketter, 17 July 2023

awk '{print $1$6$2$6$3,$11}' 5F1_noDup_counts.txt > 5F1_counts.txt
awk '{print $11}' 5F2_noDup_counts.txt > 5F2_counts.txt
awk '{print $11}' 5F3_noDup_counts.txt > 5F3_counts.txt
paste 5F1_counts.txt 5F2_counts.txt 5F3_counts.txt > 5F_counts.txt 
sed  -i '1i peak\t5F-1\t5F-2\t5F-3' 5F_counts.txt

awk '{print $1$6$2$6$3,$11}' BG1_noDup_counts.txt > BG1_counts.txt
awk '{print $11}' BG3_noDup_counts.txt > BG3_counts.txt
paste BG1_counts.txt BG3_counts.txt > BG_counts.txt 
sed  -i '1i peak\tBG-1\tBG-3' BG_counts.txt

awk '{print $1$6$2$6$3,$11}' Fen1_noDup_counts.txt > Fen1_counts.txt
awk '{print $11}' Fen2_noDup_counts.txt > Fen2_counts.txt
awk '{print $11}' Fen3_noDup_counts.txt > Fen3_counts.txt
paste Fen1_counts.txt Fen2_counts.txt Fen3_counts.txt > Fen_counts.txt 
sed  -i '1i peak\tFen-1\tFen-2\tFen-3' Fen_counts.txt

awk '{print $1$6$2$6$3,$11}' Fluni1_noDup_counts.txt > Fluni1_counts.txt
awk '{print $11}' Fluni2_noDup_counts.txt > Fluni2_counts.txt
awk '{print $11}' Fluni3_noDup_counts.txt > Fluni3_counts.txt
paste Fluni1_counts.txt Fluni2_counts.txt Fluni3_counts.txt > Fluni_counts.txt 
sed  -i '1i peak\tFluni-1\tFluni-2\tFluni-3' Fluni_counts.txt

awk '{print $1$6$2$6$3,$11}' HC1_noDup_counts.txt > HC1_counts.txt
awk '{print $11}' HC2_noDup_counts.txt > HC2_counts.txt
awk '{print $11}' HC3_noDup_counts.txt > HC3_counts.txt
paste HC1_counts.txt HC2_counts.txt HC3_counts.txt > HC_counts.txt 
sed  -i '1i peak\tHC-1\tHC-2\tHC-3' HC_counts.txt

awk '{print $1$6$2$6$3,$11}' HQ1_noDup_counts.txt > HQ1_counts.txt
awk '{print $11}' HQ2_noDup_counts.txt > HQ2_counts.txt
paste HQ1_counts.txt HQ2_counts.txt > HQ_counts.txt 
sed  -i '1i peak\tHQ-1\tHQ-2' HQ_counts.txt

awk '{print $1$6$2$6$3,$11}' Myr1_noDup_counts.txt > Myr1_counts.txt
awk '{print $11}' Myr2_noDup_counts.txt > Myr2_counts.txt
awk '{print $11}' Myr3_noDup_counts.txt > Myr3_counts.txt
paste Myr1_counts.txt Myr2_counts.txt Myr3_counts.txt > Myr_counts.txt 
sed  -i '1i peak\tMyr-1\tMyr-2\tMyr-3' Myr_counts.txt

awk '{print $1$6$2$6$3,$11}' Nerol1_noDup_counts.txt > Nerol1_counts.txt
awk '{print $11}' Nerol2_noDup_counts.txt > Nerol2_counts.txt
awk '{print $11}' Nerol3_noDup_counts.txt > Nerol3_counts.txt
paste Nerol1_counts.txt Nerol2_counts.txt Nerol3_counts.txt > Nerol_counts.txt 
sed  -i '1i peak\tNerol-1\tNerol-2\tNerol-3' Nerol_counts.txt

awk '{print $1$6$2$6$3,$11}' PBS1_noDup_counts.txt > PBS1_counts.txt
awk '{print $11}' PBS2_noDup_counts.txt > PBS2_counts.txt
awk '{print $11}' PBS3_noDup_counts.txt > PBS3_counts.txt
awk '{print $11}' PBSa1_noDup_counts.txt > PBS4_counts.txt
awk '{print $11}' PBSa2_noDup_counts.txt > PBS5_counts.txt
awk '{print $11}' PBSa3_noDup_counts.txt > PBS6_counts.txt
paste PBS1_counts.txt PBS2_counts.txt PBS3_counts.txt PBS4_counts.txt PBS5_counts.txt PBS6_counts.txt> PBS_counts.txt 
sed  -i '1i peak\tPBS-1\tPBS-2\tPBS-3\tPBS-4\tPBS-5\tPBS-6' PBS_counts.txt

join 5F_counts.txt Fen_counts.txt > combo1.txt
join Fluni_counts.txt HC_counts.txt > combo2.txt
join HQ_counts.txt Myr_counts.txt > combo3.txt
join Nerol_counts.txt PBS_counts.txt > combo4.txt

join combo1.txt combo2.txt > combo1-2.txt
join combo3.txt combo4.txt > combo3-4.txt

join combo1-2.txt combo3-4.txt > countMatrix.txt
head countMatrix.txt
