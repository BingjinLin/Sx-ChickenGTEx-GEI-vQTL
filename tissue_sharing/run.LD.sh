nohup: ignoring input
+ geno_file=/juice_data/zhh_independent/ZHH_update/Genotype/Chicken.SNP_Autosomal.recode
+ snp_file=/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/top_pairs/multi.top_vqtl.list
+ Otfile=/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl
+ echo 'Extract top vqtl'
Extract top vqtl
+ plink --bfile /juice_data/zhh_independent/ZHH_update/Genotype/Chicken.SNP_Autosomal.recode --chr-set 39 --keep-allele-order --make-bed --extract /juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/top_pairs/multi.top_vqtl.list --out /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl
PLINK v1.90b7.4 64-bit (18 Aug 2024)           www.cog-genomics.org/plink/1.9/
(C) 2005-2024 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.log.
Options in effect:
  --bfile /juice_data/zhh_independent/ZHH_update/Genotype/Chicken.SNP_Autosomal.recode
  --chr-set 39
  --extract /juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/top_pairs/multi.top_vqtl.list
  --keep-allele-order
  --make-bed
  --out /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl

128835 MB RAM detected; reserving 64417 MB for main workspace.
17489273 variants loaded from .bim file.
280 samples (0 males, 0 females, 280 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.nosex .
--extract: 67974 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 280 founders and 0 nonfounders present.
Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
Total genotyping rate is exactly 1.
67974 variants and 280 samples pass filters and QC.
Note: No phenotypes present.
--make-bed to
/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.bed +
/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.bim +
/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
+ echo 'Calculate LD top vqtl'
Calculate LD top vqtl
+ plink --bfile /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl --chr-set 39 --r2 --out /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl
PLINK v1.90b7.4 64-bit (18 Aug 2024)           www.cog-genomics.org/plink/1.9/
(C) 2005-2024 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.log.
Options in effect:
  --bfile /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl
  --chr-set 39
  --out /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl
  --r2

128835 MB RAM detected; reserving 64417 MB for main workspace.
67974 variants loaded from .bim file.
280 samples (0 males, 0 females, 280 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.nosex .
Using up to 15 threads (change this with --threads).
Before main variant filters, 280 founders and 0 nonfounders present.
Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
Total genotyping rate is exactly 1.
67974 variants and 280 samples pass filters and QC.
Note: No phenotypes present.
Running --r2 with the following filters:
  --ld-window: 10
  --ld-window-kb: 1000
  --ld-window-r2: 0.2
--r2 to /juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl.ld ...
0% [processing]writing]             done.
