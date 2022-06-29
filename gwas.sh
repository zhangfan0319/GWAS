####speedseq比对
#!/bin/bash
#PBS -q wind
#PBS -l nodes=1:ppn=6
work=/storage1/houzc/Genome/Resequencing/naked_neck/181125_A00679_0018_BHFHYWDSXX
cd /home/houzc/WORKSPACE/fan/naked_neck/realign
gatk CreateSequenceDictionary  -R /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa
for i in `cat list`
do
speedseq align -t 12 -R "@RG\tID:${i}\tSM:${i}\tLB:lib${i}"  /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa  ${work}/${i}.R1.clean.fastq.gz  ${work}/${i}.R2.clean.fastq.gz -o ${i} 
done

#########获得 gvcf 文件
for i in `cat list`
do
gatk --java-options "-Xmx128G" HaplotypeCaller  -OVI true  -R /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa  -I ${i}_1.clean.fq.gz.bam -ERC GVCF -O ./vcffile/${i}.g.vcf.gz   -stand-call-conf 30 & ##--native-pair-hmm-threads 96
done
wait

####高深度测序的gvcf文件获取
#!/bin/bash
bamfile=/storage-01/poultrylab1/yin/database/pekin_conservation/
work=/storage-01/poultrylab1/zhangfan/temp2/pppp/
interval=/storage-01/poultrylab1/zhangfan/Genome/gatk_interval
for j in `cat list`
do
mkdir pppp/${j}
for i in {00..46}
do
gatk --java-options "-Xmx1G" HaplotypeCaller \
-OVI true -R /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa \
-I ${bamfile}/${j}_1.fq.gz.bam -ERC GVCF \
-O /storage-01/poultrylab1/zhangfan/temp2/pppp/${j}/${j}_${i}.g.vcf.gz  -stand-call-conf 30 \
-L /storage-01/poultrylab1/zhangfan/Genome/gatk_interval/${i}.interval_list 2>&1 >/storage-01/poultrylab1/zhangfan/temp2/pppp/${i}.log &
done;
wait
ls
ls /storage-01/poultrylab1/zhangfan/temp2/pppp/${j}/*vcf.gz > /storage-01/poultrylab1/zhangfan/temp2/pppp/${j}/input.list
gatk  MergeVcfs -I ./pppp/${j}/input.list -O ./pppp/${j}.vcf.gz
done



####多个文件循环  split -l 50 list -d -a 2 split_file
for j in {00..18};
do
        for i in `cat split_file$j`
        do
                gatk --java-options "-Xmx6G" HaplotypeCaller  -OVI true  -R /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa  -I ${i}D.R1.fastq.gz.bam -ERC GVCF -O ${i}.g.vcf.gz  -stand-call-conf 30 &
        done
        wait
done


##合并以及基因分型   (赵毅强学生)
import sys
import os

with open(sys.argv[1],"r") as f:
        lines=f.readlines()
        for line in lines:
                chromosome=line.strip().split(",")[0]
                length=line.strip().split(",")[1]
                for i in range(1,int(length),10000000):
                        start=int(i)
                        end=int(i+10000000-1)
                        if end >= int(length):
                                print ("nohup gatk GenomicsDBImport"+" "+"--genomicsdb-update-workspace-path  /storage-01/poultrylab1/zhangfan/blup/pop3+pop4/combinefile/"+str(chromosome)+"/"+str(start)+"-"+str(length)+" --sample-name-map sampleNameMap_pop1"+" --intervals "+str(chromosome)+":"+str(start)+"-"+str(length) +" &")  
                        else:
                                print ("nohup gatk GenomicsDBImport"+" "+"--genomicsdb-update-workspace-path  /storage-01/poultrylab1/zhangfan/blup/pop3+pop4/combinefile/"+str(chromosome)+"/"+str(start)+"-"+str(end)+" --sample-name-map sampleNameMap_pop1"+" --intervals "+str(chromosome)+":"+str(start)+"-"+str(end) +" &") 

################
import sys
import os

with open(sys.argv[1],"r") as f:
        lines=f.readlines()
        for line in lines:
                chromosome=line.strip().split(",")[0]
                length=line.strip().split(",")[1]
                for i in range(1,int(length),5000000):
                        start=int(i)
                        end=int(i+5000000-1)
                        if end >= int(length):
                                print ("/public/software/Apps/gatk-4.1.9.0/gatk GenomicsDBImport"+" "+"--genomicsdb-workspace-path /public/home/ZYQ_group_level0/wangyuzhan/EC/chinaLocalPig/"+str(chromosome)+"/"+str(start)+"-"+str(length)+" --sample-name-map sampleNameMap"+" --intervals "+str(chromosome)+":"+str(start)+"-"+str(length))
                        else:
                                print ("/public/software/Apps/gatk-4.1.9.0/gatk GenomicsDBImport"+" "+"--genomicsdb-workspace-path /public/home/ZYQ_group_level0/wangyuzhan/EC/chinaLocalPig/"+str(chromosome)+"/"+str(start)+"-"+str(end)+" --sample-name-map sampleNameMap"+" --intervals "+str(chromosome)+":"+str(start)+"-"+str(end))

#########
import sys
import os

with open(sys.argv[1],"r") as f:
        lines=f.readlines()
        for line in lines:
                chromosome=line.strip().split(",")[0]
                length=line.strip().split(",")[1]
                for i in range(1,int(length),10000000):
                        start=int(i)
                        end=int(i+10000000-1)
                        if end >= int(length):
                                print ("nohup gatk GenotypeGVCFs"+" -R /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa "+" -V gendb:///storage-01/poultrylab1/zhangfan/blup/pop3+pop4/combinefile/"+str(chromosome)+"/"+str(start)+"-"+str(length)+" --allow-old-rms-mapping-quality-annotation-data"+" -O /storage-01/poultrylab1/zhangfan/blup/pop3+pop4/genotype/"+str(chromosome)+":"+str(start)+"-"+str(length)+".vcf.gz &")
                        else:
                                print ("nohup gatk GenotypeGVCFs"+" -R /storage-01/poultrylab1/zhangfan/Genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa "+" -V gendb:///storage-01/poultrylab1/zhangfan/blup/pop3+pop4/combinefile/"+str(chromosome)+"/"+str(start)+"-"+str(end)+" --allow-old-rms-mapping-quality-annotation-data"+" -O /storage-01/poultrylab1/zhangfan/blup/pop3+pop4/genotype/"+str(chromosome)+":"+str(start)+"-"+str(end)+".vcf.gz &")

os.system(“command”)

nohup gatk  MergeVcfs -I input.list -O pop134.vcf.gz &

##限制gatk使用线程数
-XX:ConcGCThreads=1
##限制内存
--Xmx

##beagle填充
zcat vcf.gz | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > out.vcf.gz
java -Xmx64g -jar /storage-01/poultrylab1/Software/beagle.28Jun21.220.jar gt=out.vcf.gz out=10_beaagle.gt


#1_TXML
gemma-0.98.1-linux-static -bfile QC2_updated -k Centered_PrunedSNP.cXX.txt -lmm 4 -c cov5.txt -n 1 -o 1_TXML

##帆
##过滤并且转化格式
vcftools --gzvcf SNP.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --minDP 3 --minQ 30  --out SNP_bi

vcftools --gzvcf rawVariants.vcf.gz --remove-indels --max-missing 0.95 --maf 0.01 --minDP 3 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filter

vcftools --gzvcf rawVariants.SNP.vcf.gz --not-chr Z --recode --out 627_rm_contig.vcf

plink --vcf 627_rm_contig.vcf.recode.vcf --snps-only --allow-extra-chr --chr-set 40 --geno 0.05 --maf 0.01  --nonfounders --make-bed --set-missing-var-ids @:#  --out QC2

plink --bfile QC2 --indep-pairwise 50 5 0.2 --allow-extra-chr --chr-set 40 --make-founders --out QC2_Indep_SNP

plink --bfile QC2 --extract QC2_Indep_SNP.prune.in --make-bed --allow-extra-chr --chr-set 40 --out IndepSNP ##提取独立SNP位点

plink --bfile QC2 --extract QC2_Indep_SNP.prune.in --pca --allow-extra-chr  --chr-set 40 --make-founders --out QC2_IndepSNP_PCA   ##主成分分析

awk '{print $3,$4,$5,$6,$7}' QC2_IndepSNP_PCA.eigenvec > cov5.txt #把前五列作为协变量

sed -i 's/^/1 /g' cov5.txt 

##替换GWAS_GEMMA_Relatedness.mds文件内容为表型
##替换QC2.fam文件内容为表型

gemma -bfile IndepSNP -gk 1 -o Centered_PrunedSNP -c cov5.txt
##gemma多核？？？
export OPENBLAS_NUM_THREADS=4
##output产生Centered_PrunedSNP.cXX.txt文件
./gemma -bfile [prefix] -lm [num] -o [prefix]   #lm
gemma -bfile QC2 -k Centered_PrunedSNP.cXX.txt -lmm 1 -c cov5.txt -n 1 -o trait  #lmm

cat 23_map,txt| awk -F" " '$4< 1/606491{print $0}' > candidata/23_wei

##gwas结果绘图
cat FI_result.assoc.txt |awk '{print $2,$1,$3,$15}' > FI_map.txt
library(rMVP)
ff <- read.table("23_days_wei_map.txt",header=T,sep=" ")
MVP.Report(ff,threshold = c(1/606491,0.05/606491),plot.type=c("m","q"),file.type = "jpg",dpi=300,outpath="/storage-01/poultrylab1/zhangfan/gwas/gemma/norm/output/23_wei")
MVP.Report(TBW,threshold = 5.036e-09,plot.type=c("m","q"),file.type = "jpg",dpi=300,outpath="/storage-01/poultrylab1/zhangfan/gwas/mtag/feng/TBW")
5.036e-09  ##mtag_峰哥的阈值

##########fastlmm#############
cd ~/zhangfan/work/CLN8
gcta64 --bfile ../../../../QC2 --make-grm --sparse-cutoff 0.05 --thread-num 150 --out sp_grm --autosome-num 40
gcta64 --bfile ../../../../QC2 --grm-sparse sp_grm --fastGWA-mlm --pheno nsfw.txt --qcovar cov.txt --covar sby.txt --thread-num 100 --out sfw_fastgmm --autosome-num 40

##########MAGMA的使用###############################################################################################
cat yeye_gene.gff |awk -v FS="\t" -v OFS="\t"  '{print $3,$1,$4,$5,$7,$9}' > gene.loc 

magma --annotate nonhuman --snp-loc QC2.bim --gene-loc gene.loc --out anno ##注释
#   --annotate window=5,5’ would set a 5kb upstream and 5kb downstream window
#example: nohup magma --bfile ../fat/QC2 --gene-annot anno.genes.annot --out ./result/23_body_weight &   
./magma --gene-model multi=all --bfile QC2 --gene-annot anno.genes.annot --covar file=covfile --pheno file=ntrait use=1 --out 42_wei  ##基因水平上的关联检验
#example: nohup  magma --bfile QC2 --gene-annot anno.genes.annot --out ./result/23_body_weight &

##cd /storage-01/poultrylab1/zhangfan/blup/pop1+pop3+pop4/analysis1/magma 
##根据p值估算基因的p值
./magma --bfile ../QC2 --pval sfw_plot.txt N=1800 --gene-annot anno.genes.annot  --out fffff


########### MTAG 使用 ###############################################
##教程 https://github.com/omeed-maghzian/mtag/wiki/Tutorial-1:-The-Basics

#去表头
sed -i 's/ /\t/g' mtag_skin_fat_wei.txt
##estimate the implicit N as N=1/[SE^2*2*p*(1-p)] where p is the allele frequency
sed -i '1d' 23_wei
cat skin_fat_wei.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_skin_fat_wei.txt
cat skin_fat_rate.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_skin_fat_rate.txt
cat shearing_force.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_shearing_force.txt
cat fatty_acid.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_fatty_acid.txt
cat abdominal_fat_rate.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_abdominal_fat_rate.txt
cat abdominal_fat_wei.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_abdominal_fat_wei.txt
cat 23_wei.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_23_wei.txt
cat 42_wei.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_42_wei.txt
cat breast_wei.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_breast_wei.txt
cat breast_rate.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_breast_rate.txt
cat breast_thickness.assoc.txt|awk '{print $2,$1,$3,$6,$5,$7,$8,sqrt($9),$14,1/($9^2*2*$7*(1-$7))}' > ../mtag_breast_thickness.txt
##snp名称更换格式
sed -r 's@(\w):@rs\1@' mtag_42_wei.txt > mtag_42_wei_t.txt
##表头  snpid chr bpos a1 a2 freq beta se pval n

python2.7 ./mtag-master/mtag.py --sumstats mtag_23_wei.txt,mtag_42_wei.txt --cores 40 \
--out ./result --n_min 0.0 --stream_stdout --use_beta_se \
--snp_name rs --bpos_name ps --a1_name allele0 --a2_name allele1 --eaf_name af --p_name p_lrt --n_name 627

python mtag/mtag.py --sumstats mtag_t_42_wei.txt,mtag_t_skin_fat_wei.txt,mtag_t_abdominal_fat_wei.txt,mtag_t_breast_wei.txt  --cores 80 --out ./result_mtag/wei_42_skin_abd_breast --n_min 0.0 --stream_stdout --use_beta_se --force

cat 23_map,txt| awk -F" " '$12< 0.05/606491{print $0}' > candidata/23_wei
############GCTA分析#########################
关联分析就是将突变与表型联系起来，如果表型是分类的，最简单的方法就是卡方检验。plink - -assoc 原理就是这样，但不能添加协变量，plink - -logistic可以做到这一点。如果表型是连续的数量表型，例如身高等，可以用 plink --linear ,也可以添加协助变量。
先将协变量信息写入到文件 SNP.covar，前两列是FID and IID，后面的列是协变量，主要有年龄、身高、体重等。表型文件 inhi.pheno1 前两列也是FID and IID，第三列是表型。

plink 关联分析

plink --bfile SNP --covar SNP.covar --logistic --hide-covar --out  bleeding 
plink --bfile SNP --covar SNP.covar --pheno inhi.pheno1 --linear --hide-covar --out  inhi1
plink --bfile SNP --covar SNP.covar --pheno inhi.pheno2 --logistic --hide-covar --out  inhi2

https://www.jianshu.com/p/06945cdf034e    
http://www.biotrainee.com/thread-1497-1-1.html
https://cnsgenomics.com/software/gcta/#MLMA  

###################GWAS后续分析################
##计算亲缘关系矩阵
gcta64 --bfile test --make-grm-gz --make-grm-alg 1 --out kinship
##LD绘图参考
../../bin/LDBlockShow -InVCF ../Example1/Test.vcf.gz -OutPut out -InGWAS gwas.pvalue   -InGFF    In.gff   -Region    chr11:24100000:24200000 
../../bin/ShowLDShow -InPreFix out -OutPut out.svg -InGWAS gwas.pvalue -Cutline 7 -InGFF In.gff -crGene yellow:lightblue:pink:orange -showNum -OutPng    -SpeSNPName Spe.snp    -ShowGWASSpeSNP 
##LD绘图实战
cat ../sfr.assoc.txt |grep '^[3:]'|awk '$3 >22493400 && $3<22493794 {print $1,$3,$12}' >LD_gwas_sfr.txt 
/storage-01/poultrylab1/Software/LDBlockShow/bin/LDBlockShow -InVCF ../../../filter.recode.vcf.gz -OutPut file -Region 11:7361442-7363442  -OutPdf  -showNum 
/storage-01/poultrylab1/Software/LDBlockShow/bin/LDBlockShow -InVCF ../../../filter.recode.vcf.gz -OutPut file -Region 11:7361442-7363442  -OutPdf -InGWAS LD_gwas_sfw.txt -Cutline 7.06
##GCTA工具计算复杂性状/特征（Complex Trait）的遗传相关性（genetic correlation）##
gcta64 --bfile test --autosome-num 40 --make-grm --out test --thread-num 100 #生成grm格式文件，方便后面的遗传相关性分析
##准备性状/表型文件，后缀为.txt格式，不需要表头，第一列为family ID, 第二列为individual ID 第三列和第四列为 phenotypes ，类似于PLINK的表型文件格式
gcta64 --reml-bivar --reml-bivar-nocove --grm test --pheno pheno.txt --reml-bivar-lrt-rg 0 --out test --autosome-num 40
#########计算性状1和2、3、4的遗传相关
for i in $(seq 1 4); do
echo $i
nohup gcta64 --reml-bivar 2 $i --reml-bivar-nocove --grm test --pheno pheno.txt  --reml-bivar-lrt-rg 0 --out trait2_${i%} &
done


####计算单个SNP的遗传力
gcta64 --grm test --pheno test.phen --mpheno 1 --reml --qcovar test_10PCs.txt --out test --thread-num 10 --autosome-num 40
####脚本
for i in {1..370}
do
sed -n ${i}p list > ${i}   ###list为要分析SNP的位置
for j in `cat ${i}`
do
plink --bfile ../QC2 --extract ${i} --make-bed --out ${j} --allow-extra-chr --chr-set 40
gcta64 --bfile ${j} --autosome-num 40 --maf 0.01 --make-grm --out ${j} --thread-num 100
gcta64 --grm ${j} --pheno nPheno.txt --reml --mpheno 1 --out ${j} --thread-num 100
rm ${i}
done
done

####ldsc计算遗传相关  /storage-01/poultrylab1/zhangfan/blup/pop1+pop3+pop4/analysis1/gcta/genetic_cor
conda activate ldsc
for i in sfr fac afw bw42 bmw sfw sft bmr afr sfr bmt sf; do cat ../../output/nsbyc5/${i}.assoc.txt |awk '{print $2,$1,$3,$5,$6,$8,$12}' > ./${i}.summary.txt & done; wait
for i in sfr fac afw bw42 bmw sfw sft bmr afr sfr bmt sf; do sed -i 's/ps allele1 allele0 beta p_wald/bp a1 a2 beta pval/g' ${i}.summary.txt  & done ;wait
for i in sfr fac afw bw42 bmw sfw sft bmr afr sfr bmt sf; do munge_sumstats.py --sumstats ${i}.summary.txt --N 17115 --out ${i}.scz & done ;wait
for q in $(seq 1 40); do plink --bfile ../../QC2 --chr $q --make-bed --allow-extra-chr --chr-set 40 --out chr$q  & done;wait
mkdir chr
### 分染色体；计算LD
for q in $(seq 1 40); do ldsc.py --bfile chr$q --l2 --ld-wind-cm 5 --yes-really --out chr/$q & done ;wait
##回归截距和遗传度
for i in sfr fac afw bw42 bmw sfw sft bmr afr sfr bmt sf; do ldsc.py --h2 ${i}.scz.sumstats.gz --ref-ld-chr chr/ --w-ld-chr chr/ --out ${i}_h2 & done; wait
for i in sfr fac afw bw42 bmw sfw sft bmr afr sfr bmt sf; do ldsc.py --rg ${i%}.scz.sumstats.gz,bw42.scz.sumstats.gz,sfw.scz.sumstats.gz,sfr.scz.sumstats.gz,sft.scz.sumstats.gz,afw.scz.sumstats.gz,afr.scz.sumstats.gz,bmw.scz.sumstats.gz,bmr.scz.sumstats.gz,bmt.scz.sumstats.gz,sf.scz.sumstats.gz,fac.scz.sumstats.gz --ref-ld-chr chr/ --w-ld-chr chr/ --out --out trait${i%}linearcorr & done; wait




import sys
import os

with open(sys.argv[1],"r") as f:
        lines=f.readlines()
        for line in lines:
                chromosome=line.strip().split(",")[0]
                length=line.strip().split(",")[1]
                for i in range(1,int(length),10000000):
                        start=int(i)
                        end=int(i+10000000-1)
                        if end >= int(length):
                                print ("nohup gatk GenomicsDBImport"+" "+"--genomicsdb-update-workspace-path  /storage-04/temp3/VcfDataBase/combinefile/"+str(chromosome)+"/"+str(start)+"-"+str(length)+" --sample-name-map sampleNameMap_pop2019"+" --intervals "+str(chromosome)+":"+str(start)+"-"+str(length) +" &")
                        else:
                                print ("nohup gatk GenomicsDBImport"+" "+"--genomicsdb-update-workspace-path  /storage-04/temp3/VcfDataBase/combinefile/"+str(chromosome)+"/"+str(start)+"-"+str(end)+" --sample-name-map   sampleNameMap_pop2019"+" --intervals "+str(chromosome)+":"+str(start)+"-"+str(end) +" &")



