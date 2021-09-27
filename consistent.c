/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//////////////////////////

//mainly check arguments provided are consistent
//also set a few variables
//note that in the code, modes are always in numerical order, but when printing, they are in readargs order

///////////////////////////

//data input - will set filenames for dtypes 11-15, plus genskip, genheaders and genprobs

if(dtype==11)	//set values corresponding to --sp
{
if(strcmp(ubimfile,"blank")!=0)
{printf("Error, you can not use \"--bim\" with \"--sp\"\n\n");exit(1);}
if(strcmp(ufamfile,"blank")!=0)
{printf("Error, you can not use \"--fam\" or \"--sample\" with \"--sp\"\n\n");exit(1);}

sprintf(udatafile,"%s.sp",unamefile);
if(strcmp(ubimfile,"blank")==0){sprintf(ubimfile,"%s.bim",unamefile);}
if(strcmp(ufamfile,"blank")==0){sprintf(ufamfile,"%s.fam",unamefile);}
famhead=0;
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=0;}
if(genprobs==-9999){genprobs=1;}
dtype=5;
}

if(dtype==12)	//set values corresponding to --sp-gz
{
if(strcmp(ubimfile,"blank")!=0)
{printf("Error, you can not use \"--bim\" with \"--sp-gz\"\n\n");exit(1);}
if(strcmp(ufamfile,"blank")!=0)
{printf("Error, you can not use \"--fam\" or \"--sample\" with \"--sp-gz\"\n\n");exit(1);}

sprintf(udatafile,"%s.sp.gz",unamefile);
if(strcmp(ubimfile,"blank")==0){sprintf(ubimfile,"%s.bim",unamefile);}
if(strcmp(ufamfile,"blank")==0){sprintf(ufamfile,"%s.fam",unamefile);}
famhead=0;
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=0;}
if(genprobs==-9999){genprobs=1;}
dtype=5;
}

if(dtype==13)	//set values corresponding to --beagle-dose
{
if(genskip==-9999){genskip=1;}
if(genheaders==-9999){genheaders=3;}
if(genprobs==-9999){genprobs=1;}
dtype=5;
}

if(dtype==14)	//set values corresponding to --beagle-probs
{
if(genskip==-9999){genskip=1;}
if(genheaders==-9999){genheaders=3;}
if(genprobs==-9999){genprobs=3;}
dtype=5;
}

if(dtype==15)	//set values corresponding to --haps
{
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=5;}
if(genprobs==-9999){genprobs=0;}
dtype=5;
}

////////

if(strcmp(datalist,"blank")!=0&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=190)
{printf("Error, you can only use \"--mbfile\", \"--msp\", \"--msped\" or \"--mspeed\" when making data or with \"--calc-sim-data\"\n\n");exit(1);}

if((mode==151||mode==152||mode==153)&&(dtype==3||dtype==5))
{printf("Error, you can not use \"--sped\" or \"--gen\" with \"--ridge\", \"--bolt\" or \"--bayesr\"; you should first remake the data using either \"--make-bed\" or \"--make-speed\"\n\n");exit(1);}

if(mode==185&&genprobs<2)
{printf("Error, you can only use \"--make-gen\" when providing genotype probabilities (\"--gen\" or \"--mgen\"); consider instead using \"--make-bed\" or \"--make-speed\"\n\n");exit(1);}

////////

if(strcmp(ubimfile,"blank")!=0&&(dtype!=5||strcmp(datalist,"blank")!=0))
{printf("Error, you can only use \"--bim\" with \"--gen\"\n\n");exit(1);}

if(strcmp(ufamfile,"blank")!=0&&(dtype!=5||strcmp(datalist,"blank")!=0))
{printf("Error, you can only use \"--fam\" or \"--sample\" with \"--gen\"\n\n");exit(1);}

if(oxchr!=-9999&&dtype!=5)
{printf("Error, you can only use \"--oxford-single-chr\" with \"--gen\" or \"--mgen\"\n\n");exit(1);}

if((genskip!=-9999||genheaders!=-9999||genprobs!=-9999)&&dtype!=5)
{printf("Error, you can only use \"--gen-skip\", \"--gen-headers\" or \"--gen-probs\" with \"--gen\" or \"--mgen\"\n\n");exit(1);}

if(nonsnp!=-9999&&dtype==-9999)
{printf("Error, you can only use \"--SNP-data\" when providing data\n\n");exit(1);}

if(nonsnp==1&&dtype==1)
{printf("Error, you can not use \"--SNP-data NO\" with \"--bfile\" or \"--mbfile\"; binary PLINK format handles only SNP data\n\n");exit(1);}

if(nonsnp==1&&genprobs==0)
{printf("Error, you can not use \"--SNP-data NO\" when providing haplotypes\n\n");exit(1);}

if(nonsnp==1&&genprobs>1)
{printf("Error, you can not use \"--SNP-data NO\" when providing genotype probabilities\n\n");exit(1);}

//need to set genprob and nonsnp now (will also set genskip and genheaders)

if(genskip==-9999)
{
if(dtype==5){genskip=0;}
}
if(genheaders==-9999)
{
if(dtype==5){genheaders=5;}
}

if(genprobs==-9999)
{
if(dtype==1||dtype==3||dtype==4){genprobs=1;}
if(dtype==5){genprobs=3;}
}

if(nonsnp==-9999&&dtype!=-9999)
{
nonsnp=0;
if(dtype==3||dtype==4||(dtype==5&&genprobs==1))
{printf("Warning, predictors are assumed to be SNP allele counts (i.e., take values within [0,2]); if this is not the case, you should use \"--SNP-data NO\"\n\n");}
}

if(dtype==5&&strcmp(datalist,"blank")==0)	//some checks for --gen (do --mgen checks when reading datalist)
{
if(genheaders<4&&strcmp(ubimfile,"blank")==0)
{printf("Error, you must provide SNP details using \"--bim\"\n\n");exit(1);}

if(oxchr!=-9999&&strcmp(ubimfile,"blank")!=0)
{printf("Error, you can not use \"--oxford-single-chr\" if using \"--bim\"\n\n");exit(1);}

if(strcmp(ubimfile,"blank")==0&&oxchr==-9999)
{printf("Error, when using \"--gen\" you must also use \"--bim\" or \"--oxford-single-chr\"\n\n");exit(1);}

if(strcmp(ufamfile,"blank")==0)
{printf("Error, when using \"--gen\" you must also use \"--fam\" or \"--sample\"\n\n");exit(1);}
}

///////////////////////////

//data filtering

if((strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0)&&mode==113)
{printf("Error, you can not use \"--keep\" or \"--remove\" with \"--join-kins\"\n\n");exit(1);}

if((strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0)&&mode==162)
{printf("Error, you can not use \"--keep\" or \"--remove\" with \"--calc-pca-loads\" (the samples will be read from the .vect file)\n\n");exit(1);}

if(num_subs!=-9999||strcmp(subpref,"blank")!=0)
{
if(mode!=102&&mode!=104&&mode!=123&&mode!=124&&mode!=171)
{printf("Error, you can only use sample subsets with \"--calc-weights\", \"--calc-weights-all\", \"--he\", \"--pcgc\" or \"--calc-stats\"\n\n");exit(1);}

if(strcmp(subpref,"blank")==0)
{printf("Error, when using \"--subset-number\", you must also use \"--subset-prefix\"\n\n");exit(1);}

if(num_subs==-9999)
{printf("Error, when using \"--subset-prefix\", you must also use \"--subset-number\"\n\n");exit(1);}

if(num_subs==1)
{
if(mode==102||mode==104)
{printf("Warning, using \"--subset-number 1\" and \"--subset-prefix %s\" is equivalent to (simply) using \"--keep %s1\"\n\n", subpref, subpref);}
if(mode==123||mode==124)
{printf("Error, to estimate heritability across and within cohorts requires at least two subsets (not 1)\n\n");exit(1);}
}
}	//end of using subsets

////////

if(mode==108&&(onechr!=-9999||strcmp(onesnp,"blank")!=0))
{printf("Error, you can not use \"--chr\" or \"--snp\" with \"--find-tags\"; you should instead use \"--extract\" and/or \"--exclude\" to specify which predictors can be used as tags\n\n");exit(1);}

if(mode==117&&(onechr!=-9999||strcmp(onesnp,"blank")!=0))
{printf("Error, you can not use \"--chr\" or \"--snp\" with \"--sub-grm\"; to subtract predictors from a kinship matrix use either \"--extract\" or \"--exclude\"\n\n");exit(1);}

if(mode==117&&strcmp(bpredfile,"blank")!=0&&strcmp(cpredfile,"blank")!=0)
{printf("Error, you can not use both \"--extract\" and \"--exclude\" with \"--sub-grm\"\n\n");exit(1);}

if(mode==122&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--calc-blups\" (the predictors will be read from the grm.details file)\n\n");exit(1);}

if((mode==146||mode==147)&&(onechr!=-9999||strcmp(onesnp,"blank")!=0))
{printf("Error, you can not use \"--chr\" or \"--snp\" with \"--sum-hers\" or \"--sum-cors\")\n\n");exit(1);}

if(mode==149&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--calc-exps\" (expectations will be calculated for all predictors in the tagging file)\n\n");exit(1);}

if(mode==162&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--calc-pca-loads\" (the predictors will be read from the grm.details file)\n\n");exit(1);}

if(strcmp(onesnp,"blank")!=0&&(strcmp(bpredfile,"blank")!=0||strcmp(cpredfile,"blank")!=0||onechr!=-9999))
{printf("Error, when using \"--snp\", you can not use \"--extract\", \"--exclude\" or \"--chr\"\n\n");exit(1);}

if(strcmp(onesnp,"blank")!=0&&mode!=131&&mode!=132&&mode!=171&&mode!=172&&mode!=173&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=190)
{printf("Warning, I think it only makes sense to use \"--snp\" with \"--linear\", \"--logistic\", \"--calc-stats\", \"--calc-scores\", \"--make-phenos\", when making data or with \"--calc-sim-data\"\nBut if you can think of another use, please write in :)\n\n");}

if(mode==174&&(minmaf!=-9999||maxmaf!=-9999))
{printf("Error, when using \"--make-snps\", you should use \"--maf-low\" and \"--maf-high\" (rather than \"--min-maf\" and \"--max-maf\")\n\n");exit(1);}

if((minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999)&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--min-maf\", \"--max-maf\", \"--min-var\" or \"--min-obs\" when making or condensing data; either first remake the data, or use \"--extract\" or \"--exclude\" to specify which predictors to analyse (you can use \"--calc-stats\" to calculate predictor allele frequencies and call-rates)\n\n");exit(1);}

if((minmaf!=-9999||maxmaf!=-9999)&&nonsnp==1)
{printf("Error, you can not use \"--min-maf\" or \"--max-maf\" with \"--SNP-data NO\"\n\n");exit(1);}

if(minmaf!=-9999&&maxmaf!=-9999&&minmaf>=maxmaf)
{printf("Error, \"--min-maf\" (%.6f) must be lower than \"--max-maf\" (%.6f)\n\n", minmaf, maxmaf);exit(1);}

if(mininfo!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--min-info\" when making data (and when providing genotype probabilities)\n\n");exit(1);}

if(mininfo!=-9999&&(genprobs==0||genprobs==1))
{printf("Error, you can only use \"--min-info\" when providing genotype probabilities (\"--gen\" or \"--mgen\"), as these are required to calculate the info score for each SNP\n\n");exit(1);}

///////////////////////////

//data scaling and coding

if(strcmp(centresfile,"blank")!=0&&mode!=112&&mode!=114&&mode!=126&&mode!=127)
{printf("Error, you can only use \"--predictor-means\" with \"--calc-kins\", \"--calc-kins-direct\", \"--fast-he\" or \"--fast-pcgc\"\n\n");exit(1);}

if((mode==101||mode==102||mode==103||mode==104)&&strcmp(weightsfile,"blank")!=0)
{printf("Sorry, it is no longer possible to calculate weightings twice; a better solution is to first thin predictors to remove duplicates (which LDAK now does by default when cutting predictors into sections)\n\n");exit(1);}

if((mode==101||mode==102||mode==103||mode==104)&&(power!=-9999||hwestand!=-9999))
{printf("Error, you can not use \"--power\" or \"--hwe-stand\" when calculating weightings n\n");exit(1);}

if((mode==117||mode==122||mode==162)&&(strcmp(weightsfile,"blank")!=0||power!=-9999||hwestand!=-9999))
{printf("Error, you can not use \"--weights\", \"--power\" or \"--hwe-stand\" with \"--calc-blups\", \"--sub-grm\" or \"--calc-pca-loads\" (the predictor scalings will be read from the grm.details file)\n\n");exit(1);}

if(mode==172&&strcmp(weightsfile,"blank")!=0)
{printf("Error, you can not use \"--weights\" with \"--calc-scores\" (if you wish to include weightings, you should incorporate them in the predictor effect sizes)\n\n");exit(1);}

if(strcmp(weightsfile,"blank")!=0&&ignoreweights==1){printf("Error, it is contradictory to use both \"--weights\" and \"--ignore-weights YES\"\n\n");exit(1);}

if(hwestand==1&&nonsnp==1)
{printf("Error, you can not use \"--hwe-stand YES\" with \"--SNP-data NO\"\n\n");exit(1);}

if(strcmp(pvafile,"blank")!=0&&mode!=106&&mode!=107&&mode!=136&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--pvalues\" with \"--cut-genes\", \"--thin\", \"--thin-tops\"  or when condensing data\n\n");exit(1);}

////////

if(encoding!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--encoding\" when making data\n\n");exit(1);}

if(encoding!=-9999&&dtype!=1)
{printf("Error, you can only use \"--encoding\" with \"--bfile\"; if your data are in a different format, you should first remake using \"--make-bed\"\n\n");exit(1);}

if((threshold!=-9999||minprob!=-9999)&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--threshold\" or \"--min-prob\" when making data\n\n");exit(1);}

if(threshold!=-9999&&nonsnp==1)
{printf("Error, you can not use \"--threshold\" with \"--SNP-data NO\"\n\n");exit(1);}

if(threshold!=-9999&&dtype==1)
{printf("Error, it does not make sense to use \"--threshold\" with \"--bfile\" (because genotypes are already 0/1/2/NA)\n\n");exit(1);}

if(threshold!=-9999&&genprobs==0)
{printf("Error, it does not make sense to use \"--threshold\" when providing haplotypes (because genotypes are already 0/1/2/NA)\n\n");exit(1);}

if(minprob!=-9999&&genprobs<2)
{printf("Error, you can only use \"--minprob\" when providing genotype probabilities (\"--gen\" or \"--mgen\")\n\n");exit(1);}

//////////////////////////

//kinships, regions, responses, summaries and fixed

if(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0)	//kinships provided (will check number correct in parsefiles.c)
{
if(mode==113)
{printf("Error, you can not use \"--grm\" or \"--mgrm\" with \"--join-kins\"; if you are looking to add kinship matrices, please use \"--add-grm\"\n\n");exit(1);}
if(mode!=115&&mode!=116&&mode!=117&&mode!=118&&mode!=119&&mode!=120&&mode!=121&&mode!=122&&mode!=123&&mode!=124&&mode!=131&&mode!=133&&mode!=138&&mode!=161&&mode!=162&&mode!=163&&mode!=164&&mode!=166&&mode!=167&&mode!=168&&mode!=169&&mode!=170&&mode!=177)
{printf("Warning, kinships are provided but will not be used\n");
strcpy(kinname,"blank");strcpy(kinlist,"blank");}
}
else	//kinships not provided
{
if(mode==115||mode==118||mode==119||mode==161||mode==163||mode==164||mode==166||mode==167||mode==168||mode==169||mode==170)
{printf("Error, you must provide a kinship matrix using \"--grm\"\n\n");exit(1);}
if(mode==116)
{printf("Error, you must provide kinship matrices to add using \"--mgrm\"\n\n");exit(1);}
if(mode==117)
{printf("Error, you must provide either multiple kinship matrices using \"--mgrm\", or one kinship matrix and a list of predictors to extract or exclude using \"--grm\" with \"--extract\" or \"--exclude\"\n\n");exit(1);}
if(mode==120)
{printf("Error, you must provide two or more kinship matrices using \"--mgrm\"\n\n");exit(1);}
if(mode==133)
{printf("Error, you must provide multiple kinship matrices using \"--mgrm\"\n\n");exit(1);}
if(mode==162)
{printf("Error, you must use \"--grm\" to provide the kinship matrix used with \"--pca\"\n\n");exit(1);}

if(strcmp(eigenfile,"blank")!=0)
{printf("Error, you must use \"--grm\" to provide the kinship matrix to which the eigen-decomposition corresponds\n\n");exit(1);}
}

if(mode==113){sprintf(kinlist,"%spartition.list", folder);}

////////

if(kindetails==0&&(mode==113||mode==116||mode==117||mode==122||mode==162))
{printf("Error, you can not use \"--kinship-details NO\" with \"--join-kins\", \"--calc-blups\", \"--add-grm\", \"--sub-grm\" or \"--calc-pca-loads\"\n\n");exit(1);}

if(kindetails==1&&(mode==118||mode==119))
{printf("Error, you can not use \"--kinship-details YES\" with \"--convert-gz\" or \"--convert-raw\"\n\n");exit(1);}

////////

if(num_regs!=-9999||strcmp(regpref,"blank")!=0)	//regions provided
{
if(mode==117)
{printf("Sorry, you can no longer use regions with \"--sub-grm\"; instead use \"--extract\" or \"--exclude\" to specify which predictors to extract or exclude from the kinship matrix\n\n");exit(1);}

if(mode==122)
{printf("Error, you can not use \"--region-number\" or \"--region-prefix\" with \"--calc-blups\" (region details are linked to in the remlfile)\n\n");exit(1);}

if(mode==138)
{printf("Error, you can not use regions with \"--calc-genes-reml\", so you should instead convert each into a kinship matrix (if you have more than one region, you should then use\"--solve-null\" to combine the kinship matrices)\n\n");exit(1);}

if(mode!=121&&mode!=123&&mode!=124)
{printf("Error, you can only use regions with \"--reml\", \"--he\" or \"--pcgc\"\n\n");exit(1);}

if(strcmp(regpref,"blank")==0)
{printf("Error, when using \"--region-number\" you must also use \"--region-prefix\"\n\n");exit(1);}

if(num_regs==-9999)
{printf("Error, when using \"--region-prefix\" you must also use \"--region-number\"\n\n");exit(1);}

if(dtype==-9999)
{printf("Error, when using regions, you must provide a set of genetic data files using \"--bfile\", \"--sp\", \"--sped\", \"--speed\" or \"--gen\"\n\n");exit(1);}
}

if(rprune!=-9999&&strcmp(regpref,"blank")==0)
{printf("Error, you can only use \"--region-prune\" when using regions\"\n\n");exit(1);}

////////

if(strcmp(respfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use both \"--pheno\" and \"--summary\"\n\n");exit(1);}

if(strcmp(respfile,"blank")==0)	//phenotypes not provided
{
if((mode==121||(mode==172&&strcmp(finalfile,"blank")!=0))&&strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--pheno\" or \"--summary\" to provide phenotypes or summary statistics\n\n");exit(1);}
if(mode==138&&strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--pheno\" or \"--summary\" to provide phenotypes or summary statistics (you can not use \"-summary\" if also providing a kinship matrix)\n\n");exit(1);}
if(mode==123||mode==124||mode==126||mode==127||mode==131||mode==132||mode==133||mode==151||mode==152||mode==153)
{printf("Error, you must use \"--pheno\" to provide phenotypes\n\n");exit(1);}
}

if((mpheno!=-9999||strcmp(npheno,"blank")!=0)&&strcmp(respfile,"blank")==0)
{printf("Error, you can only use \"--mpheno\" or \"--pheno-name\" if providing phenotypes\n\n");exit(1);}

if(mpheno!=-9999&&strcmp(npheno,"blank")!=0)
{printf("Errors, you can not use both \"--mpheno\" and \"--pheno-name\"\n\n");exit(1);}

if(mpheno==-1&&mode!=121&&mode!=123&&mode!=124&&mode!=163)
{printf("Error, you can only use \"--mpheno -1\" with \"--reml\", \"--he\", \"--pcgc\" or \"--decompose\"\n\n");exit(1);}

if(pad!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=131&&mode!=138&&mode!=163)
{printf("Error, you can only use \"--dentist\" with \"--reml\", \"--he\", \"--pcgc\", \"--linear\", \"--calc-genes-reml\" or \"--decompose\"\n\n");exit(1);}

////////

if(strcmp(sumsfile,"blank")!=0)	//summaries provided
{
if(mode!=121&&mode!=138&&mode!=146&&mode!=147&&mode!=150&&mode!=158&&mode!=159&&mode!=172)
{printf("Error, you can only use \"--summary\" with \"--reml\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--calc-posts\", \"--mega-prs\", \"--pseudo-summaries\" or \"--calc-scores\"\n\n");exit(1);}

if(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0)
{printf("Error, you can not use \"--summary\" if providing kinships\n\n");exit(1);}

if(strcmp(covarfile,"blank")!=0)
{printf("Error, you can not use \"--covar\" with \"--summary\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--enviro\" with \"--summary\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0)
{printf("Error, you can not use \"--top-preds\" with \"--summary\"\n\n");exit(1);}
}

if(strcmp(sums2file,"blank")!=0&&mode!=147&&mode!=158)
{printf("Error, you can only use \"--summary2\" with \"--sum-cors\" or \"--mega-prs\"\n\n");exit(1);}

if(fixn!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--fixed-n\" when using \"--summary\"\n\n");exit(1);}

if(amb!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--allow-ambiguous\" when using \"--summary\"\n\n");exit(1);}

if(scaling!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--scaling\" when using \"--summary\"\n\n");exit(1);}

if(scaling2!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--scaling2\" when using \"--summary2\"\n\n");exit(1);}

////////

if(prev!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=131&&mode!=133&&mode!=138&&mode!=146&&mode!=173)
{printf("Error, you can only use \"--prevalence\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-he\", \"--fast-pcgc\", \"--linear\", \"--solve-null\", \"--calc-genes-reml\", \"--sum-hers\" or \"--make-phenos\"\n\n");exit(1);}

if(ascer!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--ascertainment\" when using \"--summary\"\n\n");exit(1);}

if(prev!=-9999&&mpheno==-1)
{printf("Warning, you should only use \"--prevalence\" with \"--mpheno -1\" if the prevalence is constant across phenotypes (if it varies, you should instead analyze each phenotype separately)\n\n");}

if(prev!=-9999&&strcmp(sumsfile,"blank")!=0&&ascer==-9999)
{printf("Error, you must use \"--ascertainment\" to specify the proportion of samples who were cases when generating the summary statistics\n\n");exit(1);}

////////

if((mode==169||mode==170)&&strcmp(covarfile,"blank")!=0)
{printf("Error, when using \"--gxemm-iid\" or \"--gxemm-free\" you must use \"--enviro\" instead of \"--covar\"\n\n");exit(1);}

if(strcmp(covarfile,"blank")!=0&&mode!=121&&mode!=122&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=131&&mode!=132&&mode!=133&&mode!=138&&mode!=151&&mode!=152&&mode==153&&mode!=164&&mode!=172)
{printf("Error, you can only use \"--covar\" with \"--reml\", \"--blup\", \"--he\", \"--pcgc\", \"--fast-he\", \"--fast-pcgc\", \"--linear\", \"--logistic\", \"--solve-null\", \"--calc-genes-reml\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--adjust-grm\" or \"--calc-scores\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0&&mode!=121&&mode!=122&&mode!=123&&mode!=124&&mode!=164&&mode!=169&&mode!=170&&mode!=172)
{printf("Error, you can only use \"--enviro\" with \"--reml\", \"--blup\", \"--he\", \"--pcgc\", \"--adjust-grm\", \"--gxemm-iid\", \"--gxemm-free\" or \"--calc-scores\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode==122)
{printf("Error, you can not use \"--top-preds\" with \"--calc-blups\" (top predictors details are linked to in the remlfile)\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode==172)
{printf("Error, you can not use \"--top-preds\" with \"--calc-scores\"; the top predictors should instead be included in the score file\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=131&&mode!=132&&mode!=133&&mode!=138&&mode!=151&&mode!=152&&mode==153&&mode!=164)
{printf("Error, you can only use \"--top-preds\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-he\", \"--fast-pcgc\", \"--linear\", \"--logistic\", \"--solve-null\", \"--calc-genes-reml\", \"--ridge\", \"--bolt\", \"--bayesr\" or \"--adjust-grm\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&dtype==-9999)
{printf("Error, when using \"--top-preds\" you must provide a set of genetic data files using \"--bfile\", \"--sp\", \"--sped\", \"--speed\" or \"--gen\"\n\n");exit(1);}

//////////////////////////

//calculating weights, thinning and finding/removing tags

if(nothin!=-9999&&mode!=101)
{printf("Error, you can only use \"--no-thin\" with \"--cut-weights\"\n\n");exit(1);}

if(wprune!=-9999&&mode!=101&&mode!=106&&mode!=107)
{printf("Error, you can only use \"--window-prune\" with \"--cut-weights\", \"--thin\" or \"--thin-tops\"\n\n");exit(1);}

if(wprune!=-9999&&(nothin==1||nothin==2))
{printf("Error, it does not make sense to use \"--window-prune\" with \"--no-thin YES\" or \"--no-thin DONE\"\n\n");exit(1);}

if(window_kb!=-9999&&mode!=101&&mode!=102&&mode!=104&&mode!=106&&mode!=107&&mode!=108&&mode!=109&&mode!=141&&mode!=156&&mode!=158)
{printf("Error, you can only use \"--window-kb\" with \"--cut-weights\", \"--calc-weights\", \"--calc-weights-all\", \"--calc-tagging\", \"--calc-cors\", \"--mega-prs\", \"--thin\", \"--thin-tops\", \"--find-tags\" or \"--remove-tags\"\n\n");exit(1);}

if(window_length!=-9999&&mode!=101&&mode!=102&&mode!=104&&mode!=106&&mode!=107)
{printf("Error, you can only use \"--window-length\" with \"--cut-weights\", \"--calc-weights\", \"--calc-weights-all\", \"--thin\" or \"--thin-tops\"\n\n");exit(1);}

if(window_cm!=-9999&&mode!=101&&mode!=102&&mode!=104&&mode!=106&&mode!=107&&mode!=108&&mode!=109&&mode!=141&&mode!=156&&mode!=158)
{printf("Error, you can only use \"--window-cm\" with \"--cut-weights\", \"--calc-weights\", \"--calc-weights-all\", \"--calc-tagging\", \"--calc-cors\", \"--mega-prs\", \"--thin\", \"--thin-tops\", \"--find-tags\" or \"--remove-tags\"\n\n");exit(1);}

if(window_kb!=-9999&&window_length!=-9999)
{printf("Error, you can not use both \"--window-kb\" and \"--window-length\"\n\n");exit(1);}

if(window_cm!=-9999&&(window_kb!=-9999||window_length!=-9999))
{printf("Error, when using \"--window-cm\", you can not also use \"--window-kb\" or \"--window-length\"\n\n");exit(1);}

if(window_cm!=-9999&&dtype==5&&strcmp(ubimfile,"blank")==0)
{printf("Error, when using \"--window-cm\" with \"--gen\", you must also use \"--bim\" to provide genetic distances\n\n");exit(1);}

////////

if((section_kb!=-9999||section_length!=-9999||section_cm!=-9999)&&mode!=101)
{printf("Error, you can only use \"--section-kb\", \"--section-length\" or \"--section-cm\" with \"--cut-weights\"\n\n");exit(1);}

if(section_kb!=-9999&&section_length!=-9999)
{printf("Error, you can not use both \"--section-kb\" and \"--section-length\"\n\n");exit(1);}

if(section_cm!=-9999&&(section_kb!=-9999||section_length!=-9999))
{printf("Error, when using \"--section-cm\", you can not use \"--section-kb\" or \"--section-length\"\n\n");exit(1);}

if(section_cm!=-9999&&dtype==5&&strcmp(ubimfile,"blank")==0)
{printf("Error, when using \"--section-cm\" with \"--gen\", you must also use \"--bim\" to provide genetic distances\n\n");exit(1);}

if(section_kb!=-9999&&window_cm!=-9999)
{printf("Error, you can not use both \"--window-cm\" and \"--section-kb\" (generally, you should use either \"--window-kb\" and \"--section-kb\" or \"--window-cm\" and \"--section-cm\")\n\n");exit(1);}

if(section_cm!=-9999&&window_kb!=-9999)
{printf("Error, you can not use both \"--window-kb\" and \"--section-cm\" (generally, you should use either \"--window-kb\" and \"--section-kb\" or \"--window-cm\" and \"--section-cm\")\n\n");exit(1);}

if(section_cm!=-9999&&window_cm==-9999)
{printf("Error, when using \"--section-cm\" you must also use \"--window-cm\"\n\n");exit(1);}

////////

if((buffer_kb!=-9999||buffer_length!=-9999||buffer_cm!=-9999)&&mode!=101)
{printf("Error, you can only use \"--buffer-kb\", \"--buffer-length\" or \"--buffer-cm\" with \"--cut-weights\"\n\n");exit(1);}

if(buffer_kb!=-9999&&buffer_length!=-9999)
{printf("Error, you can not use both \"--buffer-kb\" and \"--buffer-length\"\n\n");exit(1);}

if(buffer_cm!=-9999&&(buffer_kb!=-9999||buffer_length!=-9999))
{printf("Error, when using \"--buffer-cm\", you can not also use \"--buffer-kb\" or \"--buffer-length\"\n\n");exit(1);}

if(buffer_kb!=-9999&&window_cm!=-9999)
{printf("Error, you can not use \"--buffer-kb\" with \"--window-cm\"\n\n");exit(1);}

if(buffer_cm!=-9999&&dtype==5&&strcmp(datalist,"blank")==0&&strcmp(ubimfile,"blank")==0)
{printf("Error, when using \"--buffer-cm\" you must also use \"--bim\" to provide genetic distances\n\n");exit(1);}

if(buffer_cm!=-9999&&window_cm==-9999)
{printf("Error, when using \"--buffer-cm\" you must also use \"--window-cm\"\n\n");exit(1);}

if(buffer_kb!=-9999||buffer_length!=-9999||buffer_cm!=-9999)
{printf("Warning, there is no need to use \"--buffer-kb\", \"--buffer-length\" or \"--buffer-cm\"; by default these are set from \"--window-kb\", \"--window-length\" and \"--window-cm\"\n\n");}

////////

if(section!=-9999&&mode==104)
{printf("Error, you can not use \"--section\" with \"--calc-weights-all\" (LDAK will loop through all sections)\n\n");exit(1);}

if(section!=-9999&&mode!=102)
{printf("Error, you can only use \"--section\" with \"--calc-weights\"\n\n");exit(1);}

if(section_start!=-9999&&mode!=104)
{printf("Error, you can only use \"--start-section\" with \"--calc-weights-all\"\n\n");exit(1);}

////////

if(strcmp(infosfile,"blank")!=0&&mode!=101&&mode!=102&&mode!=103&&mode!=104)
{printf("Error, you can now only use \"--infos\" when calculating weightings (if you wish to include infos when calculating kinships, you should incorporate them in the weightings)\n\n");exit(1);}

if((lddecay!=-9999||halflife!=-9999)&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--decay\" or \"--half-life\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(halflife!=-9999&&lddecay==0)
{printf("Error, you can only use \"--half-life\" with \"--decay YES\"\n\n");exit(1);}

////////

if(fudge!=-9999&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--quick-weights\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(simplex!=-9999&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--simplex\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(maxtime!=-9999&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--max-time\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(maxtime!=-9999&&simplex!=1)
{printf("Error, you can only use \"--max-time\" with \"--simplex YES\"\n\n");exit(1);}

if(spread!=-9999&&mode!=103&&mode!=104)
{printf("Error, you can only use \"--spread\" with \"--join-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

////////

if(strcmp(targetfile,"blank")!=0&&mode!=109)
{printf("Error, you can only use \"--targets\" with \"--remove-tags\"\n\n");exit(1);}

//////////////////////////

//calculating and manipulating kinships

if(part_length!=-9999&&mode!=111&&mode!=136)
{printf("Error, you can only use \"--partition-length\" with \"--cut-kins\" or \"--cut-genes\"\n\n");exit(1);}

if(bychr!=-9999&&mode!=111&&mode!=136)
{printf("Error, you can only use \"--by-chr\" with \"--cut-kins\" or \"--cut-genes\"\n\n");exit(1);}

if(part_length!=-9999&&bychr==1)
{printf("Error, you can not use both \"--partition-length\" and \"--by-chr YES\"\n\n");exit(1);}

if(num_parts!=-9999||strcmp(partpref,"blank")!=0)	//partitions provided
{
if(mode!=105&&mode!=111&&mode!=126&&mode!=127&&mode!=141)
{printf("Error, you can only use partitions with \"--cut-kins\", \"--fast-he\", \"--fast-pcgc\" or \"--calc-tagging\"\n\n");exit(1);}

if(part_length!=-9999||bychr==1)
{printf("Error, you can not use \"--partition-length\" or \"--by-chr YES\" when using partitions\n\n");exit(1);}

if(strcmp(partpref,"blank")==0)
{printf("Error, when using \"--partition-number\" you must also use \"--partition-prefix\"\n\n");exit(1);}

if(num_parts==-9999)
{printf("Error, when using \"--partition-prefix\" you must also use \"--partition-number\"\n\n");exit(1);}

if(num_parts==1&&mode==111)
{printf("Warning, using \"--partition-number 1\" and \"--partition-prefix %s\" is equivalent to (simply) using \"--calc-kins-direct\" with \"--extract %s1\"\n\n", partpref, partpref);}
}

if(checkpart!=-9999&&mode!=111)
{printf("Error, you can only use \"--check-partitions\" with \"--cut-kins\"\n\n");exit(1);}

if(checkpart!=-9999&&strcmp(partpref,"blank")==0)
{printf("Error, you can only use \"--check-partitions\" when using partitions\n\n");exit(1);}

////////

if(partition!=-9999&&mode!=112&&mode!=137&&mode!=138)
{printf("Error, you can only use \"--partition\" with \"--calc-kins\", \"--calc-genes-kins\" or \"--calc-genes-reml\"\n\n");exit(1);}

if((kingz!=-9999||kinraw!=-9999)&&mode!=112&&mode!=113&&mode!=114&&mode!=116&&mode!=117&&mode!=118&&mode!=119&&mode!=133&&mode!=137&&mode!=164&&mode!=166&&mode!=167&&mode!=168&&mode!=169&&mode!=170)
{printf("Error, you can only use \"--kinship-gz\" or \"--kinship-raw\" when producing kinship matrices\n\n");exit(1);}

if(single!=-9999&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--single\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

////////

if((maxrel!=-9999||minrel!=-9999)&&mode!=115)
{printf("Error, you can only use \"--max-rel\" or \"--min-rel\" with \"--filter\"\n\n");exit(1);}

if(maxrel!=-9999&&minrel!=-9999)
{printf("Error, you can not use both \"--max-rel\" and \"--min-rel\"\n\n");exit(1);}

if(kinstand!=-9999&&mode!=115)
{printf("Error, you can only use \"--kin-stand\" with \"--filter\"\n\n");exit(1);}

///////////////////////////

//reml, blup, he and sum-hers

if(shortcut!=-9999&&mode!=121)
{printf("Error, you can only use \"--shortcut\" with \"--reml\"\n\n");exit(1);}

if(shortcut!=-9999&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--shortcut\" with \"--summary\"\n\n");exit(1);}

if(compact!=-9999&&mode!=121)
{printf("Error, you can only use \"--compact\" with \"--reml\"\n\n");exit(1);}

if(compact==1)
{printf("Warning, you should not use \"--compact YES\" if you plan to subsequently calculate BLUPs\n\n");}

if(discenv==1&&mode!=121&&mode!=123&&mode!=124&&mode!=170)
{printf("Error, you can only use \"--subgroups\" with \"--reml\", \"--he\", \"--pcgc\" or \"--gxemm-free\"\n\n");exit(1);}

if(discenv==1&&strcmp(envfile,"blank")==0)
{printf("Error, you can only use \"--subgroups\" when using \"--enviro\"\n\n");exit(1);}

if(liab!=-9999&&mode!=121)
{printf("Error, you can only use \"--liability\" with \"--reml\"\n\n");exit(1);}

if(liab==1&&mpheno==-1)
{printf("Error, you can not use \"--mpheno -1\" with \"--liability YES\"\n\n");exit(1);}

if(liab==1&&pad==1)
{printf("Error, you can not use \"--dentist YES\" with \"--liability YES\"\n\n");exit(1);}

////////

if(strcmp(remlfile,"blank")!=0&&mode!=122)
{printf("Error, you can only use \"--remlfile\" with \"--calc-blups\"\n\n");exit(1);}

////////

if(adjusted!=-9999&&mode!=123&&mode!=124)
{printf("Error, you can only use \"--adjusted\" with \"--he\" or \"--pcgc\"\n\n");exit(1);}

////////

if(num_vects!=-9999&&mode!=126&&mode!=127)
{printf("Error, you can only use \"--repetitions\" with \"--fast-he\" or \"--fast-pcgc\"\n\n");exit(1);}

///////////////////////////

//association analysis

if(exact!=-9999&&mode!=131)
{printf("Error, you can only use \"--exact\" with \"--linear\"\n\n");exit(1);}

if(exact==0&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0)
{printf("Error, you can only use \"--exact NO\" when providing a kinship matrix\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&mode!=131)
{printf("Error, you can only use \"--sample-weights\" with \"--linear\"\n\n");exit(1);}

////////

if((strcmp(genefile,"blank")!=0||chunks!=-9999||chunksbp!=-9999)&&mode!=136&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--genefile\", \"--chunks\" or \"--chunks-bp\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if((gene_buffer!=-9999||up_buffer!=-9999||down_buffer!=-9999)&&mode!=136&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--gene-buffer\", \"--up-buffer\" or \"--down-buffer\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if((gene_buffer!=-9999||up_buffer!=-9999||down_buffer!=-9999)&&strcmp(genefile,"blank")==0)
{printf("Error, you can only use \"--gene-buffer\", \"--up-buffer\" or \"--down-buffer\" when using \"--genefile\" to provide an annotation file\n\n");exit(1);}

if((up_buffer!=-9999||down_buffer!=-9999)&&gene_buffer!=-9999)
{printf("Error, you can not use \"--up-buffer\" or \"--down-buffer\" if using \"--gene-buffer\"\n\n");exit(1);}

if(up_buffer!=-9999&&down_buffer==-9999)
{printf("Error, when using \"--up-buffer\", you must also use \"--down-buffer\"\n\n");exit(1);}

if(down_buffer!=-9999&&up_buffer==-9999)
{printf("Error, when using \"--down-buffer\", you must also use \"--up-buffer\"\n\n");exit(1);}

if(minweight!=-9999&&mode!=136&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--min-weights\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if(minweight!=-9999&&chunks!=-9999)
{printf("Error, you can not use \"--min-weight\" when using \"--chunks\" (the genome will be divided so that all chunks have equal weight)\n\n");exit(1);}

if(overlap!=-9999&&mode!=136&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--overlap\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

////////

if(gene_perms!=-9999&&mode!=138)
{printf("Error, you can only use \"--gene-permutations\" with \"--calc-genes-reml\"\n\n");exit(1);}

if(gprune!=-9999&&mode!=138)
{printf("Error, you can only use \"--gene-prune\" with \"--calc-genes-reml\"\n\n");exit(1);}

////////

if(cut1!=-9999||cut2!=-9999)
{
if(mode!=139)
{printf("Error, you can only use \"--sig1\" and \"--sig2\" with \"--join-genes-reml\"\n\n");exit(1);}

if(cut1==-9999||cut2==-9999)
{printf("Error, to identify significant genes/chunks, you must use both \"--sig1\" and \"--sig2 (typical values are 0.00001 and 0.01)\n\n");exit(1);}

if(cut2<cut1){printf("Warning, \"--sig2\" (%f) is lower than \"--sig1\" (%f), which can result in overlapping regions\n\n", cut2, cut2);}
}

if(gamp!=-9999||gam1!=-9999||gam2!=-9999)
{
if(mode!=139)
{printf("Error, you can only use \"--gamma-fraction\", \"--gamma-alpha\" or \"--gamma-beta\" with \"--join-genes-reml\"\n\n");exit(1);}
if(gamp==-9999||gam1==-9999||gam2==-9999)
{printf("Error, to specify the NULL distribution of LRT statistics, you must provide \"--gamma-fraction\", \"--gamma-alpha\" and \"--gamma-beta\"\n\n");exit(1);}
}

///////////////////////////

//sumher

if(num_anns!=-9999||strcmp(annpref,"blank")!=0)	//annotations provided
{
if(mode!=126&&mode!=127&&mode!=141)
{printf("Error, you can only use annotations with \"--fast-he\", \"--fast-pcgc\" or  \"--calc-tagging\"\n\n");exit(1);}

if(num_parts!=-9999||strcmp(partpref,"blank")!=0)
{printf("Error, you can not use both partitions and annotations; generally, you should use partitions when the lists of predictors are disjoint, and annotations when there is overlap\n\n");exit(1);}

if(strcmp(annpref,"blank")==0)
{printf("Error, when using \"--annotation-number\" you must also use \"--annotation-prefix\"\n\n");exit(1);}

if(num_anns==-9999)
{printf("Error, when using \"--annotation-prefix\" you must also use \"--annotation-number\"\n\n");exit(1);}
}

if(backpart!=-9999&&mode!=105&&mode!=126&&mode!=127&&mode!=141)
{printf("Error, you can only use \"--background\" with \"--fast-he\", \"--fast-pcgc\" or \"--calc-tagging\"\n\n");exit(1);}

if(backpart!=-9999&&strcmp(partpref,"blank")==0)
{printf("Error, you can only use \"--background\" when using partitions\n\n");exit(1);}

if(allone!=-9999&&mode!=105&&mode!=126&&mode!=127&&mode!=141)
{printf("Error, you can only use \"--all-one YES\" with \"--fast-he\", \"--fast-pcgc\" or \"--calc-tagging\"\n\n");exit(1);}

if(allone==1&&strcmp(partpref,"blank")==0)
{printf("Error, you can only use \"--all-one YES\" when using partitions\n\n");exit(1);}

if(reduce!=-9999&&mode!=105&&mode!=141)
{printf("Error, you can only use \"--reduce\" with \"--calc-tagging\"\n\n");exit(1);}

////////

if(strcmp(printfile,"blank")!=0&&mode!=141)
{printf("Error, you can only use \"--regression-predictors\" with \"--calc-tagging\"\n\n");exit(1);}

if(strcmp(herfile,"blank")!=0&&mode!=141)
{printf("Error, you can only use \"--heritability-predictors\" with \"--calc-tagging\"\n\n");exit(1);}

if(unbias!=-9999&&mode!=141)
{printf("Error, you can only use \"--unbiased\" with \"--calc-tagging\"\n\n");exit(1);}

if(savemat!=-9999&&mode!=141)
{printf("Error, you can only use \"--save-matrix\" with \"--calc-tagging\"\n\n");exit(1);}

if(cover!=-9999&&mode!=141)
{printf("Error, you can only use \"--coverage\" with \"--calc-tagging\"\n\n");exit(1);}

////////

if(strcmp(taglist,"blank")!=0&&mode!=142&&mode!=143)
{printf("Error, you can only use \"--taglist\" with \"--join-tagging\" or \"--merge-tagging\"\n\n");exit(1);}

if(strcmp(matlist,"blank")!=0&&mode!=142)
{printf("Error, you can only use \"--matlist\" with \"--join-tagging\"\n\n");exit(1);}

if(checkdups!=-9999&&mode!=142)
{printf("Error, you can only use \"--check-dups\" with \"--join-tagging\"\n\n");exit(1);}

////////

if(strcmp(tagfile,"blank")!=0&&mode!=144&&mode!=146&&mode!=147&&mode!=149)
{printf("Error, you can only use \"--tagfile\" with \"--reduce-tagging\", \"--sum-hers\", \"--sum-cors\" or \"--calc-exps\"\n\n");exit(1);}

if(strcmp(altfile,"blank")!=0&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--alternative-tags\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(strcmp(cvsfile,"blank")!=0&&mode!=146)
{printf("Error, you can only use \"--cv-predictors\" with \"--sum-hers\"\n\n");exit(1);}

if(strcmp(catfile,"blank")!=0&&mode!=144&&mode!=146)
{printf("Error, you can only use \"--categories\" with \"--reduce-tagging\" or \"--sum-hers\"\n\n");exit(1);}

if(strcmp(taufile,"blank")!=0&&mode!=146&&mode!=148&&mode!=149)
{printf("Error, you can only use \"--taus\" with \"--sum-hers\", \"--calc-ind-hers\" or \"--calc-exps\"\n\n");exit(1);}

if(strcmp(catfile,"blank")!=0&&strcmp(taufile,"blank")!=0)
{printf("Error, you can not use both \"--categories\" and \"--taus\"\n\n");exit(1);}

////////

if(checksums!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--check-sums\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(gcon!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--genomic-control\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(cept!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--intercept\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

////////

if(ldsc!=-9999&&mode!=146)
{printf("Error, you can only use \"--LDSC\" with \"--sum-hers\"\n\n");exit(1);}

if(ldsc==1&&(strcmp(catfile,"blank")!=0||strcmp(taufile,"blank")!=0))
{printf("Error, you can not use \"--LDSC YES\" with \"--categories\" or \"--taus\"\n\n");exit(1);}

if(chisol!=-9999&&mode!=146)
{printf("Error, you can only use \"--chisq-solver\" with \"--sum-hers\"\n\n");exit(1);}

if(ldsc==1&&chisol==1)
{printf("Error, you can not both \"--LDSC YES\" and \"--chisq-solver YES\"\n\n");exit(1);}

if(tagone!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--tag-one\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

////////

if(divide!=-9999&&mode!=146)
{printf("Error, you can only use \"--divisions\" with \"--sum-hers\"\n\n");exit(1);}

if(divide!=-9999&&(strcmp(catfile,"blank")!=0||strcmp(taufile,"blank")!=0))
{printf("Error, you can not use \"--divisions\" with \"--categories\" or \"--taus\"\n\n");exit(1);}

if(divide!=-9999&&ldsc==1)
{printf("Error, you can not use \"--divisions\" with \"--LDSC YES\"\n\n");exit(1);}

if(uptaus!=-9999&&mode!=146)
{printf("Error, you can only use \"--update-taus\" with \"--sum-hers\"\n\n");exit(1);}

if(uptaus!=-9999&&divide==-9999)
{printf("Error, you can only use \"--update-taus\" with \"--divisions\"\n\n");exit(1);}

if(strcmp(powfile,"blank")!=0&&mode!=146)
{printf("Error, you can only use \"--powerfile\" with \"--sum-hers\"\n\n");exit(1);}

if(strcmp(powfile,"blank")!=0&&divide==-9999)
{printf("Error, you can only use \"--powerfile\" with \"--divisions\"\n\n");exit(1);}

////////

if(plet!=-9999&&mode!=147)
{printf("Error, you can only use \"--absolutes\" with \"--sum-cors\"\n\n");exit(1);}

////////

if(strcmp(matfile,"blank")!=0&&mode!=146&&mode!=148)
{printf("Error, you can only use \"--matrix\" with \"--sum-hers\" or \"--calc-ind-hers\"\n\n");exit(1);}

if(strcmp(matfile,"blank")!=0&&divide!=-9999)
{printf("Error, you can not use \"--matrix\" with \"--divisions\" (if you wish to estimate per-predictor heritabilities, then you should run \"--calc-ind-hers\" once for each set of results from this analysis)\n\n");exit(1);}

////////

if(strcmp(propfile,"blank")!=0&&mode!=149)
{printf("Error, you can only use \"--props\" with \"--calc-exps\"\n\n");exit(1);}

if(strcmp(propfile,"blank")!=0&&her==-9999)
{printf("Error, you can only use \"--props\" with \"--her\"\n\n");exit(1);}

if(strcmp(expfile,"blank")!=0&&mode!=150)
{printf("Error, you can only use \"--expectations\" with \"--calc-posts\"\n\n");exit(1);}

//////////////////////////

//individual-level data prediction, then megaprs

if(strcmp(indhers,"blank")!=0&&mode!=151&&mode!=152&&mode!=153&&mode!=158)
{printf("Error, you can only use \"--ind-hers\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--mega-prs\"\n\n");exit(1);}

if(strcmp(indhers,"blank")!=0&&(strcmp(weightsfile,"blank")!=0||ignoreweights!=-9999||power!=-9999||her!=-9999))
{printf("Error, when using \"--ind-hers\" you should not use \"--weights\", \"--ignore-weights\", \"--power\" or \"--her\" (\"--ind-hers %s\" is equivalent to using \"--weights %s\" with \"--power -1\" and her equal to the sum of the per-predictor heritabilities)\n\n", indhers2, indhers2);exit(1);}

if(herscale!=-9999&&strcmp(indhers,"blank")!=0)
{printf("Error, you can only use \"--her-scale\" with \"--ind-hers\"\n\n");exit(1);}

////////

if(restart!=-9999&&mode!=151&&mode!=152&&mode!=153)
{printf("Error, you can only use \"--restart\" with \"--ridge\", \"--bolt\" or \"--bayesr\"\n\n");exit(1);}

if(cvprop!=-9999&&mode!=152&&mode!=153)
{printf("Error, you can only use \"--cv-proportion\" with \"--bolt\" or \"--bayesr\"\n\n");exit(1);}

if(strcmp(bvsfile,"blank")!=0&&mode!=152&&mode!=153)
{printf("Error, you can only use \"--cv-samples\" with \"--bolt\" or \"--bayesr\"\n\n");exit(1);}

if(cvprop!=-9999&&strcmp(bvsfile,"blank")!=0)
{printf("Error, you can not use both \"--cv-proportion\" and \"--cv-samples\"\n\n");exit(1);}

if(skipcv!=-9999&&mode!=152&&mode!=153)
{printf("Error, you can only use \"--skip-cv\" with \"--bolt\" or \"--bayesr\"\n\n");exit(1);}

if(skipcv==1&&(cvprop!=-9999||strcmp(bvsfile,"blank")!=0))
{printf("Error, it is contradictory to use \"--skip-cv YES\" with \"--cv-proportion\" or \"--cv-samples\"\n\n");exit(1);}

////////

if(pointmass!=-9999&&mode!=153)
{printf("Error, you can only use \"--point-mass\" with \"--bayesr\"\n\n");exit(1);}

if(fullspace!=-9999&&mode!=153)
{printf("Error, you can only use \"--full-space\" with \"--bayesr\"\n\n");exit(1);}

////////

if(strcmp(corslist,"blank")!=0&&mode!=157)
{printf("Error, you can only use \"--corslist\" with \"--join-cors\"\n\n");exit(1);}

////////

if(strcmp(corsfile,"blank")!=0&&mode!=158)
{printf("Error, you can only use \"--cors\" with \"--mega-prs\"\n\n");exit(1);}

if(segments!=-9999&&mode!=158)
{printf("Error, you can only use \"--segments\" with \"--mega-prs\"\n\n");exit(1);}

if(ptype!=-9999&&mode!=158)
{printf("Error, you can only use \"--model\" with \"--mega-prs\"\n\n");exit(1);}

if(checksize!=-9999&&mode!=158)
{printf("Error, you can only use \"--check-sizes\" with \"--mega-prs\"\n\n");exit(1);}

////////

if(multiher!=-9999&&mode!=158)
{printf("Error, you can only use \"--multi-hers\" with \"--mega-prs\"\n\n");exit(1);}

if(basicbolt!=-9999&&mode!=158)
{printf("Error, you can only use \"--basic-bolt\" with \"--mega-prs\"\n\n");exit(1);}

if(ldpred!=-9999&&mode!=158)
{printf("Error, you can only use \"--LDPred\" with \"--mega-prs\"\n\n");exit(1);}

if(basicbolt==1&&ldpred==1)
{printf("Error, you can not use \"--basic-bolt YES\" with \"--LDPred YES\"\n\n");exit(1);}

////////

if(subprop!=-9999&&mode!=159)
{printf("Error, you can only use \"--training-proportion\" with \"--pseudo-summaries\"\n\n");exit(1);}

///////////////////////////

//pca and decompose

if(axes!=-9999&&mode!=161&&mode!=167)
{printf("Error, you can only use \"--axes\" with \"--pca\"\n\n");exit(1);}

if(strcmp(pcastem,"blank")!=0&&mode!=162)
{printf("Error, you can only use \"--pcastem\" with \"--calc-pca-loads\"\n\n");exit(1);}

////////

if(eigenraw!=-9999&&mode!=163)
{printf("Error, you can only use \"--eigen-raw\" with \"--decompose\"\n\n");exit(1);}

if(strcmp(eigenfile,"blank")!=0&&mode!=121&&mode!=131&&mode!=138&&mode!=161)
{printf("Error, you can only use \"--eigen\" with \"--reml\", \"--linear\", \"--calc-genes-reml\" or \"--pca\"\n\n");exit(1);}

///////////////////////////

//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

if(strcmp(scorefile,"blank")!=0&&mode!=108&&mode!=172)
{printf("Error, you can only use \"--scorefile\" with \"--find-tags\" or \"--calc-scores\"\n\n");exit(1);}

if(strcmp(cofile,"blank")!=0&&mode==122)
{printf("Error, you can not use \"--coeffsfile\" with \"--calc-blups\" (covariates details are linked to in the remlfile)\n\n");exit(1);}

if(strcmp(cofile,"blank")!=0&&mode!=172)
{printf("Error, you can only use \"--coeffsfile\" with \"--calc-scores\"\n\n");exit(1);}

if(savecounts!=-9999&&mode!=172)
{printf("Error, you can only use \"--save-counts\" with \"--calc-scores\"\n\n");exit(1);}

if(strcmp(finalfile,"blank")!=0&&mode!=172)
{printf("Error, you can only use \"--final-effects\" with \"--calc-scores\"\n\n");exit(1);}

////////

if((num_phenos!=-9999||num_causals!=-9999)&&mode!=173)
{printf("Error, you can only use \"--num-phenos\" or \"--num-causals\" with \"--make-phenos\"\n\n");exit(1);}

if(her>1&&mode!=149&&mode!=158)
{printf("Error, you can only use \"--her-big\" with \"--calc-exps\" or \"--mega-prs\"\n\n");exit(1);}

if(her!=-9999&&mode!=149&&mode!=151&&mode!=152&&mode!=153&&mode!=158&&mode!=173)
{printf("Error, you can only use \"--her\" with \"--calc-exps\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--mega-prs\" or \"--make-phenos\"n\n");exit(1);}

if((mode==149||mode==151||mode==152||mode==153||mode==158)&&her==0)
{printf("Error, \"--her\" should be followed by a positive float (not 0)\n\n");exit(1);}

if(her!=-9999&&strcmp(taufile,"blank")!=0)
{printf("Error, you can not use both \"--taus\" and \"--her\"\n\n");exit(1);}

if(bivar!=-9999&&mode!=173)
{printf("Error, you can only use \"--bivar\" with \"--make-phenos\"\n\n");exit(1);}

if((strcmp(causalsfile,"blank")!=0||strcmp(effectsfile,"blank")!=0)&&mode!=173)
{printf("Error, you can only use \"--causals\" or \"--effects\" with \"--make-phenos\"\n\n");exit(1);}

if(strcmp(causalsfile,"blank")!=0&&extract==1)
{printf("Error, when using \"--causals\" you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\"\n\n");exit(1);}

if(strcmp(effectsfile,"blank")!=0&&bivar!=-9999)
{printf("Error, you can not use \"--effects\" with \"--bivar\"\n\n");exit(1);}

////////

if((num_inds!=-9999||num_snps!=-9999)&&mode!=174)
{printf("Error, you can only use \"--num-samples\" or \"--num-snps\" with \"--make-snps\"\n\n");exit(1);}

if((maf1!=-9999||maf2!=-9999)&&mode!=174)
{printf("Error, you can only use \"--maf-low\" or \"--maf-high\" with \"--make-snps\"\n\n");exit(1);}

if(maf1!=-9999&&maf2!=-9999&&maf1>maf2)
{printf("Error, maf-low (%.6f) can not be higher than maf-high (%.6f)\n\n", maf1, maf2);exit(1);}

////////

if(strcmp(jackfile,"blank")!=0&&mode!=176)
{printf("Error, you can only use \"--data-pairs\" with \"--jackknife\"\n\n");exit(1);}

if(strcmp(proffile,"blank")!=0&&mode!=172&&mode!=176)
{printf("Error, you can only use \"--profile\" with \"--calc-scores\" or \"--jackknife\"\n\n");exit(1);}

if(strcmp(jackfile,"blank")!=0&&strcmp(proffile,"blank")!=0)
{printf("Error, you can not use both \"--data-pairs\" and \"--profile\"\n\n");exit(1);}

if(mode==172&&strcmp(proffile,"blank")!=0&&strcmp(finalfile,"blank")==0)
{printf("Error, you can only use \"--profile\" if also using \"--final-effects\"\n\n");exit(1);}

if(auc!=-9999&&mode!=176)
{printf("Error, you can only use \"--AUC\" with \"--jackknife\"\n\n");exit(1);}

////////

if(num_folds!=-9999&&mode!=177)
{printf("Error, you can only use \"--num-folds\" with \"--cut-folds\"\n\n");exit(1);}

////////

if(strcmp(likefile,"blank")!=0&&mode!=178)
{printf("Error, you can only use \"--likelihoods\" with \"--find-gaussian\"\n\n");exit(1);}

///////////////////////////

//making and condensing data

if(usenames!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--use-names\" when making data\n\n");exit(1);}

if((comsamps!=-9999||compreds!=-9999)&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--common-samples\" or \"--common-preds\" when making data\n\n");exit(1);}

if(comsamps==1&&compreds==1)
{printf("Error, you can not use both \"--common-samples YES\" and \"--common-preds YES\" (if you wish to test concordance between two datasets, use \"--calc-sim-data\")\n\n");exit(1);}

if(compreds==0&&usenames==1)
{printf("Error, you can not use \"--common-preds NO\" with \"--use-names YES\"\n\n");exit(1);}

if((exsame!=-9999||exdups!=-9999)&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--exclude-same\" or \"--exclude-dups\" when making data\n\n");exit(1);}

if(exsame==1&&usenames==1)
{printf("Error, you can not use \"--exclude-same YES\" with \"--use-names YES\"\n\n");exit(1);}

if(passall!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--pass-all\" when making data\n\n");exit(1);}

if(passall!=-9999&&minmaf==-9999&&maxmaf==-9999&&minvar==-9999&&minobs==-9999&&mininfo==-9999)
{printf("Error, you can only \"--pass-all\" when using \"--min-maf\", \"--max-maf\", \"--min-var\", \"--min-obs\" or \"--min-info\"\n\n");exit(1);}

if(passall==1&&compreds!=1)
{printf("Error, when using \"--pass-all YES\" you must also use \"--common-preds YES\"\n\n");exit(1);}

if(speedlong!=-9999&&mode!=184&&mode!=189)
{printf("Error, you can only use \"--speed-long\" with \"--make-speed\" or \"--condense-speed\"\n\n");exit(1);}

////////

if(useminor!=-9999&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--count-minor\" when condensing data\n\n");exit(1);}

if(useminor!=-9999&&nonsnp==1)
{printf("Error, you can not use \"--count-minor\" with \"--SNP-data NO\"\n\n");exit(1);}

//////////////////////////

//common options

if(checkroot!=-9999&&mode!=122&&mode!=123&&mode!=124&&mode!=157&&mode!=162&&strcmp(corsfile,"blank")==0&&strcmp(eigenfile,"blank")==0)
{printf("Error, you can only use \"--check-root\" with \"--calc-blups\", \"--he\", \"--pcgc\", \"--join-cors\" or \"--calc-pca-loads\", or when using \"--cors\" or \"--eigen\"\n\n");exit(1);}

if(mincor!=-9999&&mode!=102&&mode!=104&&mode!=108&&mode!=109&&mode!=131&&mode!=132&&mode!=141&&mode!=156)
{printf("Error, you can only use \"--min-cor\" with \"--calc-weights\", \"--calc-weights-all\", \"--linear\", \"--logistic\", \"--calc-tagging\", \"--calc-cors\", \"--find-tags\" or \"--remove-tags\"\n\n");exit(1);}

if(cutoff!=-9999&&mode!=107&&mode!=146&&mode!=147&&mode!=166)
{printf("Error, you can only use \"--cutoff\" with \"--sum-hers\", \"--sum-cors\" or \"--thin-tops\"\n\n");exit(1);}

if(constrain!=-9999&&mode!=121&&mode!=131&&mode!=133)
{printf("Error, you can only use \"--constrain\" with \"--reml\", \"--linear\" or \"--solve-null\"\n\n");exit(1);}

if(num_blocks!=-9999&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=146&&mode!=147&&mode!=176)
{printf("Error, you can only use \"--num-blocks\" with \"--he\", \"--pcgc\", \"--fast-he\", \"--fast-pcgc\", \"--sum-hers\", \"--sum-cors\" or \"--jackknife\"\n\n");exit(1);}

if(num_blocks==-1&&mode!=123&&mode!=124)
{printf("Error, you can only use \"--num-blocks -1\" with \"--he\" or \"--pcgc\"\n\n");exit(1);}

if(permute!=-9999&&mode==133)
{printf("Error, you can not use \"--permute\" with \"--solve-null\"; you should instead permute the phenotypes manually\n\n");exit(1);}

if(permute!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=131&&mode!=132&&mode!=138)
{printf("Error, you can only use \"--permute\" with \"--reml\", \"--he\", \"--pcgc\",  \"--fast-he\", \"--fast-pcgc\", \"--linear\", \"--logistic\" or \"--calc-genes-reml\"\n\n");exit(1);}

if(shrink!=-9999&&mode!=121&&mode!=138&&mode!=158)
{printf("Error, you can only use \"--shrink\" with \"--reml\", \"--calc-genes-reml\" or \"--mega-prs\"\n\n");exit(1);}

if(strip!=-9999&&mode!=121&&mode!=138)
{printf("Error, you can only use \"--strip\" with \"--reml\" or  \"--calc-genes-reml\"\n\n");exit(1);}

if(shrink!=-9999&&strip!=-9999)
{printf("Error, you can not use both \"--shrink\" and \"--strip\"\n\n");exit(1);}

////////

if(tol!=-9999&&mode!=121&&mode!=124&&mode!=129&&mode!=131&&mode!=132&&mode!=133&&mode!=138&&mode!=146&&mode!=147&&mode!=151&&mode!=152&&mode!=153&&mode!=158&&mode!=161&&mode!=167)
{printf("Error, you can only use \"--tolerance\" with \"--reml\", \"--pcgc\", \"--fast-pcgc\", \"--linear\", \"--logistic\", \"--solve-null\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--mega-prs\" or \"--pca\"\n\n");exit(1);}

if(maxiter!=-9999&&mode!=102&&mode!=104&&mode!=121&&mode!=124&&mode!=129&&mode!=131&&mode!=132&&mode!=133&&mode!=138&&mode!=146&&mode!=147&&mode!=151&&mode!=152&&mode!=153&&mode!=158)
{printf("Error, you can only use \"--max-iter\" with \"--calc-weights\", \"--calc-weights-all\", \"--reml\",  \"--linear\", \"--logistic\", \"--solve-null\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--ridge\", \"--bolt\", \"--bayesr\" or \"--mega-prs\"\n\n");exit(1);}

if((mode==102||mode==104)&&maxiter!=-9999&&simplex!=1)
{printf("Error, you can only use \"--max-iter\" with \"--simplex YES\"\n\n");exit(1);}

if(memsave!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=133)
{printf("Error, you can only use \"--memory-save\" with \"--reml\", \"--he\", \"--pcgc\" or \"--solve-null\"\n\n");exit(1);}

//////////////////////////

