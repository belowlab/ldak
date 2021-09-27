//To-do list:

//no longer allowed to use prev with logistic

//when doing gene-based and reml with summaries, investigate values for strip and shrink
//add in dentist
//he and pcgc do not allow null model (just covars) correct?
//when writing trace, add trace check to parsefiles.c
//do not allow annotations or partitions with fast-he/pcgc - either add or make sure know which is being used (from root)
//do not allow covariates when calculating traces ... so should add in adjust-traces

//Some notes:

//we assume resp, covars, X, are not large enough to require size_t when indexing
//normal (two-sided) pvalue for T is erfc(|T| root(1/2)), chisq pvalue for S is erfc(S^.5 root(1/2))
//no longer use projections, but (most) old code in project.c

//reml integrates out fixed effs, then solve vars, then given vars, finds maximising fixed + random effs
//when vars fixed, integrating out fixed effects same (up to constant) as setting fixed effects to mle
//therefore random effects are same whether obtain with fixed effects set to mle or integrated out

/*
Copyright 2020 Doug Speed.
 
    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

//Main file - see mode.c for details of the different modes

///////////////////////////

//The pre-compiled versions use Intel MKL libraries

//I compile the serial version using
//gcc -O3 -o ldak5.1.linux source/ldak.c source/libqsopt.linux.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -m64 -I${MKLROOT}/include -lz -static

//and the parallel version using
//gcc -O3 -o ldak5.1.linux.fast source/ldak.c source/libqsopt.linux.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -m64 -I${MKLROOT}/include -lz -fopenmp -static

////////

//If you wish to compile yourself, but do not have Intel Libraries, you must change the variable MKL to zero (look a few lines below)

//On my personal (LINUX) laptop I can compile a dynamic version using
//gcc source/ldak.c source/libqsopt.linux.a -o ldak.out -lblas -llapack -lm -lz

//for a static version, consider something like
//gcc source/ldak.c source/libqsopt.linux.a -o ldak.out -lc libraries/liblapack.a libraries/libblas.a libraries/libc.a libraries/libgfortran.a libraries/libm.a libraries/libz.a -static

//On a MAC, compile using 
//gcc source/ldak.c source/libqsopt.mac.a -o ldak.out -lblas -llapack -lm -lz
//you may have to add --framework accelerate to this this command, and/or first run xcode-select --install

///////////////////////////

#define MET 0	//0 for qsopt, 1 for glpk (for glpk, edit compilation line to include source/glpk.h)
#define MKL 0	//1 to compile with mkl, 0 to compile without mkl (changes one line below and one in defaults.c)

//library includes

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <zlib.h>
#include <limits.h>

#if MET==0
#include "qsopt.h"
#else
#include"glpk.h"
#endif

#if MKL==1
#include<mkl.h>
#endif


//will export functions from fortran

extern void dgemm_();
extern void sgemm_();
extern void dgemv_();
extern void dpotrf_();
extern void spotrf_();
extern void dpotri_();
extern void spotri_();
extern void dpotrs_();
extern void spotrs_();
extern void dsytrf_();
extern void ssytrf_();
extern void dsytri_();
extern void ssytri_();
extern void dsytrs_();
extern void ssytrs_();
extern void dsyev_();
extern void ssyev_();
extern void dsyevx_();
extern void dgesvd_();
extern void dgesdd_();

struct sorting_string{char *ptr;int index;};
struct sorting_double{double value;int index;};

//these defines are used in nnls.c
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

#include "norm.c"
#include "gamma.c"
#include "nm_gamma.c"
#include "nnls.c"
#include "erfcx.c"

#include "permute.c"
#include "oddsnends.c"
#include "sort.c"
#include "decomp.c"
#include "regfuns.c"

#include "filedata.c"
#include "dataops.c"
#include "filemain.c"
#include "filekins.c"

#include "weights.c"
#include "kinfuns.c"
#include "genefuns.c"

#include "remlmulti.c"
#include "remlfast.c"
#include "blup.c"
#include "he.c"

#include "remllinear.c"
#include "remlgene.c"
#include "remladv.c"

#include "sumfuns.c"

///////////////////////////

int main (int argc, const char * argv[])
{
//this line makes the buffer flush
setlinebuf(stdout);

if(MKL==0){printf("Stop - DOUG HAS FAILED TO COMPILE THIS WITH MKL - PLEASE LET HIM KNOW - POSSIBLE REWARD :)\n\n");}

//set random seed
srand((unsigned)clock()+(unsigned)time(NULL)+(unsigned)getpid());

//set seed for normal generator
zigset(rand());

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\nLDAK - Software for obtaining Linkage Disequilibrium Adjusted Kinships and Loads More\nVersion 5.1 - Help pages at http://www.ldak.org\n-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");

if(sizeof(unsigned char)!=1||sizeof(unsigned short)!=2||sizeof(unsigned)!=4||sizeof(float)!=4||sizeof(size_t)!=8)
{printf("Error, this compilation uses %zd / %zd / %zd / %zd / %zd bytes to save an unsigned char / unsigned short / unsigned integer / float / size_t (LDAK requires 1 / 2 / 4 / 4 / 8)\n", sizeof(unsigned char), sizeof(unsigned short), sizeof(unsigned), sizeof(float), sizeof(size_t));exit(1);}

//declare variables
#include "declare.c"

//deal with command line arguments (come in pairs)
#include "readargs.c"
//this sets mode - see modes.c for a list of modes

//append filenames and check exist
#include "append.c"

//check whether arguments are consistent
#include "consistent.c"

//read section, partition or gene details - will set a few more parameters
if(mode==102||mode==103||mode==104||mode==112||mode==113||mode==137||mode==138||mode==139)
{
#include "readdetails.c"
}

//set a few parameters - in particular those required for parsefiles.c
#include "defaults.c"

//check data/regions/kinships etc
#include "parsefiles.c"

//mode-specific checks
#include "required.c"

//tell the user what will happen
#include "param.c"

///////////////////////////

if(use_data!=5)	//get num_samples_use, and get size of data
{
#include "getnums.c"

if(use_data==1)	//using data (normally)
{
printf("Data contain %d samples and %d predictors", num_samples, num_preds);
if(num_samples_use<num_samples||num_preds_use<num_preds)
{printf("; will be using %d and %d", num_samples_use, num_preds_use);}
printf("\n\n");
}
}
//will deal with use_data=5 (merging) later

///////////////////////////

if(strcmp(topfile,"blank")!=0)	//using tops - get annotations - already have tkeeppreds
{
tchr=malloc(sizeof(double)*num_tops);
tbp=malloc(sizeof(double)*num_tops);
tpreds=malloc(sizeof(char *)*num_tops);
tal1=malloc(sizeof(char)*num_tops);
tal2=malloc(sizeof(char)*num_tops);
for(j=0;j<num_tops;j++)
{
tchr[j]=allchr[tkeeppreds[j]];
tbp[j]=allbp[tkeeppreds[j]];
copy_string(tpreds,j,allpreds[tkeeppreds[j]]);
tal1[j]=allal1[tkeeppreds[j]];
tal2[j]=allal2[tkeeppreds[j]];
}
}	//end of using tops

if(num_regs>0)	//sort regions - must have modes 121/123/124
{
rkeeppreds=malloc(sizeof(int)*num_preds);
regindex=malloc(sizeof(int*)*num_regs);
rnum_preds_use=read_regions(regindex, rkeeppreds, num_regs, regpref, num_preds, allpreds, predorder, num_preds_use, keeppreds);
rpreds=malloc(sizeof(char*)*rnum_preds_use);
ral1=malloc(sizeof(char)*rnum_preds_use);
ral2=malloc(sizeof(char)*rnum_preds_use);
for(j=0;j<rnum_preds_use;j++)
{
copy_string(rpreds,j,allpreds[rkeeppreds[j]]);
ral1[j]=allal1[rkeeppreds[j]];
ral2[j]=allal2[rkeeppreds[j]];
}
}	//end of num_regs>0

///////////////////////////

if(use_data==1)	//reduce predictors to data_length and keeppreds_use
{
#include "setdl.c"
}

////////

if(use_data==1)	//(normal case) sort scalings and pvalues, perhaps reduce, then sort summaries
{
centres=malloc(sizeof(double)*data_length);
mults=malloc(sizeof(double)*data_length);
sqdevs=malloc(sizeof(double)*data_length);
weights=malloc(sizeof(double)*data_length);
pvalues=malloc(sizeof(double)*data_length);

for(j=0;j<data_length;j++){centres[j]=-9999;mults[j]=-9999;sqdevs[j]=-9999;weights[j]=1;pvalues[j]=2;}

if(strcmp(centresfile,"blank")!=0)
{
printf("Reading centres from %s\n", centresfile);
read_centresfile(centresfile, centres, data_length, preds, al1, al2, bimfile);
printf("First few centres are: %s %.4e", preds[0], centres[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4e", preds[j], centres[j]);}}
printf("\n\n");
}

if(strcmp(weightsfile,"blank")!=0)	//read weights from weightsfile
{
printf("Reading weights from %s\n", weightsfile);
read_weightsfile(weightsfile, weights, data_length, preds, keeppreds_use, num_preds, allpreds, bimfile);
printf("First few weights are: %s %.4f", preds[0], weights[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4f", preds[j], weights[j]);}}
printf("\n\n");
}

if(strcmp(pvafile,"blank")!=0)
{
printf("Reading p-values from %s\n", pvafile);
read_pvafile(pvafile, pvalues, data_length, preds);
printf("First few p-values are: %s %.4e", preds[0], pvalues[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4e", preds[j], pvalues[j]);}}
printf("\n\n");
}

if((mode==105&&reduce==1)||mode==112||mode==114||mode==126||mode==127||mode==136||(mode==141&&reduce==1)||mode==151||mode==152||mode==153||mode==158)	//can reduce
{
count=reduce_values(weights,data_length,keeppreds_use,chr,preds,cm,bp,cmbp,al1,al2,centres,pvalues,0,NULL);
if(count==0)
{printf("Error, there are no predictors with non-zero weight\n\n");exit(1);}
if(count<data_length)
{printf("By ignoring those with weight zero, number of predictors reduced to %d", count);
if(mode==105||mode==141){printf(" (to avoid this, use \"--reduce NO\")");}
printf("\n\n");
}
data_length=count;
}

if(strcmp(indhers,"blank")!=0)	//read per-snp heritabilities (and reduce and set her)
{
printf("Reading per-predictor heritabilities from %s\n", indhers);
read_indhers(indhers, weights, data_length, preds);
count=reduce_values(weights,data_length,keeppreds_use,chr,preds,cm,bp,cmbp,al1,al2,centres,pvalues,0,NULL);
if(count==0)
{printf("Error, there are no predictors with non-zero heritability\n\n");exit(1);}
if(count<data_length)
{printf("By ignoring those with non-positive heritability, number of predictors reduced to %d\n", count);}
data_length=count;

if(herscale!=-9999)	//must be mode 158
{
for(j=0;j<data_length;j++){weights[j]*=herscale;}
}

her=0;for(j=0;j<data_length;j++){her+=weights[j];}
printf("Sum of per-predictor heritabilities: %.6f\n", her);
if(her>=1)
{
printf("Warning, this must be less than one; will rescale pre-predictor heritabilities so that they sum to 0.99\n\n");
for(j=0;j<data_length;j++){weights[j]=weights[j]/her*.99;}
her=0.99;
}
printf("\n");
}

if(strcmp(sumsfile,"blank")!=0)	//must be mode 131, 132, 138, 158 or 172
{
nss=malloc(sizeof(double)*data_length);
chis=malloc(sizeof(double)*data_length);
rhos=malloc(sizeof(double)*data_length);

printf("Reading summary statistics from %s\n", sumsfile);
read_sumsfile(sumsfile, nss, chis, rhos, data_length, preds, al1, al2, bimfile, amb, fixn, scaling, 0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis[0], nss[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis[j], nss[j]);}}
printf("\n\n");
}

if(strcmp(sums2file,"blank")!=0)	//must be mode 158
{
nss2=malloc(sizeof(double)*data_length);
chis2=malloc(sizeof(double)*data_length);
rhos2=malloc(sizeof(double)*data_length);

printf("Reading training summary statistics from %s\n", sums2file);
read_sumsfile(sums2file, nss2, chis2, rhos2, data_length, preds, al1, al2, bimfile, amb, fixn, scaling, 0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis2[0], nss2[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis2[j], nss2[j]);}}
printf("\n\n");

if(checksize==1)
{
count=0;for(j=0;j<data_length;j++){count+=(nss2[j]>=nss[j]);}
if(count>0){printf("Error, %d of the predictors have equal or larger sample sizes in %s than in %s; if you are sure %s contains summary statistics for training samples, use \"--check-sizes NO\"\n\n", count, sums2file, sumsfile, sums2file);exit(1);}
}
}
}	//end of use_data=1

///////////////////////////

if(num_resps_use>0)	//get resp (possibly a fake response)
{
resp=malloc(sizeof(double)*num_samples_use*num_resps_use);
respcounts=malloc(sizeof(int)*num_resps_use);

if(strcmp(respfile,"blank")!=0)	//phenotype provided
{read_respfile(respfile, resp, num_samples_use, ids3, num_resps_use, keepresps, num_resps, respcounts, missingvalue, (prev!=-9999||mode==132), pad);}
else	//make a fake phenotype (must have num_resps_use=1, mode 115, 121, 122, 138 or 172)
{
if(mode==121||mode==138)	//need a realistic phenotype (might add more later)
{
if(prev!=-9999)	//binary - will have ascer
{
value=round(ascer*num_samples_use);
for(i=0;i<num_samples_use;i++){resp[i]=(i<value)-value/num_samples_use;}
}
else	//quantitative - use normal quantiles
{
for(i=0;i<num_samples_use;i++){value=(double)(i+1)/(num_samples_use+1);resp[i]=normal_inv(value);}
}
respcounts[0]=num_samples_use;
}
else	//make blank one
{
for(i=0;i<num_samples_use;i++){resp[i]=missingvalue;}
respcounts[0]=0;
}
}	//end of no phenotypes provided
}	//end of using resp

////////

if(num_fixed>0)	//using covar - modes 121/122/123/124/126/127/131/132/133/138/145/151/152/153/169/170/172
{
covar=malloc(sizeof(double)*num_samples_use*num_fixed);
for(i=0;i<num_samples_use;i++){covar[i]=1;}

if(strcmp(covarfile,"blank")!=0)	//normal covariates
{read_covarfile(covarfile, covar+num_samples_use, num_samples_use, ids3, num_covars-1, missingvalue);}

if(strcmp(envfile,"blank")!=0)	//environmental variables
{read_envfile(envfile, covar+num_covars*num_samples_use, num_samples_use, ids3, num_envs, discenv, ids1, ids2);}

if(strcmp(topfile,"blank")!=0)	//have tops
{
tcentres=malloc(sizeof(double)*num_tops);
tvars=malloc(sizeof(double)*num_tops);
read_tops(datafile, covar+(num_covars+num_envs)*num_samples_use, num_samples_use, ids3, num_tops, tkeeppreds, tcentres, tvars, famfile, famhead, dtype, binary, genskip, genheaders, genprobs, num_preds, allpreds, missingvalue, nonsnp, sumsfile);

if(strcmp(sumsfile,"blank")!=0)	//read in summaries for tops - currently redudant, as not allowed
{
tnss=malloc(sizeof(double)*num_tops);
tchis=malloc(sizeof(double)*num_tops);
trhos=malloc(sizeof(double)*num_tops);

printf("Reading summary statistics for top predictors from %s\n", sumsfile);
read_sumsfile(sumsfile, tnss, tchis, trhos, num_tops, tpreds, tal1, tal2, bimfile, amb, fixn, scaling, 0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", tpreds[0], tchis[0], tnss[0]);
for(j=1;j<num_tops;j++){if(j<3){printf(" | %s %.3f %.1f", tpreds[j], tchis[j], tnss[j]);}}
printf("\n\n");
}
}

if(strcmp(cofile,"blank")!=0)	//mode must be 122 or 172
{
printf("Reading coefficients for %d covariates from %s\n\n", num_fixed, cofile);
thetas=malloc(sizeof(double)*num_fixed);
if(countcols(cofile)==1){read_values(cofile,thetas,num_fixed,NULL,1,0,0);}
else{read_values(cofile,thetas,num_fixed,NULL,2,1,0);}
}
}

///////////////////////////

if(num_regs>0)	//deal with regions - sort scalings and summaries, then read data, set rdata_length
{
rcentres=malloc(sizeof(double)*rnum_preds_use);
rmults=malloc(sizeof(double)*rnum_preds_use);
rsqdevs=malloc(sizeof(double)*rnum_preds_use);
rweights=malloc(sizeof(double)*rnum_preds_use);

for(j=0;j<rnum_preds_use;j++){rcentres[j]=-9999;rmults[j]=-9999;rweights[j]=1;}

if(strcmp(weightsfile,"blank")!=0)	//read weights from weightsfile and reduce
{
printf("Reading weights for region predictors from %s\n", weightsfile);
read_weightsfile(weightsfile, rweights, rnum_preds_use, rpreds, rkeeppreds, num_preds, allpreds, bimfile);
printf("First few weights are: %s %.6f", rpreds[0], rweights[0]);
for(j=1;j<rnum_preds_use;j++){if(j<5){printf(" | %s %.6f", rpreds[j], rweights[j]);}}
printf("\n\n");

count=reduce_values(rweights,rnum_preds_use,rkeeppreds,NULL,rpreds,NULL,NULL,NULL,ral1,ral2, NULL, NULL, num_regs, regindex);
//reduce_values will report an error if count==0
if(count<rnum_preds_use)
{printf("By ignoring those with weight zero, the number of (unique) region predictors is reduced to %d\n\n", count);}
rnum_preds_use=count;
}

if(strcmp(sumsfile,"blank")!=0)	//read in summaries for regions (must all be present)
{
rnss=malloc(sizeof(double)*rnum_preds_use);
rchis=malloc(sizeof(double)*rnum_preds_use);
rrhos=malloc(sizeof(double)*rnum_preds_use);

printf("Reading summary statistics for %d region predictors from %s\n", rnum_preds_use, sumsfile);
read_sumsfile(sumsfile, rnss, rchis, rrhos, rnum_preds_use, rpreds, ral1, ral2, bimfile, amb, fixn, scaling, 0, -9999);

printf("First few stats and ns are: %s %.3f %.1f", rpreds[0], rchis[0], rnss[0]);
for(j=1;j<rnum_preds_use;j++){if(j<3){printf(" | %s %.3f %.1f", rpreds[j], rchis[j], rnss[j]);}}
printf("\n\n");
}

//read in regions data and standardize
rdata_warn(rnum_preds_use,num_samples_use);
rdata=malloc(sizeof(double)*num_samples_use*rnum_preds_use);

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
(void)read_data_fly(datafile, dtype, rdata, NULL, num_samples_use, keepsamps, 0, rnum_preds_use, rkeeppreds, datainputgz, 0, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
if(binary==0){gzclose(datainputgz);}

stand_data(rdata, rcentres, rmults, rsqdevs, num_samples_use, rnum_preds_use, missingvalue, power, 0, hwestand, rweights, 1, rpreds);

//adjust regindex in order to remove trivial predictors, and possibly prune 
rdata_length=prune_regions(num_regs, regindex, rmults, rprune, rdata, num_samples_use);

if(mode==121&&shortcut==1)	//warn if too many region predictors
{
if(rdata_length>num_samples_use)
{printf("Warning, there are more region predictors (%d) than samples (%d); usually it is better to include them as kinships instead\n\n", rdata_length, num_samples_use);}
if(rdata_length>2000&&rdata_length<num_samples_use)
{printf("Warning, the number of region predictors is quite high (%d); consider using \"--region-prune\" to thin these or including them as kinships instead\n\n", rdata_length);}
}
}	//end of num_regs>0

///////////////////////////

if(num_kins>0&&(mode==120||mode==121||mode==123||mode==124||mode==131||mode==133||mode==138))	//deal with kinships - will store in mkins or mkins_single
{
if(memsave==0)	//for calc-sim-grms, he and pcgc don't need double precision
{
if(mode==120||mode==123||mode==124){kin_warn(num_kins, num_samples_use, 0, 1);}
else{kin_warn(num_kins, num_samples_use, 0, 0);}
}

if(mode==120||mode==123||mode==124){mkins_single=malloc(sizeof(float *)*num_kins);}
else{mkins=malloc(sizeof(double *)*num_kins);}

kintraces=malloc(sizeof(double)*num_kins);
for(k=0;k<num_kins;k++)
{
if(memsave==0)	//read kinships (and get traces direct)
{
if(mode==120||mode==123||mode==124)
{
mkins_single[k]=malloc(sizeof(float*)*num_samples_use*num_samples_use);
read_kins(kinstems[k], NULL, mkins_single[k], 1.0, num_samples_use, ids3, 0);
sum=0;for(i=0;i<num_samples_use;i++){sum+=mkins_single[k][(size_t)i*num_samples_use+i];}
kintraces[k]=sum/num_samples_use;
}
else
{
mkins[k]=malloc(sizeof(double*)*num_samples_use*num_samples_use);
read_kins(kinstems[k], mkins[k], NULL, 1.0, num_samples_use, ids3, 0);
sum=0;for(i=0;i<num_samples_use;i++){sum+=mkins[k][(size_t)i*num_samples_use+i];}
kintraces[k]=sum/num_samples_use;
}
}
else	//just get traces
{kintraces[k]=read_kin_trace(kinstems[k], num_samples_use, ids3, 0);}

if(fabs(kintraces[k]-1)>.01)
{printf("Warning, Kinship %d has average diagonal value %.6f (usually this is 1)\n", k+1, kintraces[k]);}
}
printf("\n");
}

////////

if(num_kins==1&&((mode==121&&shortcut==1)||mode==131||mode==138))	//will be using a decomposition
{
decomp_warn(num_samples_use);
U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);

if(strcmp(eigenfile,"blank")==0)	//decompose
{
if(memsave==0)	//are there cases when we decompose and must keep kinship? - I think not
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){U[(size_t)i2*num_samples_use+i]=mkins[0][(size_t)i2*num_samples_use+i];}
}
}
else
{
kin_warn(1, num_samples_use, 0, 0);
read_kins(kinstems[0], U, NULL, 1.0, num_samples_use, ids3, 0);
printf("\n");
}

printf("Performing an eigen-decomposition; this can take a while\nIf preferred, you can perform the analysis in two steps; first replace the main argument with \"--decompose\" (this will save the eigen-decomposition) then restart adding \"--eigen\"\n\n");

lwork=-1;
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);
}
else	//read it from file
{read_eigens(eigenfile, U, E, num_samples_use, ids3, -9999);}
}

///////////////////////////
//start of work
///////////////////////////

if(mode==101)	//cut-weights
{
if(nothin==0)	//reduce data_length to thinned SNPs 
{
//usedpreds starts as -1 (can be considered) or -2 (filtered out)
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=-2;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds_use[j]]=-1;}

if(window_kb>25||window_length>25)	//do in two steps
{
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("Will prune in two passes; first with windows of size 0.01cM, then %.4fcM\n\n", window_cm);}
else
{printf("Will prune in two passes; first with windows of size 10kb, then %.2fkb\n\n", window_kb);}
}
else
{printf("Will prune in two passes; first with windows of length 10, then %d\n\n", window_length);}
bitsize2=bitsize;window_kb2=window_kb;window_length2=window_length;
if(window_kb!=-9999){window_kb=10;}
if(window_length!=-9999){window_length=10;}
flag=1;
#include "prune.c"
bitsize=bitsize2;window_kb=window_kb2;window_length=window_length2;flag=2;
flag=2;
#include "prune.c"
}
else	//do in one
{
flag=0;
#include "prune.c"
}

free(usedpreds);
}
if(nothin==1)	//probably pointless, but print "duplicates"
{
sprintf(filename,"%sthin.dups",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor Replacement\n"); 
for(j=0;j<num_preds_use;j++){fprintf(output, "%s %s\n", allpreds[keeppreds[j]], allpreds[keeppreds[j]]);}
fclose(output);
}
//for nothin=2, nothing to do

cut_sections(data_length, chr, cmbp, section_kb, section_length, buffer_kb, buffer_length, num_samples_use, folder, window_kb, window_length, window_cm, datafile, bimfile, extract, nothin);
}	//end of mode=101

if(mode==102||mode==104)	//calc-weights
{
#include "solveweights.c"
}

if(mode==103||mode==104)	//join-weights
{
join_sections(folder, num_sections, sstarts, sends, sstarts2, sends2, allpreds, num_preds, keeppreds, num_preds_use, keeppreds_use, spread);
}

if(mode==105)	//adjust-weights
{
#include "tagging.c"
}

////////

if(mode==106||mode==107)	//thin or thin-tops
{
//usedpreds starts as -1 (can be considered) or -2 (filtered out)
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=-2;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds_use[j]]=-1;}

if(window_kb>25||window_length>25)	//do in two steps
{
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("Will prune in two passes; first with windows of size 0.01cM, then %.4fcM\n\n", window_cm);}
else
{printf("Will prune in two passes; first with windows of size 10kb, then %.2fkb\n\n", window_kb);}
}
else
{printf("Will prune in two passes; first with windows of length 10, then %d\n\n", window_length);}
bitsize2=bitsize;window_kb2=window_kb;window_length2=window_length;
if(window_kb!=-9999){window_kb=10;}
if(window_length!=-9999){window_length=10;}
flag=1;
#include "prune.c"
bitsize=bitsize2;window_kb=window_kb2;window_length=window_length2;flag=2;
flag=2;
#include "prune.c"
}
else
{
flag=0;
#include "prune.c"
}

free(usedpreds);
}	//end of mode=106/107

if(mode==108||mode==109)	//find-tags or remove-tags - each target SNP becomes a gene
{
#include "finding.c"
}

///////////////////////////

if(mode==111)	//cut-kins - make kinships partitions
{
cut_partitions(data_length, chr, preds, part_length, bychr, num_parts, partpref, checkpart, folder, datafile, bimfile, extract);
}

if(mode==112||mode==114)	//calculate kinships
{
#include "kins.c"
}

if(mode==113)	//join-kins - have set kinstems (in consistent.c)
{
sprintf(outfile,"%skinships.all",folder);
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 0, NULL);
}

////////

if(mode==115)	//filter
{
#include "filter.c"
}

if(mode==116)	//add-grm
{
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 1, NULL);
}

if(mode==117)	//sub-grm
{
if(extract==1)
{
#include "subkin.c"
}	//end of subtracting predictors
else
{
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 2, NULL);
}
}

if(mode==118)	//convert-gz
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kinsgz(kinstems[0], kins_single, num_samples_use, ids3);
write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 1);
free(kins_single);
}

if(mode==119)	//convert-raw
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kinsraw(kinstems[0], kins_single, num_samples_use, ids3);
write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 1);
free(kins_single);
}

////////

if(mode==120)	//computing correlations between kinships - must have single kins and memsave=0
{
cors=malloc(sizeof(double)*num_kins*num_kins);

for(k=0;k<num_kins;k++)
{
cors[k+k*num_kins]=1;
for(k2=k+1;k2<num_kins;k2++)
{
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;count=0;
for(i2=1;i2<num_samples_use;i2++)
{
for(i=0;i<i2;i++)
{
sum+=mkins_single[k][(size_t)i2*num_samples_use+i];
sum2+=mkins_single[k2][(size_t)i2*num_samples_use+i];
sumsq+=pow(mkins_single[k][(size_t)i2*num_samples_use+i],2);
sumsq2+=pow(mkins_single[k2][(size_t)i2*num_samples_use+i],2);
sumsq3+=mkins_single[k][(size_t)i2*num_samples_use+i]*mkins_single[k2][(size_t)i2*num_samples_use+i];
count++;
}
}
mean=sum/count;mean2=sum2/count;
value=(sumsq3-count*mean*mean2)/pow(sumsq-count*mean*mean,.5)/pow(sumsq2-count*mean2*mean2,.5);
cors[k+k2*num_kins]=value;cors[k2+k*num_kins]=value;
printf("Correlation between Matrices %d and %d is %.6f\n", k+1, k2+1, value);
}	//end of k2 loop
}
printf("\n");

sprintf(filename,"%s.cors.grms",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(k=0;k<num_kins;k++)
{
for(k2=0;k2<num_kins;k2++){fprintf(output,"%.10f ", cors[k+k2*num_kins]);}
fprintf(output,"\n");
}
fclose(output);

printf("Correlations saved in %s\n\n", filename);

free(cors);
}	//end of mode=120

///////////////////////////

if(mode==121)	//reml
{
//allocate - begin by guessing memory used by allocations within remlmulti.c
if(shortcut==0){anal_warn((3+memsave)*num_samples_use, num_samples_use);}

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

if(num_regs>0)	//will be using X
{
X=malloc(sizeof(double)*num_samples_use*rdata_length);
Xstarts=malloc(sizeof(int)*num_regs);
Xends=malloc(sizeof(int)*num_regs);
Xsums=malloc(sizeof(double)*num_regs);
if(strcmp(sumsfile,"blank")!=0)
{
Xnss=malloc(sizeof(double)*rdata_length);
Xrhos=malloc(sizeof(double)*rdata_length);
}
}

for(m=0;m<num_resps_use;m++)	//loop through phenotypes
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Z (maybe permuted) and maybe X (not permuted) - will fill Y (maybe permuted) below
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}
if(strcmp(sumsfile,"blank")==0)
{Xtotal=fill_X(X, Xstarts, Xends, Xsums, NULL, NULL, num_samples_use, num_regs, regindex, rdata, NULL, NULL);}
else
{Xtotal=fill_X(X, Xstarts, Xends, Xsums, Xnss, Xrhos, num_samples_use, num_regs, regindex, rdata, rnss, rrhos);}

if(mpheno!=-1)	//single phenotype
{
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]];}
strcpy(filename,outfile);
}
else	//analysing phenotypes in turn
{
printf("Analysing Response %d of %d\n", m+1, num_resps_use);
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]+m*num_samples_use];}
sprintf(filename, "%s.%d", outfile, m+1);
}

if(num_kins==0&&num_regs>0&&shortcut==1)	//do quick version - might have summary stats
{
if(strcmp(sumsfile,"blank")==0){fast_reml(num_samples_use, num_covars, num_tops, num_envs, num_regs, Y, Z, X, Xtotal, Xstarts, Xends, Xsums, NULL, NULL, prev, respcounts[m], compact, constrain, shrink, strip, tol, maxiter, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, rkeeppreds, regindex, rcentres, rmults, rweights, tkeeppreds, tcentres);}
else{fast_reml(num_samples_use, num_covars, 0, 0, num_regs, Y, Z, X, Xtotal, Xstarts, Xends, Xsums, Xnss, Xrhos, prev, -9999, compact, constrain, shrink, strip, tol, maxiter, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, rkeeppreds, regindex, rcentres, rmults, rweights, tkeeppreds, tcentres);}
}
else
{multi_reml(num_samples_use, num_covars, num_envs, num_tops, num_kins, num_regs, Y, Z, mkins, kintraces, kinsums, X, Xtotal, Xstarts, Xends, Xsums, prev, respcounts[m], shortcut, U, E, compact, discenv, liab, constrain, tol, maxiter, memsave, kinstems, ids3, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, rkeeppreds, regindex, rcentres, rmults, rweights, tkeeppreds, tcentres, NULL);}
}	//end of m loop

free(order);free(Y);free(Z);
if(num_regs>0){free(X);free(Xstarts);free(Xends);free(Xsums);}
if(num_regs>0&&strcmp(sumsfile,"blank")!=0){free(Xnss);free(Xrhos);}
}

////////

if(mode==122)	//calc-blups
{
#include "timesa.c"
}

////////

if(mode==123||mode==124)	//he and pcgc regression
{
//allocate
order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

if(num_regs>0)	//will be using X
{
X=malloc(sizeof(double)*num_samples_use*rdata_length);
Xstarts=malloc(sizeof(int)*num_regs);
Xends=malloc(sizeof(int)*num_regs);
Xsums=malloc(sizeof(double)*num_regs);
}

if(num_blocks==-1||num_blocks>num_samples_use){num_blocks=num_samples_use;}

for(m=0;m<num_resps_use;m++)	//loop through phenotypes
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Z (maybe permuted) and maybe X (not permuted) - will fill Y (maybe permuted) below
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}
Xtotal=fill_X(X, Xstarts, Xends, Xsums, NULL, NULL, num_samples_use, num_regs, regindex, rdata, NULL, NULL);

if(mpheno!=-1)	//single phenotype
{
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]];}
strcpy(filename,outfile);
}
else	//analysing phenotypes in turn
{
printf("Analysing Response %d of %d\n", m+1, num_resps_use);
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]+m*num_samples_use];}
sprintf(filename, "%s.%d", outfile, m+1);
}

he_reg(num_samples_use, num_subs, num_covars, num_envs, num_tops, num_kins, num_regs, subindex, Y, Z, mkins_single, kintraces, kinsums, X, Xtotal, Xstarts, Xends, Xsums, prev, respcounts[m], discenv, num_blocks, tol, memsave, kinstems, ids3, filename, (mode==124));
}	//end of m loop

free(order);free(Y);free(Z);
if(num_regs>0){free(X);free(Xstarts);free(Xends);free(Xsums);}
}

////////

if(mode==126||mode==127)	//fast he and pcgc regression
{
#include "fast.c"
}

///////////////////////////

if(mode==131||mode==132)	//linear (num_kins will be 0 or 1) and logistic
{
#include "single.c"
}

////////

if(mode==133)	//solve-null - will always have shortcut=0 and num_resps_use=1
{
anal_warn((3+memsave)*num_samples_use, num_samples_use);

Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

vstarts=malloc(sizeof(double)*num_kins);

//fill Y and Z
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[i];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//estimate the variances
printf("Estimating the genetic variances\n");
multi_reml(num_samples_use, num_covars, 0, num_tops, num_kins, 0, Y, Z, mkins, kintraces, kinsums, NULL, 0, NULL, NULL, NULL, prev, respcounts[0], 0, NULL, NULL, 0, 0, 0, constrain, tol, maxiter, memsave, kinstems, ids3, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, NULL, NULL, NULL, NULL, NULL, tkeeppreds, tcentres, vstarts);

//add up kinships and save (for simplicity, will read kinships from file)
printf("Creating the merged kinship matrix\n\n");
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 3, vstarts);
printf("Merged kinship matrix saved with stem %s; complete the mixed-model linear regression by using \"--linear\" with \"--grm %s\" (with the same covariates, top predictors and sample filterings used now)\n\n", outfile, outfile);

free(Y);free(Z);
free(vstarts);
}	//end of solve-null

////////

if(mode==136)	//cut-genes
{
//get an upper limit on number of genes/chunks
if(strcmp(genefile,"blank")!=0)	//using genefile - easy
{num_genes=countrows(genefile)-check_head(genefile,"Gene","Name",1);}

if(chunks!=-9999)	//using weights
{
sum=weights[0];count=1;
for(j=1;j<data_length;j++)
{
sum+=weights[j];
if(chr[j]!=chr[j-1]){count++;}
}
num_genes=sum/chunks*(1+overlap)+2*count;
}

if(chunksbp==1)	//each basepair observed becomes a gene - easy
{num_genes=data_length;}

if(chunksbp>1)
{
sum=0;count=1;
for(j=1;j<data_length;j++)
{
if(chr[j]!=chr[j-1]){sum+=bp[j-1];count++;}
}
sum+=bp[data_length-1];

if(sum/chunksbp*(1+overlap)+2*count>pow(2,31))
{printf("Error, there are too many chunks; to continue, you must increase chunksbp (currently %d)\n\n", chunksbp);exit(1);}
num_genes=sum/chunksbp*(1+overlap)+2*count;
}

gchr=malloc(sizeof(int)*num_genes);
gbp1=malloc(sizeof(double)*num_genes);
gbp2=malloc(sizeof(double)*num_genes);
gnames=malloc(sizeof(char*)*num_genes);
gstrand=malloc(sizeof(int)*num_genes);
gstarts=malloc(sizeof(int)*num_genes);
gends=malloc(sizeof(int)*num_genes);
for(g=0;g<num_genes;g++)
{
if(chunks!=-9999||chunksbp!=-9999){gnames[g]=malloc(sizeof(char)*100);}
gstarts[g]=-9999;gends[g]=-9999;
}

cut_genes(gchr, gbp1, gbp2, gnames, gstrand, gstarts, gends, data_length, chr, preds, bp, weights, genefile, chunks, chunksbp, up_buffer, down_buffer, minweight, overlap, pvafile, pvalues, part_length, bychr, 0, folder, datafile, bimfile, extract, -9999);

for(g=0;g<num_genes;g++){free(gnames[g]);}free(gnames);
free(gchr);free(gbp1);free(gbp2);free(gstrand);free(gstarts);free(gends);
}

////////

if(mode==137||mode==138)	//calc-genes-kins and calc-genes-reml
{
#include "genes.c"
}

if(mode==139)	//join-genes-reml
{
join_genes_reml(folder, num_genes, gnames, gstarts, gends, gparts, gchr, gbp1, gbp2, gamp, gam1, gam2, cut1, cut2, data_length, preds, keeppreds_use);
}

///////////////////////////

if(mode==141)	//calc-tagging
{
#include "tagging.c"
}

if(mode==142)	//join-tagging
{
#include "join.c"
}

if(mode==143)	//merge-tagging
{
#include "combine.c"
}

if(mode==144)	//reduce-tagging
{
#include "reduce.c"
}

////////

if(mode==146)	//sum-hers
{
#include "sumsa.c"
}

if(mode==147)	//sum-cors
{
#include "sumsb.c"
}

if((mode==146&&strcmp(matfile,"blank")!=0)||mode==148||mode==149)	//calc-ind-hers and calc-exps
{
#include "exps.c"
}

if(mode==150)	//calc-posts
{
#include "posts.c"
}

///////////////////////////

if(mode==151||mode==152||mode==153)	//blup, bolt and bayesr
{
#include "boltbayesr.c"
}

if(mode==156)	//calc-cors
{
#include "cors.c"
}

if(mode==157)	//join-cors
{
#include "concat.c"
}

if(mode==158)	//ridge-lite, bolt-lite, bayesr-lite and lasso-lite
{
#include "lites.c"
}

if(mode==159)	//pseudo-summaries
{
#include "timesd.c"
}

///////////////////////////

if(mode==161)	//pca
{
if(strcmp(eigenfile,"blank")==0)
{
kin_warn(1, num_samples_use, 0, 0);
kins=malloc(sizeof(double*)*num_samples_use*num_samples_use);
}

U=malloc(sizeof(double)*num_samples_use*axes);
E=malloc(sizeof(double)*axes);

if(strcmp(eigenfile,"blank")==0)	//must compute pcs
{
read_kins(kinstems[0], kins, NULL, 1.0, num_samples_use, ids3, 0);
printf("\n");

printf("Performing the decomposition; this can take a while\n\n");

vl=-1;vu=-1;	//vl and vu required but not used
count=num_samples_use-axes+1;
iwork=malloc(sizeof(int)*5*num_samples_use);
ifail=malloc(sizeof(int)*num_samples_use);
lwork=-1;
dsyevx_("V", "I", "U", &num_samples_use, kins, &num_samples_use, &vl, &vu, &count, &num_samples_use, &tol, &found, E, U, &num_samples_use, &wkopt, &lwork, iwork, ifail, &info);
if(info!=0){printf("PCA error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyevx_("V", "I", "U", &num_samples_use, kins, &num_samples_use, &vl, &vu, &count, &num_samples_use, &tol, &found, E, U, &num_samples_use, work, &lwork, iwork, ifail, &info);
if(info!=0){printf("PCA error 2; please tell Doug\n\n");exit(1);}

free(kins);
free(iwork);free(ifail);free(work);
}
else	//can get from eigen-decomposition
{read_eigens(eigenfile, U, E, num_samples_use, ids3, axes);}

sprintf(filename, "%s.vect", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output, "%s %s ", ids1[i], ids2[i]);
for(k=0;k<axes;k++){fprintf(output, "%.6f ", U[(size_t)(axes-1-k)*num_samples_use+i]);}
fprintf(output, "\n");
}
fclose(output);

sprintf(filename2, "%s.values", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(k=0;k<axes;k++){fprintf(output2, "%.6f\n", E[axes-1-k]);}
fclose(output2);

sprintf(filename3, "%s.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
fprintf(output3, "Kinship %s\n", kinstems[0]);
fclose(output3);

printf("First %d eigenvectors saved in %s, with eigenvalues saved in %s\n\n", axes, filename, filename2);

free(U);free(E);
}	//end of mode=161

////////

if(mode==162)	//calc-blups
{
#include "timesb.c"
}

if(mode==163)	//decompose and save
{
U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);

kin_warn(1, num_samples_use, 0, 1);
read_kins(kinstems[0], U, NULL, 1.0, num_samples_use, ids3, 0);
printf("\n");

printf("Performing the decomposition; this can take a while\n\n");

lwork=-1;
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug\n\n");exit(1);}

write_eigen(outfile, U, E, num_samples_use, ids1, ids2, kinstems[0], eigenraw);

free(U);free(E);free(work);
}	//end of mode=163

if(mode==164)	//adjust-grm
{
kin_warn(1, num_samples_use, 0, 0);
kins=malloc(sizeof(double*)*num_samples_use*num_samples_use);

Z=malloc(sizeof(double)*num_samples_use*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);

kinZ=malloc(sizeof(double)*num_samples_use*num_fixed);
ZTkinZ=malloc(sizeof(double)*num_fixed*num_fixed);
W=malloc(sizeof(double)*num_samples_use*num_fixed);
WZTkinZ=malloc(sizeof(double)*num_samples_use*num_fixed);

read_kins(kinstems[0], kins, NULL, 1.0, num_samples_use, ids3, 0);
printf("\n");

//fill Z
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//get inverseZTZ
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &num_samples_use, &alpha, Z, &num_samples_use, Z, &num_samples_use, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

//then kinZ, ZTkinZ, W=ZinvZTZ and WZTkinZ
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_samples_use, &num_fixed, &num_samples_use, &alpha, kins, &num_samples_use, Z, &num_samples_use, &beta, kinZ, &num_samples_use);
dgemm_("T", "N", &num_fixed, &num_fixed, &num_samples_use, &alpha, Z, &num_samples_use, kinZ, &num_samples_use, &beta, ZTkinZ, &num_fixed);
dgemm_("N", "N", &num_samples_use, &num_fixed, &num_fixed, &alpha, Z, &num_samples_use, ZTZ, &num_fixed, &beta, W, &num_samples_use);
dgemm_("N", "N", &num_samples_use, &num_fixed, &num_fixed, &alpha, W, &num_samples_use, ZTkinZ, &num_fixed, &beta, WZTkinZ, &num_samples_use);

//ready to adjust kins
printf("Regressing the kinship matrix on the fixed effects\n\n");

//first subtract WZTkin and kinZWT
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &num_fixed, &alpha, W, &num_samples_use, kinZ, &num_samples_use, &beta, kins, &num_samples_use);
dgemm_("N", "T", &num_samples_use, &num_samples_use, &num_fixed, &alpha, kinZ, &num_samples_use, W, &num_samples_use, &beta, kins, &num_samples_use);

//then add WZTkinZ WT
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &num_fixed, &alpha, WZTkinZ, &num_samples_use, W, &num_samples_use, &beta, kins, &num_samples_use);

write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 1);

if(kindetails==1)	//copy over details and adjust
{
sprintf(cmd, "cp %s.grm.details %s.grm.details", kinstems[0], outfile);
system(cmd);
sprintf(cmd, "cp %s.grm.adjust %s.grm.adjust", kinstems[0], outfile);
system(cmd);
}

//make root
sprintf(filename,"%s.grm.root",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
fprintf(output,"Kinship %s\n", kinstems[0]);
if(strcmp(covarfile,"blank")!=0){fprintf(output,"Covariates %s\n", covarfile);}
else{fprintf(output,"Covariates none\n");}
if(strcmp(topfile,"blank")!=0){fprintf(output,"Top_Predictors %s\n", topfile);}
else{fprintf(output,"Top_Predictors none\n");}
if(strcmp(envfile,"blank")!=0){fprintf(output,"Environments %s\n", envfile);}
else{fprintf(output,"Environments none\n");}
fclose(output);

free(kins);
free(Z);free(ZTZ);free(ZTZ2);free(ZTZ3);
free(kinZ);free(ZTkinZ);free(W);free(WZTkinZ);
}	//end of mode=164

////////

if(mode==166||mode==167||mode==168||mode==169||mode==170)	//truncate-grm, pca-grm, square-grm and gxemms
{
#include "others.c"
}
///////////////////////////

if(mode==171)	//calc-stats
{
#include "stats.c"
}

if(mode==172)	//calc-blups
{
#include "timesc.c"
}

if(mode==173)	//make-phenos
{
#include "phens.c"
}

if(mode==174)	//make-snps
{
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.bed", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
gen=108;fwrite(&gen, sizeof(unsigned char), 1, output2);
gen=27;fwrite(&gen, sizeof(unsigned char), 1, output2);
gen=1;fwrite(&gen, sizeof(unsigned char), 1, output2);

bittotal=(num_snps-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_snps){bitend=num_snps;}
bitlength=bitend-bitstart;

printf("Making SNPs for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Making SNPs for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

for(j=0;j<bitlength;j++)
{
gen=0;
maf=(double)rand()/RAND_MAX*(maf2-maf1)+maf1;
for(i=0;i<num_inds;i++)
{
unifrand=(double)rand()/RAND_MAX;
gen2=0;					//genotype 2
if(unifrand<1-pow(maf,2)){gen2=2;}	//genotype 1
if(unifrand<pow(1-maf,2)){gen2=3;}	//genotype 0
gen+=gen2<<(2*(i%4));
if(i%4==3||i==num_inds-1){fwrite(&gen, sizeof(unsigned char), 1, output2);gen=0;}
}
}
}	//end of bit loop
printf("\n");
fclose(output2);

sprintf(filename3, "%s.bim", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
for(j=0;j<num_snps;j++){fprintf(output3,"1 SNP%d 0 %d A B\n", j+1, 1000*(j+1));}
fclose(output3);

sprintf(filename4, "%s.fam", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
for(i=0;i<num_inds;i++){fprintf(output4,"IND%d IND%d 0 0 0 0\n", i+1, i+1);}
fclose(output4);

printf("Simulated data saved in %s, %s and %s\n\n", filename2, filename3, filename4);
}	//end of mode=174

if(mode==176)	//jackknife
{
#include "jackknife.c"
}

if(mode==177)	//cut-folds
{
usedids=malloc(sizeof(int)*num_samples_use);
order=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){order[i]=i;}
permute_int(order,num_samples_use);

count=(num_samples_use-1)/num_folds+1;
printf("Dividing the %d samples into %d folds with size (approximately) %d\n\n", num_samples_use, num_folds, count);

for(k=0;k<num_folds;k++)
{
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(i=0;i<num_samples_use;i++)
{
if(i>=k*count&&i<(k+1)*count){usedids[order[i]]=1;}
}

sprintf(filename,"%s.train%d", outfile, k+1);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
sprintf(filename2,"%s.test%d", outfile, k+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==1){fprintf(output2,"%s %s\n", ids1[i], ids2[i]);}
else{fprintf(output,"%s %s\n", ids1[i], ids2[i]);}
}
fclose(output);
fclose(output2);
}	//end of k loop

printf("Training folds saved in %s.train1, ..., %s.train%d;  testing folds in %s.test1, .., %s.test%d\n\n", outfile, outfile, num_folds, outfile, outfile, num_folds);

free(usedids);free(order);
}	//end of mode=177

if(mode==178)	//find-gaussian
{
#include "grid.c"
}

///////////////////////////

if(mode==181||mode==182||mode==183||mode==184||mode==185)	//make-data
{
#include "merge.c"
}

if(mode==186||mode==187||mode==188||mode==189)	//condense data
{
#include "condense.c"
}	//end of condense data

if(mode==190)	//calc-sim-data
{
#include "merge.c"
}

///////////////////////////
//free stuff
///////////////////////////

//allocations from readdetails.c

if(mode==102||mode==103||mode==104){free(sstarts);free(sends);free(sstarts2);free(sends2);}
if(mode==112){free(pstarts);free(pends);}
if(mode==137||mode==138||mode==139)
{
for(g=0;g<num_genes;g++){free(gnames[g]);}free(gnames);
free(gstarts);free(gends);free(gparts);free(gchr);free(gbp1);free(gbp2);free(gpvas);
}

//datafile, kinships and resp allocations (from parsefiles.c and top of ldak.c)

if(use_data!=0)
{
for(k=0;k<num_files;k++){free(datastems[k]);free(bimstems[k]);free(famstems[k]);}
free(datastems);free(bimstems);free(famstems);
}

if(num_kins>0)
{
for(k=0;k<num_kins;k++){free(kinstems[k]);}free(kinstems);
free(kinnums);free(kinsums);

if(mode==120||mode==123||mode==124)
{for(k=0;k<num_kins;k++){if(memsave==0){free(mkins_single[k]);}}
free(mkins_single);free(kintraces);}
if(mode==121||mode==131||mode==133||mode==138)
{for(k=0;k<num_kins;k++){if(memsave==0){free(mkins[k]);}}free(mkins);free(kintraces);}
}

if(num_resps_use>0){free(keepresps);free(resp);free(respcounts);}

//more allocations from parsefiles.c

if(strcmp(taglist,"blank")!=0)
{
for(k=0;k<num_tags;k++){free(tagstems[k]);}
free(tagstems);
}

if(strcmp(matlist,"blank")!=0)
{
for(k=0;k<num_tags;k++){free(matstems[k]);}
free(matstems);
}

if(strcmp(catfile,"blank")!=0){free(keepparts);free(keepparts2);}

if(strcmp(corslist,"blank")!=0)
{
for(k=0;k<num_cors;k++){free(corstems[k]);}
free(corstems);
}

//allocations from getnums.c 

if(use_data!=5&&num_samples!=-9999)
{
for(i=0;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);free(ids3[i]);}
free(ids1);free(ids2);free(ids3);

if(num_subs>0)
{
for(s=0;s<num_subs;s++){free(subindex[s]);}free(subindex);
}

if(use_data==1||use_data==2||use_data==4){free(keepsamps);}
}

if(use_data==1||use_data==2||use_data==3)
{
free(allchr);free(allcm);free(allbp);free(allal1);free(allal2);
for(j=0;j<num_preds;j++){free(allpreds[j]);}free(allpreds);
free(predorder);free(keeppreds);
}

//global allocations from setdl.c - will have done mode-specific ones above

if(use_data==1)
{
free(keeppreds_use);
free(chr);free(cm);free(bp);free(cmbp);free(al1);free(al2);
for(j=0;j<data_length;j++){free(preds[j]);}free(preds);
}

//tops and regions (from top of ldak.c)

if(strcmp(topfile,"blank")!=0)
{
free(tkeeppreds);
free(tchr);free(tbp);free(tal1);free(tal2);
for(j=0;j<num_tops;j++){free(tpreds[j]);}free(tpreds);
free(tcentres);free(tvars);
if(strcmp(sumsfile,"blank")!=0){free(tnss);free(tchis);free(trhos);}
}

if(num_regs>0)
{
free(rkeeppreds);
for(r=0;r<num_regs;r++){free(regindex[r]);}free(regindex);
free(ral1);free(ral2);
for(j=0;j<rnum_preds_use;j++){free(rpreds[j]);}free(rpreds);
free(rcentres);free(rmults);free(rsqdevs);free(rweights);
if(strcmp(sumsfile,"blank")!=0){free(rnss);free(rchis);free(rrhos);}
free(rdata);
}

//more allocations from top of ldak.c

if(use_data==1)
{
free(centres);free(mults);free(sqdevs);free(weights);free(pvalues);
if(strcmp(sumsfile,"blank")!=0){free(nss);free(chis);free(rhos);}
if(strcmp(sums2file,"blank")!=0){free(nss2);free(chis2);free(rhos2);}
}

if(num_fixed>0)
{
free(covar);
if(strcmp(cofile,"blank")!=0){free(thetas);}
}

if(num_kins==1&&((mode==121&&shortcut==1)||mode==131||mode==138)){free(U);free(E);}

///////////////////////////

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\nMission completed. All your basepair are belong to us :)\n-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n");

return(0);
}	//end of main

///////////////////////////

//old code from mode 146
/*
if(strcmp(topfile,"blank")!=0)	//check predictors, then get topher
{
wcount=0;
for(j=0;j<num_tops;j++)
{
for(j2=0;j2<count2;j2++)
{
if(strcmp(tpreds[j],spreds[j2])==0)
{printf("Error, %s is both a top predictor and in the tagging file\n", tpreds[j]);wcount++;break;}
}}
if(wcount>0){printf("\n");exit(1);}

ZTZ=malloc(sizeof(double)*num_tops*num_tops);
ZTZ2=malloc(sizeof(double)*num_tops);
ZTZ3=malloc(sizeof(double)*num_tops*num_tops);
ZTY=malloc(sizeof(double)*num_tops);

//get ZTZ, ZTY/YTY^.5 (stored in ZTY) and inv ZTZ - can ignore intercept
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_tops, &num_tops, &num_samples_use, &alpha, covar+num_samples_use, &num_samples_use, covar+num_samples_use, &num_samples_use, &beta, ZTZ, &num_tops);
for(j=0;j<num_tops;j++){ZTY[j]=trhos[j]*pow(ZTZ[j+j*num_tops],.5);}
(void)eigen_invert(ZTZ, num_tops, ZTZ2, -1, ZTZ3, 1);

//toppher is ZTY ZTZ ZTY
topher=0;
for(j=0;j<num_tops;j++)
{
for(j2=0;j2<num_tops;j2++){topher+=ZTY[j]*ZTZ[j+j2*num_tops]*ZTY[j2];}
}
free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTY);

if(num_tops==1){printf("The proportion of variance explained by the top predictor is %.4f\n\n", topher);}
if(num_tops>1){printf("The proportion of variance explained by the %d top predictor is %.4f\n\n", num_tops, topher);}
}
else{topher=0;}
*/

