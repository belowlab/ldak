/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Data operations

///////////////////////////

void open_datagz(gzFile *inputgz, char *datafile, int num_samples, int genskip, int genheaders, int genprobs)
{
int j, count, size, size2;
char readchar, readchar2;
char *gzbuffer;


size=30000000+num_samples*(genprobs*20+(genprobs==0)*4);
gzbuffer=malloc(sizeof(char)*size);

//check size - open and read in until newline (or gzbuffer read)
if((*inputgz=gzopen(datafile,"r"))==NULL)
{printf("Error opening %s\n\n",datafile);exit(1);}

if(gzgets(*inputgz,gzbuffer,size)==NULL)
{printf("Error reading %s; file appears to be empty\n\n", datafile);exit(1);}
size2=strlen(gzbuffer);

if(gzbuffer[0]==9||gzbuffer[0]==10||gzbuffer[0]==32)
{printf("Error, %s starts with a space or empty row\n", datafile);exit(1);}
if(size2==size-1)
{printf("Error reading %s; Row 1 is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", datafile, size2, size);exit(1);}

//get number of elements in row
count=1;readchar2=gzbuffer[0];
for(j=1;j<size2-1;j++)
{
readchar=gzbuffer[j];
if((readchar2==9||readchar2==32)&&readchar!=9&&readchar!=32){count++;}
readchar2=readchar;
}
gzclose(*inputgz);

if(genprobs>0&&count!=genheaders+num_samples*genprobs)
{printf("Error reading %s; should contain %d columns, not %d\n\n", datafile, genheaders+num_samples*genprobs, count);exit(1);}

if(genprobs==0&&count!=genheaders+num_samples*2)
{printf("Error reading %s; should contain %d columns, not %d\n\n", datafile, genheaders+num_samples*2, count);exit(1);}

//open and skip
if((*inputgz=gzopen(datafile,"r"))==NULL)
{printf("Error opening %s\n\n",datafile);exit(1);}
count=0;
while(count<genskip)
{
if(gzgets(*inputgz,gzbuffer,size)==NULL)
{printf("Error reading Header Row %d of %s\n\n", count+1, datafile);exit(1);}
size2=strlen(gzbuffer);

if(size2>=size-1)
{printf("Error reading %s; Row %d is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", datafile, count+1, size2, size);exit(1);}
printf("Skipping Row %d of %s, I hope this is a header\n(Row begins: %c", count+1, datafile, gzbuffer[0]);
for(j=1;j<size2-1;j++)
{
if(j<50){printf("%c", gzbuffer[j]);}
}
printf(")\n");
count++;
}
if(count>0){printf("\n");}

free(gzbuffer);
}	//end of open_datagz_fly

////////

int read_data_fly(char *datafile, int dtype, double *data, float **probs, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds_use, gzFile inputgz, size_t current, int num_samples, int num_preds, int genskip, int genheaders, int genprobs, double missingvalue, double threshold, double minprob, int nonsnp)
{
if(start>=end){return(current);}

if(dtype==1)
{read_bed_fly(datafile, data, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue);}

//if(dtype==2)
//{current=read_sp_fly(datafile, data, num_samples_use, keepsamps, start, end, keeppreds_use, inputgz, current, num_samples, num_preds, missingvalue, threshold);}

if(dtype==3)
{read_sped_fly(datafile, data, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue, threshold, nonsnp);}

if(dtype==4)
{read_speed_fly(datafile, data, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue, threshold, nonsnp);}

if(dtype==5)
{current=read_gen_fly(datafile, data, probs, num_samples_use, keepsamps, start, end, keeppreds_use, inputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, threshold, minprob, nonsnp);}

return(current);
}

///////////////////////////

void stand_data(double *data, double *centres, double *mults, double *sqdevs, int ns, int length, double missingvalue, double power, int gotcent, int hwestand, double *weights, int type, char **preds)
//type=0 - just calc centre and mult
//type=1 - calc centre and mult, and standardize
//type=2 - standardize using provided centre and mult
{
int i, j, indcount, wcount;
double sum, sumsq, mean, var, value;


wcount=0;
#pragma omp parallel for private(j,sum,sumsq,indcount,i,mean,var,value) schedule (static)
for(j=0;j<length;j++)
{
if(type==0||type==1)	//calculate centre, mult and sqdev
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<ns;i++)
{
if(data[(size_t)j*ns+i]!=missingvalue){sum+=data[(size_t)j*ns+i];sumsq+=pow(data[(size_t)j*ns+i],2);indcount++;}
}

if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=0;var=0;}

if(gotcent==0){centres[j]=mean;}
else{mean=centres[j];}
if(var>0)
{
if(hwestand==1){mults[j]=pow(mean*(1-mean/2),power/2);}
else{mults[j]=pow(var*indcount/ns,power/2);}
}
else	//trivial
{
mults[j]=-9999;
if(wcount<5){printf("Warning, Predictor %s is trivial (takes at most one non-missing value) and will be ignored\n", preds[j]);}
wcount++;
}
sqdevs[j]=var*indcount/ns;
}	//end of calculating centre, mult and sqdevs

if(type==1||type==2)	//standardize (with weights if provided)
{
if(weights!=NULL)	//using weights
{
if(mults[j]!=-9999&&weights[j]>0)	//not trivial
{
value=mults[j]*pow(weights[j],.5);
for(i=0;i<ns;i++)
{
if(data[(size_t)j*ns+i]!=missingvalue){data[(size_t)j*ns+i]=(data[(size_t)j*ns+i]-centres[j])*value;}
else{data[(size_t)j*ns+i]=0;}
}
}
else	//trivial
{
for(i=0;i<ns;i++){data[(size_t)j*ns+i]=0;}
}
}
else	//not using weights
{
if(mults[j]!=-9999)	//not trivial
{
for(i=0;i<ns;i++)
{
if(data[(size_t)j*ns+i]!=missingvalue){data[(size_t)j*ns+i]=(data[(size_t)j*ns+i]-centres[j])*mults[j];}
else{data[(size_t)j*ns+i]=0;}
}
}
else	//trivial
{
for(i=0;i<ns;i++){data[(size_t)j*ns+i]=0;}
}
}
}	//end of type=1/2
}	//end of j loop

if(wcount>5){printf("In total, %d predictors are trivial\n", wcount);}
if(wcount>0){printf("\n");}
}	//end of stand_data

///////////////////////////

int prune_regions(int num_regs, int **regindex, double *rmults, double rprune, double *rdata, int ns)
{
int i, j, j2, r, count, count2, count3;
int *retain;
double sum, sumsq, mean, var, value, alpha, beta;
double *rdata2, *cors;


count2=0;for(r=0;r<num_regs;r++){count2+=regindex[r][0];}

//first deal with trivial
for(r=0;r<num_regs;r++)
{
count=0;
for(j=0;j<regindex[r][0];j++)
{
j2=regindex[r][1+j];
if(rmults[j2]!=-9999){regindex[r][1+count]=j2;count++;}
}
if(count==0)
{printf("Error, all %d predictors in Region %d are trivial\n", regindex[r][0], r+1);exit(1);}
regindex[r][0]=count;
}

count3=0;for(r=0;r<num_regs;r++){count3+=regindex[r][0];}
if(count3<count2)
{printf("After removing trivial predictors, the number of region predictors is reduced to %d\n", count3);}

if(rprune>1){return(count3);}

////////

printf("Pruning region predictors based on a correlation squared threshold of %f\nTo avoid this, set \"--region-prune\" to a value greater than 1\n", rprune);

for(r=0;r<num_regs;r++)
{
retain=malloc(sizeof(int)*regindex[r][0]);
for(j=0;j<regindex[r][0];j++){retain[j]=1;}

rdata2=malloc(sizeof(double)*ns*regindex[r][0]);
cors=malloc(sizeof(double)*regindex[r][0]*regindex[r][0]);

for(j=0;j<regindex[r][0];j++)
{
j2=regindex[r][1+j];
sum=0;sumsq=0;
for(i=0;i<ns;i++)
{sum+=rdata[(size_t)j2*ns+i];sumsq+=pow(rdata[(size_t)j2*ns+i],2);}
mean=sum/ns;
var=sumsq/ns-pow(mean,2);
value=pow(var,-.5);
for(i=0;i<ns;i++){rdata2[(size_t)j*ns+i]=(rdata[(size_t)j2*ns+i]-mean)*value;}
}

alpha=1.0/ns;beta=0.0;
dgemm_("T", "N", regindex[r], regindex[r], &ns, &alpha, rdata2, &ns, rdata2, &ns, &beta, cors, regindex[r]);

for(j=0;j<regindex[r][0];j++)
{
if(retain[j]==1)
{
for(j2=j+1;j2<regindex[r][0];j2++)
{
if(retain[j2]==1)
{
if(pow(cors[(size_t)j2*regindex[r][0]+j],2)>rprune){retain[j2]=0;}
}}
}}

count=0;
for(j=0;j<regindex[r][0];j++)
{
if(retain[j]==1){regindex[r][1+count]=regindex[r][1+j];count++;}
}
regindex[r][0]=count;

free(retain);
free(rdata2);free(cors);
}	//end of r loop

count2=0;for(r=0;r<num_regs;r++){count2+=regindex[r][0];}
if(count2<count3){printf("Total number of region predictors reduced to %d\n\n", count2);}
else{printf("No predictors removed\n\n");}

return(count2);
}	//end of prune_regions

///////////////////////////

int fill_X(double *X, int *Xstarts, int *Xends, double *Xsums, double *Xnss, double *Xrhos, int ns, int num_regs, int **regindex, double *rdata, double *rnss, double *rrhos)
{
int i, j, j2, r, count;
double sumsq;


count=0;
for(r=0;r<num_regs;r++)
{
Xstarts[r]=count;
sumsq=0;
for(j=0;j<regindex[r][0];j++)
{
j2=regindex[r][1+j];
for(i=0;i<ns;i++)
{
X[i+count*ns]=rdata[(size_t)j2*ns+i];
sumsq+=pow(X[i+count*ns],2);
}
if(Xnss!=NULL){Xnss[count]=rnss[j2];Xrhos[count]=rrhos[j2];}
count++;
}
Xends[r]=count;
Xsums[r]=sumsq/ns;
}

return(count);
}

///////////////////////////

int extraction(int *usedpreds, int length, char **preds, int *chr, int *predorder, char *bpredfile, char *cpredfile, int onechr, char *onesnp, char *filename)
{
int j, count, count2;
int *indexer;
char **wantpreds;


count=0;
if(strcmp(bpredfile,"blank")!=0)
{
count=countrows(bpredfile);
printf("Reading list of %d predictors to extract from %s\n", count, bpredfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(bpredfile, wantpreds, count, NULL, 1, 0);
}
if(strcmp(onesnp,"blank")!=0)
{
count=1;
wantpreds=malloc(sizeof(char*));
wantpreds[0]=malloc(sizeof(char)*(strlen(onesnp)+1));
strcpy(wantpreds[0],onesnp);
}

if(count>0)
{
indexer=malloc(sizeof(int)*count);
count2=find_strings(preds, length, wantpreds, count, indexer, NULL, NULL, NULL, predorder, NULL, 3);
if(count2==0)
{
if(strcmp(bpredfile,"blank")!=0){printf("Error, none of these are in %s\n\n", filename);exit(1);}
else{printf("Error, Predictor %s is not in %s\n\n", onesnp, filename);exit(1);}
}

if(count2<count){printf("Warning, only %d of these are in %s\n", count2, filename);}
for(j=0;j<count2;j++){usedpreds[indexer[j]]++;}
for(j=0;j<length;j++){usedpreds[j]=(usedpreds[j]==2);}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}

if(strcmp(cpredfile,"blank")!=0)	//exclude predictors
{
count=countrows(cpredfile);
printf("Reading list of %d predictors to exclude from %s (takes priority over any other predictor filtering)\n", count, cpredfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(cpredfile, wantpreds, count, NULL, 1, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(preds, length, wantpreds, count, indexer, NULL, NULL, NULL, predorder, NULL, 3);
if(count2==0){printf("Warning, none of these are in %s\n", filename);}
if(count2>0&&count2<count){printf("Warning, only %d of these are in %s\n", count2, filename);}
for(j=0;j<count2;j++){usedpreds[indexer[j]]=0;}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}

if(onechr>0)
{
count=0;
for(j=0;j<length;j++)
{
if(chr[j]==onechr){count++;}
else{usedpreds[j]=0;}
}
if(count==0)
{printf("Error, there are no Chromosome %d predictors in %s\n\n", onechr, filename);exit(1);}
printf("There are %d Chromosome %d predictors\n", count, onechr);
}
if(onechr==-1)
{
count=0;
for(j=0;j<length;j++)
{
if(chr[j]>0&&chr[j]<23){count++;}
else{usedpreds[j]=0;}
}
if(count==0){printf("Error, there are no autosomal predictors in %s\n\n", filename);exit(1);}
printf("There are %d autosomal predictors\n", count);
}

count=0;for(j=0;j<length;j++){count+=usedpreds[j];}
if(count==0){printf("Error, no predictors remain\n\n");exit(1);}

return(count);
}	//end of extraction

///////////////////////////

double get_sum_stats(double *stags, double *snss, double *schis, double *srhos, int length, int type, double cutoff)
{
int j;
double sum, sum2, mean, value, value2;
struct sorting_double *dptrs;


//get statistics - get sum of weights and mean, then max variance
sum=0;for(j=0;j<length;j++){sum+=pow(stags[j],-1);}
sum2=0;for(j=0;j<length;j++){sum2+=schis[j]/stags[j];}
mean=sum2/sum;
value=0;
for(j=0;j<length;j++)
{
if(pow(srhos[j],2)>value){value=pow(srhos[j],2);}
}

//now sort and get median
dptrs=malloc(sizeof(struct sorting_double)*length);
for(j=0;j<length;j++){dptrs[j].value=schis[j];dptrs[j].index=j;}
qsort(dptrs, length, sizeof(struct sorting_double), compare_sorting_double);

sum2=0;
for(j=0;j<length;j++)
{
value2=dptrs[j].value;
sum2+=pow(stags[dptrs[j].index],-1);
if(sum2>sum/2){break;}
}

//finally get average sample size
sum2=0;for(j=0;j<length;j++){sum2+=snss[j]/stags[j];}

if(type==0){printf("Max ");}
else{printf("Trait %d: max ", type);}
printf("variance explained %.4f, max test statistic %.2f, weighted mean %.3f, weighted median %.3f, weighted GIF %.3f, weighted sample size %.2f\n\n", value, dptrs[length-1].value, mean, value2, value2/0.4549364, sum2/sum);

if(cutoff==-9999&&value>.01){printf("Warning, estimates of heritability can be biased by the presence of strong-effect loci, so consider using \"--cutoff\" to exclude them (for example, \"--cutoff 0.01\" will exclude predictors which individually explain at least 1%% of phenotypic variation)\n\n");}

free(dptrs);

return(value2/0.4549364);
}

///////////////////////////

void load_data_chunk(unsigned char **data_char, int dtype, double *data, int ns, int start, int end, int *order, double **bytelookup, float *speedstarts, float *speedscales, double *centres)
{
int i, j, j2, readint;

for(j=start;j<end;j++)
{
j2=j;
if(order!=NULL){j2=order[j];}
if(dtype==1)	//can use lookups
{
for(i=0;i<ns;i+=4){memcpy(data+(j-start)*ns+i,bytelookup[j2]+4*(int)data_char[j2][i/4],sizeof(double)*4);}
}
else	//do manually
{
for(i=0;i<ns;i++)
{
readint=(int)data_char[j2][i];
if(readint!=255){data[(size_t)(j-start)*ns+i]=(speedstarts[j2]+speedscales[j2]*readint-centres[j2]);}
else{data[(size_t)(j-start)*ns+i]=0;}
}
}
}

}

///////////////////////////

