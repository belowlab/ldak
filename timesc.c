/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//do work for calculating scores

///////////////////////////

if(bitsize>data_length){bitsize=data_length;}

//allocate variables
if(strcmp(proffile,"blank")==0)
{
data_warn2(bitsize,(1+savecounts)*num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
if(savecounts==1){data2=malloc(sizeof(double)*num_samples_use*bitsize);}

nummiss=malloc(sizeof(int)*data_length);
blupmids=malloc(sizeof(double)*bitsize*num_scores);
blupprojs=malloc(sizeof(double)*num_samples_use*num_scores);
}

guesses=malloc(sizeof(double*)*(num_scores+1));	//last one is covs
for(k=0;k<num_scores+1;k++){guesses[k]=malloc(sizeof(double)*num_samples_use);}

if(savecounts==1)
{
nums=malloc(sizeof(int*)*num_scores);
for(k=0;k<num_scores;k++){nums[k]=malloc(sizeof(int)*num_samples_use);}
}

if(strcmp(respfile,"blank")!=0||strcmp(sumsfile,"blank")!=0){predcors=malloc(sizeof(double)*num_scores);}

if(strcmp(finalfile,"blank")!=0)
{
predsds=malloc(sizeof(double)*num_scores);
sP=malloc(sizeof(double)*num_samples_use*num_scores);
sPTP=malloc(sizeof(double)*num_scores*num_scores);
sPTP2=malloc(sizeof(double)*num_scores*num_scores);
sPTY=malloc(sizeof(double)*num_scores);
}

if(strcmp(proffile,"blank")!=0)	//have previously computed profiles, so read from file
{
count=countcols(proffile);
if(count!=2*num_scores+4)
{printf("Error, %s should have %d columns (not %d); make sure it was produced using the score file \"%s\"\n\n", proffile, 2*num_scores+4, count, scorefile);exit(1);}
count=countrows(proffile);
if(count!=1+num_samples_use)
{printf("Error, %s should have %d rows (not %d); make sure you are using the same data and sample filterings now as when you constructed the file\n\n", proffile, 1+num_samples_use, count);exit(1);}

printf("Reading observed and predicted phenotypes for %d samples from %s\n\n", count2, proffile);
if((input=fopen(proffile,"r"))==NULL)
{printf("Error opening %s\n\n",proffile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

count=0;
for(i=0;i<num_samples_use;i++)
{
if(fscanf(input, "%s %s %s %lf ", readstring, readstring2, readstring3, guesses[num_scores]+i)!=4)
{printf("Error reading first four values of Row %d of %s\n\n", j+1, proffile);exit(1);}
if(strcmp(readstring,ids1[i])!=0||strcmp(readstring2,ids2[i])!=0)
{printf("Error, Row %d of %s should correspond to Samples %s and %s (not %s and %s); make sure you are using the same data and sample filterings now as when you constructed the file\n\n", i+2, proffile, ids1[i], ids2[i], readstring, readstring2);exit(1);}

for(k=0;k<num_scores;k++)
{
if(fscanf(input, "%lf %s ", guesses[k]+i, readstring)!=2)
{printf("Error reading Profile %d from Row %d of %s\n\n", k+1, i+1, proffile);exit(1);}
}
}
fclose(input);
}
else	//must compute from scratch
{
//set guesses and possibly nums to zero
for(k=0;k<num_scores+1;k++)
{
for(i=0;i<num_samples_use;i++){guesses[k][i]=0;}
}
if(savecounts==1)
{
for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){nums[k][i]=0;}
}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating scores for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating scores for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

//read data
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);

//standardize, keeping track of number of missing values
wcount=0;
for(j=0;j<bitlength;j++)
{
//fill nummiss, centres, mults and sqdevs
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{sum+=data[(size_t)j*num_samples_use+i];sumsq+=pow(data[(size_t)j*num_samples_use+i],2);indcount++;}
}
if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=0;var=0;}

nummiss[bitstart+j]=num_samples_use-indcount;
if(blupcentres[0][bitstart+j]==-9999){centres[bitstart+j]=mean;}
else{centres[bitstart+j]=blupcentres[0][bitstart+j];mean=blupcentres[0][bitstart+j];}
if(var>0)
{
if(hwestand==1){mults[bitstart+j]=pow(mean*(1-mean/2),power/2);}
else{mults[bitstart+j]=pow(var*indcount/num_samples_use,power/2);}
}
else	//trivial
{
if(power==0&&blupcentres[0][bitstart+j]!=-9999&&strcmp(sumsfile,"blank")==0)	//can still use
{mults[bitstart+j]=1;}
else
{
mults[bitstart+j]=-9999;
if(wcount<5){printf("Warning, Predictor %s is trivial (takes at most one non-missing value) and will be ignored\n", preds[bitstart+j]);}
wcount++;
}
}
sqdevs[bitstart+j]=var*indcount/num_samples_use;

if(savecounts==1)	//must copy data into data2
{
for(i=0;i<num_samples_use;i++){data2[(size_t)j*num_samples_use+i]=data[(size_t)j*num_samples_use+i];}
}

//now standardize
if(mults[bitstart+j]!=-9999)	//not trivial
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{data[(size_t)j*num_samples_use+i]=(data[(size_t)j*num_samples_use+i]-centres[bitstart+j])*mults[bitstart+j];}
else
{data[(size_t)j*num_samples_use+i]=0;}
}
}
else	//trivial
{
for(i=0;i<num_samples_use;i++){data[(size_t)j*num_samples_use+i]=0;}
}
}	//end of j loop

//load scaled effect sizes into blupmids
for(k=0;k<num_scores;k++)
{
for(j=0;j<bitlength;j++){blupmids[j+k*bitlength]=blupfactors[k][bitstart+j];}
}

//get contribution to projections
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_samples_use, &num_scores, &bitlength, &alpha, data, &num_samples_use, blupmids, &bitlength, &beta, blupprojs, &num_samples_use);

//add these to guesses
for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){guesses[k][i]+=blupprojs[i+k*num_samples_use];}
}

if(savecounts==1)	//deal with nums
{
//first increase nums assuming all relevant values are present for each score
for(k=0;k<num_scores;k++)
{
count=0;for(j=0;j<bitlength;j++){count+=(blupfactors[k][bitstart+j]!=0&&mults[bitstart+j]!=-9999);}
for(i=0;i<num_samples_use;i++){nums[k][i]+=count;}
}

//now loop through values subtracting from nums for each missing
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999&&nummiss[bitstart+j]>0)	//there are some missing values
{
for(i=0;i<num_samples_use;i++)
{
if(data2[(size_t)j*num_samples_use+i]==missingvalue)
{
for(k=0;k<num_scores;k++){nums[k][i]-=(blupfactors[k][bitstart+j]!=0);}
}}
}}
}
}	//end of bit loop
printf("\n");

if(strcmp(cofile,"blank")!=0)	//get contribution of covariates
{
alpha=1.0;beta=0.0;
dgemv_("N", &num_samples_use, &num_covars, &alpha, covar, &num_samples_use, thetas, &one, &beta, guesses[num_scores], &one);
}
//else will remain zero

//save
if(savecounts==0){write_scores(guesses, NULL, num_samples_use, num_scores, ids1, ids2, outfile, resp, missingvalue);}
else{write_scores(guesses, nums, num_samples_use, num_scores, ids1, ids2, outfile, resp, missingvalue);}
}	//end of computing profiles

////////

if(strcmp(respfile,"blank")!=0||strcmp(sumsfile,"blank")!=0)	//get correlations
{
for(r=0;r<num_scores;r++)
{
if(strcmp(respfile,"blank")!=0)
{
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(resp[i]!=missingvalue)
{
sum+=guesses[r][i];sum2+=resp[i]-guesses[num_scores][i];
sumsq+=pow(guesses[r][i],2);sumsq2+=pow(resp[i]-guesses[num_scores][i],2);
sumsq3+=guesses[r][i]*(resp[i]-guesses[num_scores][i]);
indcount++;
}}

if(sumsq>sum*sum/indcount)	//score might be trivial (but have checked phenotype is not, so indcount>1)
{
mean=sum/indcount;mean2=sum2/indcount;
predcors[r]=(sumsq3-indcount*mean*mean2)/pow(sumsq-indcount*mean*mean,.5)/pow(sumsq2-indcount*mean2*mean2,.5);
}
else{predcors[r]=-9999;}
}
else	//so have summaries
{
//get cov(score,Y)/SD(Y) = weighted sum of cov(Xj,Y)/SD(Y) - cov(Xj,Y)/SD(Y) = rhos[j]*SD(Xj)
value=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){value+=blupfactors[r][j]*rhos[j]*pow(sqdevs[j],.5);}
}
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=guesses[r][i];sumsq+=pow(guesses[r][i],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

if(var>0){predcors[r]=value*pow(var,-.5);}
else{predcors[r]=-9999;}
}
}	//end of r loop

sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
for(r=0;r<num_scores;r++)
{
if(predcors[r]!=-9999){fprintf(output2,"Score_%d %.6f\n", r+1, predcors[r]);}
else{fprintf(output2,"Score_%d NA\n", r+1);}
}
fclose(output2);

if(num_scores==1){printf("Correlation between score and phenotype is %.4f, saved in %s\n\n", predcors[0], filename2);}
else
{
min=predcors[0];max=predcors[0];
for(r=1;r<num_scores;r++)
{
if(predcors[r]<min){min=predcors[r];}
if(predcors[r]>max){max=predcors[r];}
}
printf("Correlations between %d scores and phenotype range from %.4f to %.4f, saved in %s\n\n", num_scores, min, max, filename2);
}
}	//end of getting correlations for scores

////////

if(strcmp(finalfile,"blank")!=0)	//get weights (must also have phenotype or summaries)
{
//first get best (checking at least one non-trivial score)
best=-9999;value=-9999;
for(r=0;r<num_scores;r++)
{
if(predcors[r]>value){best=r;value=predcors[r];}
}
if(best==-9999){printf("Error, all %d scores are trivial\n\n", num_scores);exit(1);}

//create new effects

count=countrows(finalfile)-1;
printf("Reading effects for %d predictors from %s\n\n", count, finalfile);

if((input=fopen(finalfile,"r"))==NULL)
{printf("Error opening %s\n\n",finalfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

sprintf(filename3,"%s.effects.best",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Centre Effect_Best\n");

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %s %s %s ", readstring, readstring2, readstring3, readstring4)!=4)
{printf("Error reading first four values of Row %d of %s\n\n", j+1, finalfile);exit(1);}

for(r=0;r<num_scores;r++)
{
if(fscanf(input, "%lf ", &readdouble)!=1)
{printf("Error reading Effect %d from Row %d of %s\n\n", r+1, j+1, finalfile);exit(1);}
if(r==best){value=readdouble;}
}

fprintf(output3, "%s %s %s %s %.6f\n", readstring, readstring2, readstring3, readstring4, value);
}
fclose(input);
fclose(output3);

printf("The best-fitting model is saved in %s\n\n", filename3); 
}

////////

//free allocations from setdl.c
for(k=0;k<num_scores;k++){free(blupcentres[k]);free(blupfactors[k]);}free(blupcentres);free(blupfactors);

//frees from above
if(strcmp(proffile,"blank")==0)
{
free(data);
if(savecounts==1){free(data2);}
free(nummiss);free(blupmids);free(blupprojs);
}
for(k=0;k<num_scores+1;k++){free(guesses[k]);}free(guesses);
if(savecounts==1)
{
for(k=0;k<num_scores;k++){free(nums[k]);}free(nums);
}

if(strcmp(respfile,"blank")!=0||strcmp(sumsfile,"blank")!=0){free(predcors);}
if(strcmp(finalfile,"blank")!=0){free(predsds);free(sP);free(sPTP);free(sPTP2);free(sPTY);}
if(strcmp(proffile,"blank")==0&&binary==0){gzclose(datainputgz);}

///////////////////////////

