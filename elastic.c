/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//elastic net (and lasso, ridge regression) - when changing this, consider changing also bolt.c and bayesr.c

//this was written, assuming priors divided by e = varphen(1-her) - must change

///////////////////////////

//set num_train, num_test and keepboths (index training then test samples, relative to ids, then allids)
keepboth=malloc(sizeof(int)*num_samples_use);
keepboth2=malloc(sizeof(int)*num_samples_use);

if(skipcv==0)	//sort cv samples
{
if(cvprop!=-9999)
{
num_test=cvprop*num_samples_use;
num_train=num_samples_use-num_test;
for(i=0;i<num_samples_use;i++){keepboth[i]=i;}
permute_int(keepboth,num_samples_use);
}

if(strcmp(bvsfile,"blank")!=0)
{
count=countrows(bvsfile);
printf("Reading list of %d cross-validation samples from %s\n", count, bvsfile);
wantids=malloc(sizeof(char*)*count);
read_ids(bvsfile, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
num_test=find_strings(ids3, num_samples_use, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(num_test==0){printf("Error, can not find any of these samples in the data\n\n");exit(1);}
if(num_test<count){printf("Warning, only %d of these are in the data\n", num_test);}
num_train=num_samples_use-num_test;

usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(i=0;i<num_test;i++){usedids[indexer[i]]=1;}
count2=0;
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==0){keepboth[count2]=i;count2++;}
}
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==1){keepboth[count2]=i;count2++;}
}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(usedids);free(indexer);
}

printf("Will be using %d samples to train and %d to test\n\n", num_train, num_test);
if(num_train<3){printf("Error, unable to continue with fewer than three training samples\n\n");exit(1);}
}
else
{
num_train=num_samples_use;
num_test=0;
for(i=0;i<num_samples_use;i++){keepboth[i]=i;}
}

for(i=0;i<num_samples_use;i++){keepboth2[i]=keepsamps[keepboth[i]];}

//work out sequence of heritabilities
num_hers=1+(hermax-1e-4)/herstep;
tryhers=malloc(sizeof(double)*num_hers);
for(j=0;j<num_hers;j++){tryhers[j]=hermax+(j+1-num_hers)*herstep;}

//allocate other variables
if(dtype==1){count=(num_samples_use-1)/4+1;}
else{count=num_samples_use;}
value=(double)data_length/1024/1024/1024*count;
if(value>1){printf("Warning, to read the data requires %.1f Gb (this can not be reduced)\n\n", value);}

data_char=malloc(sizeof(unsigned char*)*data_length);
for(j=0;j<data_length;j++){data_char[j]=malloc(sizeof(unsigned char)*count);}

if(dtype==4){speedstarts=malloc(sizeof(float)*data_length);speedscales=malloc(sizeof(float)*data_length);}

Y=malloc(sizeof(double)*num_samples_use);
Y2=malloc(sizeof(double)*num_samples_use);
YTdata=malloc(sizeof(double)*data_length);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
datasqs=malloc(sizeof(double)*data_length);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

exps=malloc(sizeof(double)*data_length);
order=malloc(sizeof(int)*data_length);

lambdas=malloc(sizeof(double)*data_length);
lambdas2=malloc(sizeof(double)*data_length);

evers=malloc(sizeof(int)*data_length);
strongs=malloc(sizeof(int)*data_length);
effs=malloc(sizeof(double)*data_length);
residuals=malloc(sizeof(double)*num_samples_use);

if(skipcv==0){effs2=malloc(sizeof(double)*data_length);}

//add variables for speeding up code here

////////

//we begin using all samples (training then test)

//read in the data - will be using either bed format or short SPEED format (dtype 1 or 4)
if(dtype==1)
{(void)read_bed_full(datafile, data_char, num_samples_use, keepboth2, data_length, keeppreds_use, num_samples, num_preds);}
else
{(void)read_speed_full(datafile, speedstarts, speedscales, data_char, NULL, num_samples_use, keepboth2, data_length, keeppreds_use, num_samples, num_preds, nonsnp);}

//fill Y and Z, get thetas and residuals
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[keepboth[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[keepboth[i]+j*num_samples_use];}
}
reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Y2, 0, NULL, NULL);

//get variance of residuals (Y2 will have mean zero)
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(Y2[i],2);}
varphen=sumsq/num_samples_use;

////////

//get means and variances - mults indicates which predictors non-trivial (will subsequently use exps)
printf("Computing predictor means and variances\n\n");
wcount=0;
for(j=0;j<data_length;j++)
{
if(dtype==1)	//bed format - actual values stored (3 means missing)
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
readint=(int)((data_char[j][i/4] >> (2*(i%4))) & 3);
if(readint!=3){sum+=readint;sumsq+=pow(readint,2);indcount++;}
}
}
else	//short speed format - values are speedstart + value * speedscale
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
readint=(int)data_char[j][i];
if(readint!=255){sum+=readint;sumsq+=pow(readint,2);indcount++;}
}
sumsq=indcount*pow(speedstarts[j],2)+2*speedstarts[j]*speedscales[j]*sum+pow(speedscales[j],2)*sumsq;
sum=indcount*speedstarts[j]+speedscales[j]*sum;
}

if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=-9999;var=0;}

centres[j]=mean;
if(var>0){mults[j]=1;}
else	//trivial
{
mults[j]=-9999;
if(wcount<5){printf("Warning, Predictor %s is trivial (takes at most one non-missing value) and will be ignored\n", preds[j]);}
wcount++;
}
sqdevs[j]=var*indcount/num_samples_use;
}	//end of j loop
if(wcount==data_length){printf("Error, all predictors are trivial\n\n");exit(1);}
if(wcount>5){printf("In total, %d predictors are trivial\n", wcount);}
if(wcount>0){printf("\n");}

//make lookup table

////////

//get expected squared effect size of each predictor, divided by varphen * hers

//start by getting expected heritabilities (set to zero for trivial predictors)
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

//scale exps so they sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}

//now divide by (expected) predictor variance
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(hwestand==1){exps[j]*=pow(centres[j]*(1-centres[j]/2),-1);}
else{exps[j]*=pow(sqdevs[j],-1);}
}
}

//finally, get reverse order (will not change when we later multiply by varphens and hers)
dptrs=malloc(sizeof(struct sorting_double)*data_length);
for(j=0;j<data_length;j++){dptrs[j].value=exps[j];dptrs[j].index=j;}
qsort(dptrs,data_length,sizeof(struct sorting_double),compare_sorting_double_rev);
for(j=0;j<data_length;j++){order[j]=dptrs[j].index;}
free(dptrs);

////////

//first solve using only training samples (although use all samples when reading data and computing residuals)

//get XTX for each predictor
for(j=0;j<data_length;j++)
{
datasqs[j]=0;
if(exps[j]>0)
{
if(dtype==1)	//bed format - actual values stored (3 means missing)
{
for(i=0;i<num_train;i++)
{
readint=(int)((data_char[j][i/4] >> (2*(i%4))) & 3);
if(readint!=3){datasqs[j]+=pow(readint-centres[j],2);}
}
}
else	//short speed format - values are speedstart + value * speedscale
{
for(i=0;i<num_train;i++)
{
readint=(int)data_char[j][i];
if(readint!=255){datasqs[j]+=pow(speedstarts[j]+speedscales[j]*readint-centres[j],2);}
}
}
}}

if(enalpha!=0)	//for elastic net / lasso, set lambdas2 to a value where no predictors enter model
{
value=tryhers[0]/2;
set_lambdas(lambdas2, exps, data_length, order, value, varphen, enalpha);

//get YTdata - could optimize and use dgemv
for(j=0;j<data_length;j++)
{
if(exps[j]>0){YTdata[j]=get_YTdata(j, Y2, dtype, data_char, num_train, speedstarts, speedscales, centres);}
else{YTdata[j]=0;}
}

while(1)
{
flag=0;
for(j=0;j<data_length;j++)
{
if(fabs(YTdata[j])>enalpha*lambdas2[j]){flag=1;break;}
}
if(flag==0){break;}
value=value/2;
set_lambdas(lambdas2, exps, data_length, order, value, varphen, enalpha);
}
}

//start with effect sizes set at zero, residuals those after regressing out covariates, and ever set empty
for(j=0;j<data_length;j++){effs[j]=0;}
for(i=0;i<num_samples_use;i++){residuals[i]=Y2[i];}
total=0;for(j=0;j<data_length;j++){evers[j]=0;}	
flag=1;	//indicates whether YTdata correct (only used with elastic net / ridge)
minmse=-9999;	//minimum mse from cross-valdiation

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Loop\tHeritability\tNum_Non_Zero\tNum_Strong_Set\tNum_Iters_One\tNum_Iters_Two\tConverged\tDifference\tMSE\n");
fclose(output);

////////

//ready to do loops of iterations

for(p=0;p<num_hers;p++)	//solve for tryhers[p]
{
printf("Loop %d of %d (heritability %.6f)\n", p+1, num_hers, tryhers[p]);

//set lambdas (for elastic net / lasso, lambdas2 already set)
set_lambdas(lambdas, exps, data_length, order, tryhers[p], varphen, enalpha);

//get pen and (partial) likelihood - at bottom is code to get likelihood constant
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[p]);

if(enalpha!=0)	//lasso / elastic net - using posterior modes
{
if(flag==0)	//residuals have changed, so update YTdata
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0){YTdata[j]=get_YTdata(j, residuals, dtype, data_char, num_train, speedstarts, speedscales, centres);}
else{YTdata[j]=0;}
}
flag=1;
}

//fill strong vector (indicates which predictors we will loop over) - |YTdata| < alpha (2lam -lam2)
total2=0;
for(j=0;j<data_length;j++)	//trivial predictors will have evers and YTdata zero
{
if(evers[j]==1||fabs(YTdata[j])>=enalpha*(2*lambdas[j]-lambdas2[j])){strongs[j]=1;total2++;}
else{strongs[j]=0;}
}

if(total>0)	//update effects for predictors ever in model
{
printf("Updating effects for %d predictors that have been in at least one model\n", total);
count=0;
while(1)
{
count++;
for(j=0;j<data_length;j++)
{
if(evers[j]==1)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_train, speedstarts, speedscales, centres);
value=0;
if(sum>enalpha*lambdas[j]){value=(sum-enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}
if(-sum>enalpha*lambdas[j]){value=(sum+enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}

if(value!=effs[j])	//effect size changed, so update residuals (for all samples) and flag
{
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
flag=0;
effs[j]=value;
}
}}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[p]);

if(total==1){break;}
if(count>1)
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;

if(count>=20){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count, diff, tol);break;}
}
}
else{count=0;printf("No predictors in the model\n");}

if(total2>0)	//update effects for predictors in strong set
{
printf("Updating effects for %d predictors in strong set\n", total2);
count2=0;
cflag=0;
while(1)	
{
count2++;
for(j=0;j<data_length;j++)
{
if(strongs[j]==1)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_train, speedstarts, speedscales, centres);
value=0;
if(sum>enalpha*lambdas[j]){value=(sum-enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}
if(-sum>enalpha*lambdas[j]){value=(sum+enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}

if(value!=effs[j])	//effect size changed, so update residuals (for all samples) and flag
{
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
flag=0;
effs[j]=value;
}
}}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[p]);

if(total2==1){diff=0;break;}
if(count2>1)
{
diff=like-likeold;
if(fabs(diff)<tol){cflag=1;break;}
}
likeold=like;

if(count2>=maxiter){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count2, diff, tol);break;}
}
}
else{count2=0;diff=0;cflag=1;printf("No predictors in the strong set\n");}

//could perform kkt test for all predictors, but experience indicates unnecessary

//update ever active set
total=0;
for(j=0;j<data_length;j++)
{
if(effs[j]!=0){evers[j]=1;}
total+=evers[j];
}

//put lambdas into lambdas2 ready for next iterations
for(j=0;j<data_length;j++){lambdas2[j]=lambdas[j];}
}	//end of enalpha!=0
else	//ridge - using posterior means
{
//update effects for all non-trivial predictors
count=0;
cflag=0;
while(1)	
{
count++;
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_train, speedstarts, speedscales, centres);
value=sum/(datasqs[j]+lambdas[j]);
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
effs[j]=value;
}
}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[p]);

if(count>1)
{
diff=like-likeold;
if(fabs(diff)<tol){cflag=1;break;}
}
likeold=like;

if(count>=maxiter){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count, diff, tol);break;}
}
}

//get size of model
found=0;for(j=0;j<data_length;j++){found+=(effs[j]!=0);}

if(skipcv==0)	//measure accuracy
{
sumsq=0;for(i=num_train;i<num_samples_use;i++){sumsq+=pow(residuals[i],2);}
value=sumsq/num_test;

if(p==0||value<minmse)	//save values
{
best=p;
minmse=value;
for(j=0;j<data_length;j++){effs2[j]=effs[j];}
}
}

if(enalpha!=0){printf("Performed %d + %d iterations, current model contains %d predictors", count, count2, found);}
else{printf("Performed %d iterations, currrent model contains %d predictors", count, found);}
if(skipcv==0){printf(" (mean squared error %.4f)", value);}
printf("\n\n");

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%.4f\t%d\t", p+1, tryhers[p], found);
if(enalpha!=0){fprintf(output, "%d\t%d\t%d\t", total2, count, count2);}
else{fprintf(output, "NA\t%d\tNA\t", count);}
if(cflag==1){fprintf(output, "YES\t%.6f\t", diff);}
else{fprintf(output, "NO\t%.6f\t", diff);}
if(skipcv==0){fprintf(output, "%.6f\n", value);}
else{fprintf(output, "NA\n");}
fclose(output);

//save
sprintf(filename2,"%s.effects.%d",outfile, p+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<num_tops;j++)
{fprintf(output2, "%s %c %c %.6f %.6f\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j]);}
for(j=0;j<data_length;j++)
{
if(effs[j]!=0){fprintf(output2, "%s %c %c %.6f %.6f\n", preds[j], al1[j], al2[j], centres[j], effs[j]);}
}
fclose(output2);
}	//end of p loop

////////

if(skipcv==0)	//now work out best model
{
printf("Solving for best fitting model (heritability %.6f, mean squared error %.4f))\n", tryhers[best], minmse);
for(j=0;j<data_length;j++){effs[j]=effs2[j];}
set_lambdas(lambdas, exps, data_length, order, tryhers[best], varphen, enalpha);
for(i=0;i<num_samples_use;i++){residuals[i]=Y2[i];}
for(j=0;j<data_length;j++){update_residuals(j, effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);}

//recompute XTX (using all samples)
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
datasqs[j]=0;
if(dtype==1)	//bed format - actual values stored (3 means missing)
{
for(i=0;i<num_samples_use;i++)
{
readint=(int)((data_char[j][i/4] >> (2*(i%4))) & 3);
if(readint!=3){datasqs[j]+=pow(readint-centres[j],2);}
}
}
else	//short speed format - values are speedstart + value * speedscale
{
for(i=0;i<num_samples_use;i++)
{
readint=(int)data_char[j][i];
if(readint!=255){datasqs[j]+=pow(speedstarts[j]+speedscales[j]*readint-centres[j],2);}
}
}
}}

if(enalpha!=0)
{
//set lambdas2 to heritability prior to best and solve (approximate solution will suffice)
if(best!=0)
{
set_lambdas(lambdas2, exps, data_length, order, tryhers[best-1], varphen, enalpha);

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas2[j]*enalpha*fabs(effs[j])+.5*lambdas2[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[best-1]);

if(total>0)
{
printf("Updating effects for %d predictors that have been in at least one model\n", total);
count3=0;
while(1)
{
count3++;
for(j=0;j<data_length;j++)
{
if(evers[j]==1)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
value=0;
if(sum>enalpha*lambdas2[j]){value=(sum-enalpha*lambdas2[j])/(datasqs[j]+(1-enalpha)*lambdas2[j]);}
if(-sum>enalpha*lambdas2[j]){value=(sum+enalpha*lambdas2[j])/(datasqs[j]+(1-enalpha)*lambdas2[j]);}

if(value!=effs[j])	//effect size changed, so update residuals
{
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
effs[j]=value;
}
}}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas2[j]*enalpha*fabs(effs[j])+.5*lambdas2[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[best-1]);

if(total==1){break;}
if(count3>1)
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;

if(count3>=20){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count3, diff, tol);break;}
}
}
}	//end of best!=0
else	//seems very unlikely we will have best=0, so hopefully the following rubbish fudge will not matter
{
set_lambdas(lambdas2, exps, data_length, order, tryhers[best]/2, varphen, enalpha);
for(j=0;j<data_length;j++){effs[j]=0;}
for(i=0;i<num_samples_use;i++){residuals[i]=Y2[i];}
count3=0;
}

//get YTdata and strong set
for(j=0;j<data_length;j++)
{
if(exps[j]>0){YTdata[j]=get_YTdata(j, residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);}
else{YTdata[j]=0;}
}

total2=0;
for(j=0;j<data_length;j++)
{
if(evers[j]==1||fabs(YTdata[j])>=enalpha*(2*lambdas[j]-lambdas2[j])){strongs[j]=1;total2++;}
else{strongs[j]=0;}
}

//considered resetting effs to effs2, but seems fine to leave

if(total>0)	//update effects for predictors ever in model
{
printf("Updating effects for %d predictors that have been in at least one model (again)\n", total);
count=0;
while(1)
{
count++;
for(j=0;j<data_length;j++)
{
if(evers[j]==1)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
value=0;
if(sum>enalpha*lambdas[j]){value=(sum-enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}
if(-sum>enalpha*lambdas[j]){value=(sum+enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}

if(value!=effs[j])	//effect size changed, so update residuals
{
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
effs[j]=value;
}
}}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[best]);

if(total==1){break;}
if(count>1)
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;

if(count>=20){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count, diff, tol);break;}
}
}
else{count=0;printf("No predictors in the model\n");}

if(total2>0)	//update effects for predictors in strong set
{
printf("Updating effects for %d predictors in strong set\n", total2);
count2=0;
cflag=0;
while(1)	
{
count2++;
for(j=0;j<data_length;j++)
{
if(strongs[j]==1)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
value=0;
if(sum>enalpha*lambdas[j]){value=(sum-enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}
if(-sum>enalpha*lambdas[j]){value=(sum+enalpha*lambdas[j])/(datasqs[j]+(1-enalpha)*lambdas[j]);}

if(value!=effs[j])	//effect size changed, so update residuals
{
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
effs[j]=value;
}
}}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[best]);

if(total2==1){diff=0;break;}
if(count2>1)
{
diff=like-likeold;
if(fabs(diff)<tol){cflag=1;break;}
}
likeold=like;

if(count2>=maxiter){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count2, diff, tol);break;}
}
}
else{count2=0;diff=0;cflag=1;}
}	//end of enalpha!=0
else	//ridge - using posterior means
{
//update effects for all non-trivial predictors
count=0;
cflag=0;
while(1)	
{
count++;
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
sum=effs[j]*datasqs[j]+get_YTdata(j, residuals, dtype, data_char, num_train, speedstarts, speedscales, centres);
value=sum/(datasqs[j]+lambdas[j]);
update_residuals(j, value-effs[j], residuals, dtype, data_char, num_samples_use, speedstarts, speedscales, centres);
effs[j]=value;
}
}

//get pen and (partial) likelihood
pen=0;for(j=0;j<data_length;j++){pen+=lambdas[j]*enalpha*fabs(effs[j])+.5*lambdas[j]*(1-enalpha)*pow(effs[j],2);}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
like=-(.5*sumsq+pen)/varphen/(1-tryhers[best]);

if(count>1)
{
diff=like-likeold;
if(fabs(diff)<tol){cflag=1;break;}
}
likeold=like;

if(count>=maxiter){printf("Warning, coordinate descent did not converge after %d iterations (diff in likelihood %f, tolerance %f)\n", count, diff, tol);break;}
}
}

//get size of model
found=0;for(j=0;j<data_length;j++){found+=(effs[j]!=0);}

if(enalpha!=0){printf("Performed %d + %d + %d iterations, final model contains %d predictors\n\n", count3, count, count2, found);}
else{printf("Performed %d iterations, final model contains %d predictors\n\n", count, found);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Final\t%.4f\t%d\t", tryhers[best], found);
if(enalpha!=0){fprintf(output, "%d\t%d\t%d\t", total2, count, count2);}
else{fprintf(output, "NA\t%d\tNA\t", count);}
if(cflag==1){fprintf(output, "YES\t%.6f\tNA\n", diff);}
else{fprintf(output, "NO\t%.6f\tNA\n", diff);}
fclose(output);

//save
sprintf(filename2,"%s.effects.final",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<num_tops;j++)
{fprintf(output2, "%s %c %c %.6f %.6f\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j]);}
for(j=0;j<data_length;j++)
{
if(effs[j]!=0){fprintf(output2, "%s %c %c %.6f %.6f\n", preds[j], al1[j], al2[j], centres[j], effs[j]);}
}
fclose(output2);
}	//end of skipcv=0

//can save coefficients
sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Component Effect SD P\n");
fprintf(output3, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);} 
for(j=0;j<num_envs;j++){fprintf(output3, "Enviromental_%d %.6f %.6f %.4e\n",j, thetas[num_covars+j], thetasds[num_covars+j], thetapvas[num_covars+j]);}
fclose(output3);

if(num_hers==1){printf("Model saved in %s.effects.1\n", outfile);}
else{printf("Models saved in %s.effects.1 ... %s.effects.%d\n", outfile, outfile, num_hers);}
if(skipcv==0){printf("Model with best accuracy saved in %s.effects.final\n", outfile);}
printf("\n");

free(keepboth);free(keepboth2);
free(tryhers);
for(j=0;j<data_length;j++){free(data_char[j]);}free(data_char);
if(dtype==4){free(speedstarts);free(speedscales);}
free(Y);free(Y2);free(YTdata);free(Z);free(datasqs);
free(thetas);free(thetasds);free(thetapvas);
free(exps);free(order);
free(lambdas);free(lambdas2);free(strongs);free(evers);free(effs);free(residuals);
if(skipcv==0){free(effs2);}

///////////////////////////

/*
//get likelihood constant (will get final likelihood by adding on -.5 RSS/e -P/e, where e=varphen(1-her)
var=varphen*(1-tryhers[p]);
if(enalpha!=1){value2=enalpha*pow(2*(1-enalpha)*var,-.5);}
sum=0;sum2=0;value=0;total=0;
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
sum+=log(lambdas[j]);
sum2+=lambdas[j];
if(enalpha!=1)	//want log(erfc(x)), but more stable to get log(erfcx(x))-x^2
{value+=log(my_erfcxf(value2*pow(lambdas[j],.5)))-pow(value2,2)*lambdas[j];}
total++;
}
}

if(enalpha==0)	//ridge
{likecon=-.5*num_samples_use*log(2*M_PI*var)+.5*sum-.5*total*log(2*M_PI*var);}
if(enalpha==1)	//lasso
{likecon=-.5*num_samples_use*log(2*M_PI*var)+sum-total*log(2*var);}
if(enalpha>0&&enalpha<1)	//elastic
{likecon=-.5*num_samples_use*log(2*M_PI*var)+.5*sum-.5*pow(enalpha,2)/(1-enalpha)/var*sum2-.5*total*log(2*M_PI*var/(1-enalpha))-value;}
*/

