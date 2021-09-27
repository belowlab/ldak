/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Perform ridge, bolt and bayesr

///////////////////////////

//set cross-validation values
if(mode==151)	//ridge - just one component
{num_try=1;}
if(mode==152)	//bolt - small then big component
{
num_try=18;
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
for(p=0;p<num_try;p++){tryps[p]=loads[p%6];tryp2s[p]=1-loads[p%6];tryf2s[p]=loads2[(p/6)];}
free(loads);free(loads2);
}
if(mode==153)	//bayesr
{
if(pointmass==1)	//three gaussians and point mass (can not have all ps zero) - add in ridge
{
if(fullspace==0){num_try=35;}
else{num_try=125;}
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);

loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.001;loads[2]=.005;loads[3]=.01;loads[4]=.02;
count=0;
for(j=0;j<5;j++){
for(j2=0;j2<5;j2++){
for(j3=0;j3<5;j3++){
if(j+j2+j3>0&&((j2>=j&&j3>=j2)||fullspace==1))
{
tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
tryps[count]=0;tryp2s[count]=0;tryp3s[count]=0;tryp4s[count]=1;count++;
free(loads);
}
else	//four gaussians
{
if(fullspace==0){num_try=35;}
else{num_try=125;}
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);

loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.001;loads[2]=.005;loads[3]=.01;loads[4]=.02;
count=0;
for(j=0;j<5;j++){
for(j2=0;j2<5;j2++){
for(j3=0;j3<5;j3++){
if((j2>=j&&j3>=j2)||fullspace==1)
{
tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}
}

////////

if(restart==0)	//normal start - make a blank .save file
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=-1;fwrite(&readint, sizeof(int), 1, output2);
readint=0;fwrite(&readint, sizeof(int), 1, output2);
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);

//convenient to set restage and recount
restage=-1;recount=0;
}
else	//restarting from previous run
{
//set readstring, for now and later
if(mode==151){sprintf(readstring,"\"--ridge\"");}
if(mode==152){sprintf(readstring,"\"--bolt\"");}
if(mode==153){sprintf(readstring,"\"--bayesr\"");}

printf("Reading details from previous run\n");
//check progress file exists
sprintf(filename,"%s.progress", outfile);
if(just_check(filename)!=0)
{printf("Error reading %s; make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename, readstring);exit(1);}

//now check .save file exists and get restage and recount
sprintf(filename2,"%s.save.bin", outfile);
if(just_check(filename2)!=0)
{printf("Error reading %s; make sure the output filename is the same as when you first ran %s\n\n", filename2, readstring);exit(1);}

if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&restage, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&recount, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(restage==-1&&recount==0){printf("Error, the previous run had not completed at least ten iterations (if running on a cluster, make sure to allocate sufficient memory and time)\n\n");exit(1);}
if(restage!=-1&&skipcv==1){printf("Error, the previous run performed cross-validation; make sure the cross-validation options are the same as when you first ran %s\n\n", readstring);exit(1);}

//check size
if(restage==-1)	//solving for num_try models
{
fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use)*num_try)
{printf("Error, %s should have size %jd (not %jd); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, (off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use)*num_try, ftello(input2), readstring);exit(1);}
}
else	//solving for a single model
{
fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use))
{printf("Error, %s should have size %jd (not %jd); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, (off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use), ftello(input2), readstring);exit(1);}
}
fclose(input2);

if(restage!=-1)	//check mse file exists
{
sprintf(filename2,"%s.mse", outfile);
if(just_check(filename2)!=0)
{printf("Error reading %s; make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, readstring);exit(1);}
}

printf("Will be restarting from Iteration %d\n\n", recount+1);
}

////////

//set num_train, num_test and keepboths (index training then test samples, relative to ids, then allids)
keepboth=malloc(sizeof(int)*num_samples_use);
keepboth2=malloc(sizeof(int)*num_samples_use);

if(skipcv==0)	//sort cv samples
{
if(cvprop!=-9999)
{
num_test=cvprop*num_samples_use;
num_train=num_samples_use-num_test;

if(restart==0)	//pick samples at random then save
{
for(i=0;i<num_samples_use;i++){keepboth[i]=i;}
permute_int(keepboth,num_samples_use);

sprintf(filename2,"%s.cv.samples",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(i=num_train;i<num_samples_use;i++){fprintf(output2,"%s %s\n", ids1[keepboth[i]], ids2[keepboth[i]]);}
fclose(output2);

sprintf(filename2,"%s.cv.index",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(i=0;i<num_samples_use;i++){fprintf(output2,"%d\n", keepboth[i]);}
fclose(output2);
}
else	//read samples used in first run - will have already set readstring
{
sprintf(filename2,"%s.cv.index", outfile);
if(just_check(filename2)!=0)
{printf("Error reading %s; make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, readstring);exit(1);}
if(countrows(filename2)!=num_samples_use)
{printf("Error, %s should have %d rows (not %d); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, num_test, countrows(filename2), readstring);exit(1);}

if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{
if(fscanf(input2,"%d ", keepboth+i)!=1)
{printf("Error reading row %d of %s\n\n", i+1, filename2);exit(1);}
}
fclose(input2);
}
}	//end of cvprop!=-9999

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
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);free(usedids);
}

printf("Will be using %d samples to train and %d to test\n\n", num_train, num_test);
if(num_train<3){printf("Error, unable to continue with fewer than three training samples\n\n");exit(1);}
if(num_test<3){printf("Error, unable to continue with fewer than three test samples\n\n");exit(1);}
}
else
{
num_train=num_samples_use;
num_test=0;
for(i=0;i<num_samples_use;i++){keepboth[i]=i;}
}

for(i=0;i<num_samples_use;i++){keepboth2[i]=keepsamps[keepboth[i]];}

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
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
datasqs=malloc(sizeof(double)*data_length);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

lambdas=malloc(sizeof(double)*data_length*num_try);
lambdas2=malloc(sizeof(double)*data_length*num_try);
lambdas3=malloc(sizeof(double)*data_length*num_try);
lambdas4=malloc(sizeof(double)*data_length*num_try);

effs=malloc(sizeof(double)*data_length*num_try);
residuals=malloc(sizeof(double)*num_samples_use*num_try);

exps=malloc(sizeof(double)*data_length);

pens=malloc(sizeof(double)*num_try);
likes=malloc(sizeof(double)*num_try);
likesold=malloc(sizeof(double)*num_try);

//variables for speeding up code
anal_warn(bitsize,(num_samples_use+data_length));
data=malloc(sizeof(double)*(num_samples_use*bitsize+4));
cors=malloc(sizeof(double)*bitsize*data_length);
YTdata=malloc(sizeof(double)*bitsize*num_try);
changes=malloc(sizeof(double)*bitsize*num_try);
if(dtype==1)
{
bytelookup=malloc(sizeof(double*)*data_length);
for(j=0;j<data_length;j++){bytelookup[j]=malloc(sizeof(double)*256*4);}
}

//will be reading data in bits
bittotal=(data_length-1)/bitsize+1;

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

if(dtype==1)	//get look up for each predictor
{
printf("Constructing lookup table for each predictor\n\n");
for(j=0;j<data_length;j++)
{
for(k=0;k<256;k++)
{
readint=(k & 3);
if(readint!=3){bytelookup[j][k*4]=readint-centres[j];}
else{bytelookup[j][k*4]=0;}
readint=((k >> 2) & 3);
if(readint!=3){bytelookup[j][k*4+1]=readint-centres[j];}
else{bytelookup[j][k*4+1]=0;}
readint=((k >> 4) & 3);
if(readint!=3){bytelookup[j][k*4+2]=readint-centres[j];}
else{bytelookup[j][k*4+2]=0;}
readint=((k >> 6) & 3);
if(readint!=3){bytelookup[j][k*4+3]=readint-centres[j];}
else{bytelookup[j][k*4+3]=0;}
}
}
}

////////

//get expected squared effect size of each predictor, divided by (varphen * hers)

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

//set lambdas
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++){lambdas[j+p*data_length]=0;lambdas2[j+p*data_length]=0;lambdas3[j+p*data_length]=0;lambdas4[j+p*data_length]=0;}

if(mode==151)	//ridge
{set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, varphen, her, NULL, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, 3);}
if(mode==152)	//bolt
{set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, varphen, her, NULL, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, 4);}
if(mode==153)	//bayesr
{set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, varphen, her, NULL, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, 5+(pointmass==0));}
}

////////

if(restage==-1)	//first solve using only training samples (although use all samples when reading data and computing residuals)
{
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

if(restart==0)	//start with effect sizes at zero
{
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++){effs[j+p*data_length]=0;}
for(i=0;i<num_samples_use;i++){residuals[i+p*num_samples_use]=Y2[i];}
}

//null likelihood is same for all p (computed using training individuals)
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i],2);}
for(p=0;p<num_try;p++){likes[p]=-.5*num_train*log(2*M_PI*(1-her)*varphen)-.5*sumsq/(1-her)/varphen;}
}
else	//read effect sizes, residuals and likelihoods from file
{
sprintf(filename2,"%s.save.bin", outfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(fread(likes, sizeof(double), num_try, input2)!=num_try)
{printf("Error reading model likelihoods from %s\n\n", filename2);exit(1);}
for(p=0;p<num_try;p++)
{
if(fread(effs+p*data_length, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading effects sizes for Model %d from %s\n\n", p+1, filename2);exit(1);}
if(fread(residuals+p*num_samples_use, sizeof(double), num_samples_use, input2)!=num_samples_use)
{printf("Error reading residuals for Model %d from %s\n\n", p+1, filename2);exit(1);}
}
fclose(input2);
}

//get predictor-predictor correlations
printf("Computing predictor-predictor correlations\n\n");
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_train, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);
}

//screen and file print (if restart=1, have set filename and checked progress file exists)
if(skipcv==0)	//doing cv - so must be bolt or bayesr and have multiple models
{printf("Estimating effect sizes for %d models using training samples\nIter\tNum_Con\tMax_Diff\tTolerance\n", num_try);}
else	//not doing cv - probably ridge, but can be bolt or bayesr
{
if(num_try==1){printf("Estimating effect sizes\nIter\tDiff\tTolerance\n");}
else{printf("Estimating effect sizes for %d models\nIter\tNum_Con\tMax_Diff\tTolerance\n", num_try);}
}

if(restart==0)
{
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
if(skipcv==0){fprintf(output,"Solving using training samples\n");}
if(num_try==1){fprintf(output,"Iter\tDifference\tTolerance\n");}
else{fprintf(output,"Iter\tNum_Converged\tMax_Difference\tTolerance\n");}
fclose(output);
}

//iterate effect sizes
count=recount;
while(1)
{
count++;

for(p=0;p<num_try;p++){pens[p]=0;}	//use this to keep track of likelihood terms
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//get t(X) residuals
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_try, &num_train, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitsize);

#pragma omp parallel for private(p,j,j2,sum,postmean) schedule(dynamic)
for(p=0;p<num_try;p++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XTresiduals
sum=effs[j+p*data_length]*datasqs[j]+YTdata[(j-bitstart)+p*bitsize];

if(mode==151)	//ridge
{postmean=get_postmean(sum, lambdas[j+p*data_length], -9999, -9999, -9999, datasqs[j], (1-her)*varphen, -9999, -9999, -9999, -9999, pens+p, 3);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], -9999, -9999, datasqs[j], (1-her)*varphen, tryps[p], tryp2s[p], -9999, -9999, pens+p, 4);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], (1-her)*varphen, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], pens+p,  5+(pointmass==0));}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitsize]=postmean-effs[j+p*data_length];
effs[j+p*data_length]=postmean;
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitsize]-=changes[j-bitstart+p*bitsize]*cors[(size_t)j*bitsize+j2-bitstart];}
}

postmean=0;
}	//end of j loop
}	//end of p loop

//update residuals (for all samples)
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_try, &bitlength, &alpha, data, &num_samples_use, changes, &bitsize, &beta, residuals, &num_samples_use);
}	//end of bit loop

for(p=0;p<num_try;p++)	//approx likelihoods
{
likesold[p]=likes[p];
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
likes[p]=-.5*num_train*log(2*M_PI*(1-her)*varphen)-.5*sumsq/(1-her)/varphen-pens[p];
}

if(count%10==0)	//save in case need to restart
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=-1;fwrite(&readint, sizeof(int), 1, output2);
readint=count;fwrite(&readint, sizeof(int), 1, output2);
fwrite(likes, sizeof(double), num_try, output2);
for(p=0;p<num_try;p++)
{
fwrite(effs+p*data_length, sizeof(double), data_length, output2);
fwrite(residuals+p*num_samples_use, sizeof(double), num_samples_use, output2);
}
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);
}

//print update (require number converged and largest difference in likelihood)
cflag=0;for(p=0;p<num_try;p++){cflag+=(fabs(likes[p]-likesold[p])<tol);}
for(p=0;p<num_try;p++)
{
if(p==0){diff=fabs(likes[p]-likesold[p]);}
if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

if(num_try==1)
{
if(count==1){printf("%d\tNA\t%.4f\n", count, tol);}
else{printf("%d\t%.4f\t%.4f\n", count, diff, tol);}
}
else
{
if(count==1){printf("%d\t%d\tNA\t%.4f\n", count, cflag, tol);}
else{printf("%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);}
}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(num_try==1)
{
if(count==1){fprintf(output,"%d\tNA\t%.4f\n", count, tol);}
else{fprintf(output,"%d\t%.4f\t%.4f\n", count, diff, tol);}
}
else
{
if(count==1){fprintf(output,"aa%d\t%d\tNA\t%.4f\n", count, cflag, tol);}
else{fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);}
}
fclose(output);

if(diff<tol){break;}
if(count>=maxiter){printf("\nWarning, Variational Bayes did not converge after %d iterations (%d out of %d converged, largest difference in likelihood %f, tolerance %f); consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance\n", count, cflag, num_try, diff, tol);break;}
}	//end of while loop
printf("\n");

for(p=0;p<num_try;p++)
{
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
likes[p]=-.5*num_train*log(2*M_PI*(1-her)*varphen)-.5*sumsq/(1-her)/varphen-pens[p];
}

if(skipcv==1)	//have final models, so save
{
if(num_try==1)
{
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<num_tops;j++)
{fprintf(output2, "%s %c %c %.6f %.4e\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j]);}
for(j=0;j<data_length;j++)
{
if(exps[j]>0){fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[j]);}
}
fclose(output2);
printf("Model saved in %s\n\n", filename2);
}
else
{
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre");
for(p=0;p<num_try;p++){fprintf(output2," Model%d", p+1);}
fprintf(output2,"\n");
for(j=0;j<num_tops;j++)
{
fprintf(output2, "%s %c %c %.6f", tpreds[j], tal1[j], tal2[j], tcentres[j]);
for(p=0;p<num_try;p++){fprintf(output2," %.4e", thetas[num_covars+j]);}
fprintf(output2,"\n");
}
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=0;p<num_try;p++){fprintf(output2," %.4e", effs[j+p*data_length]);}
fprintf(output2,"\n");
}
}
fclose(output2);
printf("All models saved in %s\n\n", filename2);
}
}
}	//end of restart=-1

////////

if(skipcv==0)	//now focus on best-fitting model (must be bolt or bayesr)
{
if(restage==-1)	//get best model based on training samples
{
//measure and print accuracy, recording best
printf("Measuring accuracy of each model\n");
sprintf(filename2,"%s.mse",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(mode==152){fprintf(output2,"Model\tp\tf2\tMean_Squared_Error\tApprox_Likelihood\n");}
else{fprintf(output2,"Model\tp1\tp2\tp3\tMean_Squared_Error\tApprox_Likelihood\n");}

for(p=0;p<num_try;p++)
{
sumsq=0;for(i=num_train;i<num_samples_use;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
value=sumsq/num_test;

if(mode==152)
{printf("Model %d: p %.4f, f2 %.4f, mean squared error %.4f\n", p+1, tryps[p], tryf2s[p], value);
fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryps[p], tryf2s[p], value, likes[p]);}
else
{printf("Model %d: p1 %.4f, p2 %.4f, p3 %.4f, mean squared error %.4f, approx likelihood %.4f\n", p+1, tryp4s[p], tryp3s[p], tryp2s[p], value, likes[p]);
fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryp4s[p], tryp3s[p], tryp2s[p], value, likes[p]);}

if(p==0){minmse=value;best=0;}
if(value<minmse){minmse=value;best=p;}
}
printf("\n");
fclose(output2);
}
else	//already know best, can read accuracy from file
{
best=restage;
sprintf(filename2,"%s.mse", outfile);
read_values(filename2,&minmse,1,NULL,5,1+best,0);
}

//now solve for best model, using all samples

//can extract XTX for each predictor from sqdevs
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*num_samples_use;}

if(restage==-1)	//already have effects and residuals, so just renew likelihood (will be approximate, as using existing pen)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(residuals[i+best*num_samples_use],2);}
likes[best]=-.5*num_samples_use*log(2*M_PI*(1-her)*varphen)-.5*sumsq/(1-her)/varphen-pens[best];
}
else	//read effect sizes, residuals and likelihood from file
{
sprintf(filename2,"%s.save.bin", outfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(fread(likes+best, sizeof(double), 1, input2)!=1)
{printf("Error reading likelihood for Model %d from %s\n\n", best+1, filename2);exit(1);}
if(fread(effs+best*data_length, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading effects sizes for Model %d from %s\n\n", best+1, filename2);exit(1);}
if(fread(residuals+best*num_samples_use, sizeof(double), num_samples_use, input2)!=num_samples_use)
{printf("Error reading residuals for Model %d from %s\n\n", best+1, filename2);exit(1);}

fclose(input2);
}

//get predictor-predictor correlations
printf("Computing predictor-predictor correlations\n\n");
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);
}

//screen and file print
if(mode==152){printf("Estimating effect sizes for best-fitting model (p %.4f, f2 %.4f, mean squared error %.4f)\n", tryps[best], tryf2s[best], minmse);}
else{printf("Estimating effect sizes for best-fitting model (p1 %.4f, p2 %.4f, p3 %.4f, mean squared error %.4f)\n", tryp4s[best], tryp3s[best], tryp2s[best], minmse);}

printf("Iter\tNum_Con\tDiff\tTolerance\n");
if(restage==-1)
{
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Solving best-fitting model using all samples\n");
fprintf(output,"Iter\tConverged\tDifference\tTolerance\n");
fclose(output);
}

//iterate effect sizes
if(restage==-1){count=0;}
else{count=recount;}
while(1)
{
count++;

pens[best]=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//get t(X) residuals
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, residuals+best*num_samples_use, &one, &beta, YTdata+best*bitsize, &one);

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XTresiduals
sum=effs[j+best*data_length]*datasqs[j]+YTdata[(j-bitstart)+best*bitsize];

if(mode==151)	//ridge
{postmean=get_postmean(sum, lambdas[j+best*data_length], -9999, -9999, -9999, datasqs[j], (1-her)*varphen, -9999, -9999, -9999, -9999, pens+best, 3);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[j+best*data_length], lambdas2[j+best*data_length], -9999, -9999, datasqs[j], (1-her)*varphen, tryps[best], tryp2s[best], -9999, -9999, pens+best, 4);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[j+best*data_length], lambdas2[j+best*data_length], lambdas3[j+best*data_length], lambdas4[j+best*data_length], datasqs[j], (1-her)*varphen, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], pens+best, 5+(pointmass==0));}

//save, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+best*bitsize]=postmean-effs[j+best*data_length];
effs[j+best*data_length]=postmean;
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+best*bitsize]-=changes[j-bitstart+best*bitsize]*cors[(size_t)j*bitsize+j2-bitstart];}
}
}	//end of j loop

//update residuals
alpha=-1.0;beta=1.0;
dgemv_("N", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, changes+best*bitsize, &one, &beta, residuals+best*num_samples_use, &one);
}	//end of bit loop

//approx likelihood
likesold[best]=likes[best];
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(residuals[i+best*num_samples_use],2);}
likes[best]=-.5*num_samples_use*log(2*M_PI*(1-her)*varphen)-.5*sumsq/(1-her)/varphen-pens[best];

if(count%10==0)	//save in case need to restart
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=best;fwrite(&readint, sizeof(int), 1, output2);
readint=count;fwrite(&readint, sizeof(int), 1, output2);
fwrite(likes+best, sizeof(double), 1, output2);
fwrite(effs+best*data_length, sizeof(double), data_length, output2);
fwrite(residuals+best*num_samples_use, sizeof(double), num_samples_use, output2);
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);
}

//see whether converged and get difference in likelihood
cflag=(fabs(likes[best]-likesold[best])<tol);
diff=fabs(likes[best]-likesold[best]);

printf("%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);
fclose(output);

if(diff<tol){break;}
if(count>=maxiter){printf("\nWarning, Variational Bayes did not converge after %d iterations (difference in likelihood %f, tolerance %f); consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance\n", count, diff, tol);break;}
}	//end of while loop
printf("\n");

//save model
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<num_tops;j++)
{fprintf(output2, "%s %c %c %.6f %.4e\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j]);}
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[j+best*data_length]);}
}
fclose(output2);

printf("Best-fitting model saved in %s\n\n", filename2);
}	//end of skipcv=0

////////

if(mode==152){free(tryps);free(tryp2s);free(tryf2s);}
if(mode==153){free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);}
free(keepboth);free(keepboth2);
for(j=0;j<data_length;j++){free(data_char[j]);}free(data_char);
if(dtype==4){free(speedstarts);free(speedscales);}
free(Y);free(Y2);free(Z);free(datasqs);
free(thetas);free(thetasds);free(thetapvas);
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(effs);free(residuals);
free(exps);
free(pens);free(likes);free(likesold);
free(data);free(cors);free(YTdata);free(changes);
if(dtype==1){for(j=0;j<data_length;j++){free(bytelookup[j]);}free(bytelookup);}

///////////////////////////

