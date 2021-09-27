/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Approximate he and pcgc

///////////////////////////

//get the heritability model
model_warn(data_length*3/2, num_parts+1);
pindexes=malloc(sizeof(int *)*(num_parts+1));
pweights=malloc(sizeof(double *)*(num_parts+1));

addpart=get_her_model(num_parts, partpref, pindexes, pweights, data_length, keeppreds_use, num_preds, allpreds, predorder, parttype, backpart, allone);

if(addpart==2)	//were some redundant predictors, so squeeze down, and reset addpart
{
count=0;
for(j=0;j<data_length;j++)
{
flag=0;for(q=0;q<num_parts;q++){flag+=pindexes[q][j];}
if(flag>0)
{
if(count!=j)
{
chr[count]=chr[j];cmbp[count]=cmbp[j];al1[count]=al1[j];al2[count]=al2[j];
free(preds[count]);copy_string(preds,count,preds[j]);
keeppreds_use[count]=keeppreds_use[j];
weights[count]=weights[j];pvalues[count]=pvalues[j];
for(q=0;q<num_parts;q++){pindexes[q][count]=pindexes[q][j];pweights[q][count]=pweights[q][j];}
}
count++;
}}
for(j=count;j<data_length;j++){free(preds[j]);}
data_length=count;

addpart=0;
}

//find start and end of jackknife blocks
bstarts=malloc(sizeof(int)*num_blocks);
bends=malloc(sizeof(int)*num_blocks);
for(p=0;p<num_blocks;p++)
{bstarts[p]=(double)p/num_blocks*data_length;bends[p]=(double)(p+1)/num_blocks*data_length;}

//check each kinship matrix not entirely contained in a block
for(q=0;q<num_parts+addpart;q++)
{
count=0;
for(j=0;j<data_length;j++){count+=(pweights[q][j]!=0);}
for(p=0;p<num_blocks;p++)
{
count2=0;
for(j=bstarts[p];j<bends[p];j++){count2+=(pweights[q][j]!=0);}
if(count2==count){printf("Error, all %d predictors that contribute to Kinship Matrix %d are contained within Jackknife Block %d; most likely, you should remove this kinship matrix from the model\n\n", count, q+1, p+1);exit(1);}
}
}

////////

//set total, total3 total3, and bitsize
total=num_parts+addpart;
total2=num_vects+1;
total3=(num_vects+1)*(num_parts+addpart);
if(bitsize>data_length){bitsize=data_length;}

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

value=(double)num_samples_use/1024/1024/1024*8*(2*total+total2+2*total3)+(double)bitsize/1024/1024/1024*8*(total2+total3);
if(value>1){printf("Warning, to perform the analysis will require approximately %.1f Gb\n\n", value);}

kinsums=malloc(sizeof(double)*total);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Y2=malloc(sizeof(double)*num_samples_use);

R=malloc(sizeof(double)*num_samples_use*total2);
RTdata=malloc(sizeof(double)*total2*bitsize);
RTdata2=malloc(sizeof(double)*total3*bitsize);

mkinsD=malloc(sizeof(double)*num_samples_use*total);
mkinsD2=malloc(sizeof(double)*num_samples_use*total);
mkinsR=malloc(sizeof(double)*num_samples_use*total3);
mkinsR2=malloc(sizeof(double)*num_samples_use*total3);

KKtraces=malloc(sizeof(double)*total*total);
KKtraces2=malloc(sizeof(double)*total*total);
KYtraces=malloc(sizeof(double)*total);
KYtraces2=malloc(sizeof(double)*total);

hers=malloc(sizeof(double)*(total+1));
hersds=malloc(sizeof(double)*(total+1));
shares=malloc(sizeof(double)*total);
sharesds=malloc(sizeof(double)*total);
cohers=malloc(sizeof(double)*total*total);
jacks=malloc(sizeof(double)*(total+1+total)*num_blocks);

////////

//fill Y and Z
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

//solve covariates (get covher and topher) and fill Y2
if(mode==126)	//linear model - Y2 contains standardized residuals
{reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Y2, 1, &covher, &topher);}
else	//logistic model - Y2 contains pcgc residuals
{reg_covar_log(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Y2, &covher, &topher, prev, 0.0001, 20);}

if(num_covars>1){printf("Proportion of variance explained by the %d covariates: %.4f\n", num_covars, covher);}
if(num_tops==1){printf("Proportion of variance explained by the top predictor: %.4f\n", topher);}
if(num_tops>1){printf("Proportion of variance explained by the %d top predictor: %.4f\n", num_tops, topher);}
if(num_fixed>1){printf("\n");}

//fill start of R with random values, then end with adjusted phenotypes
for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){R[i+g*num_samples_use]=RNOR;}
}
for(i=0;i<num_samples_use;i++){R[i+num_vects*num_samples_use]=Y2[i];}

//set mkinsD, mkinsD2, mkinsR and mkinsR2 to zero
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++){mkinsD[(size_t)q*num_samples_use+i]=0;mkinsD2[(size_t)q*num_samples_use+i]=0;}
}
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++){mkinsR[(size_t)k*num_samples_use+i]=0;mkinsR2[(size_t)k*num_samples_use+i]=0;}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//open temp file
sprintf(filename2,"%s.temp", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

////////

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating traces for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating traces for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, power, strcmp(centresfile,"blank")!=0, hwestand, weights+bitstart, 1, preds+bitstart);

if(num_fixed>1)	//adjust data for covariates
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

//compute t(R) data
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &bitlength, &num_samples_use, &alpha, R, &num_samples_use, data, &num_samples_use, &beta, RTdata, &total2);

//load RTdata2, which has RTdata Wk for each vector then response, for kinship 1, then kinship 2, etc
#pragma omp parallel for private(j,q,g) schedule (static)
for(j=0;j<bitlength;j++)
{
for(q=0;q<total;q++)
{
for(g=0;g<total2;g++){RTdata2[q*total2+g+j*total3]=RTdata[g+j*total2]*pweights[q][bitstart+j];}
}
}

for(p=0;p<num_blocks;p++)	//loop through blocks, to find one(s) being used
{
if(bstarts[p]>=bitend){break;}

if(bends[p]>bitstart)	//some predictors within this block
{
start=bstarts[p];
if(start<bitstart){start=bitstart;}
end=bends[p];
if(end>bitend){end=bitend;}

//get contribution to diagonal terms
for(j=start;j<end;j++)
{
for(q=0;q<total;q++)
{
if(pweights[q][j]!=0)
{
for(i=0;i<num_samples_use;i++)
{mkinsD2[(size_t)q*num_samples_use+i]+=pow(data[(size_t)(j-bitstart)*num_samples_use+i],2)*pweights[q][j];}
}}
}

//get contribution to traces
if(parttype==0)	//multiply for all kinships at once
{
token=end-start;
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &total3, &token, &alpha, data+(size_t)(start-bitstart)*num_samples_use, &num_samples_use, RTdata2+(start-bitstart)*total3, &total3, &beta, mkinsR2, &num_samples_use);
}
else	//only multiply for kinships that contain predictors
{
for(q=0;q<total;q++)
{
count=0;for(j=start;j<end;j++){count+=(pweights[q][j]!=0);}
if(count>0)
{
token=end-start;
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &total2, &token, &alpha, data+(size_t)(start-bitstart)*num_samples_use, &num_samples_use, RTdata2+(start-bitstart)*total3+q*total2, &total3, &beta, mkinsR2+(size_t)q*total2*num_samples_use, &num_samples_use);
}
}
}

if(bends[p]<=bitend)	//finished with this block
{
//save mkinsD2 and mkinsR2
for(q=0;q<total;q++){fwrite(mkinsD2+(size_t)q*num_samples_use, sizeof(double), num_samples_use, output2);}
for(k=0;k<total3;k++){fwrite(mkinsR2+(size_t)k*num_samples_use, sizeof(double), num_samples_use, output2);}

//add contribution of mkinsD2 to mkinsD then set to zero
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++)
{mkinsD[(size_t)q*num_samples_use+i]+=mkinsD2[(size_t)q*num_samples_use+i];mkinsD2[(size_t)q*num_samples_use+i]=0;}
}

//add contribution of mkinsR2 to mkinsR then set to zero
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++)
{mkinsR[(size_t)k*num_samples_use+i]+=mkinsR2[(size_t)k*num_samples_use+i];mkinsR2[(size_t)k*num_samples_use+i]=0;}
}
}
}	//end of working on block p

}	//end of p loop
}	//end of bit loop
printf("\n");

fclose(output2);

//get the traces of the kinship matrices
for(q=0;q<total;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=mkinsD[(size_t)q*num_samples_use+i];}
kinsums[q]=sum/num_samples_use;
if(kinsums[q]==0){printf("Error, Kinship Matrix %d has trace zero, so it is not possible to continue\n\n", q+1);exit(1);}
}

////////

//ready to estimate heritabilities, first for all predictors, then excluding a block

//open temp file for reading
if((output2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//will save trace products in case wish to repeat
sprintf(filename3,"%s.traces",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

//set KKtraces to minus the cross-product of diagonals
alpha=-1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, mkinsD, &num_samples_use, mkinsD, &num_samples_use, &beta, KKtraces, &total);

//add average contribution from vectors
for(g=0;g<num_vects;g++)
{
token=num_samples_use*total2;
alpha=1.0/num_vects;beta=1.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, mkinsR+(size_t)g*num_samples_use, &token, mkinsR+(size_t)g*num_samples_use, &token, &beta, KKtraces, &total);
}

//KYtraces is phenotype terms from mkinsR x Y, minus mkinsD x Y^2
token=num_samples_use*total2;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &total, &alpha, mkinsR+(size_t)num_vects*num_samples_use, &token, Y2, &one, &beta, KYtraces, &one);
for(q=0;q<total;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=mkinsD[(size_t)q*num_samples_use+i]*pow(Y2[i],2);}
KYtraces[q]-=sum;
}

//divide by traces and by two
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){KKtraces[q+q2*total]=KKtraces[q+q2*total]/kinsums[q]/kinsums[q2]/2;}
KYtraces[q]=KYtraces[q]/kinsums[q]/2;
}

//save traces, and KYtraces in KYtraces2
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){fprintf(output3,"%.6f ", KKtraces[q+q2*total]);}
fprintf(output3,"%.6f\n", KYtraces[q]);
KYtraces2[q]=KYtraces[q];
}

//get estimates
(void)eigen_invert(KKtraces, total, KKtraces2, 1, KYtraces, 1);
sum=0;for(q=0;q<total;q++){sum+=KYtraces[q];}
for(q=0;q<total;q++){hers[q]=KYtraces[q];shares[q]=KYtraces[q]/sum;}
hers[total]=sum;

//get likelihoods and lrts
sum=0;sum2=0;for(i=0;i<num_samples_use;i++){sum+=pow(Y2[i],2);sum2+=pow(Y2[i],4);}
sumsq=(pow(sum,2)-sum2)/2;
smax=(size_t)num_samples_use*(num_samples_use-1)/2;
likenull=-.5*smax*(1+log(2*M_PI*sumsq/smax));
for(q=0;q<total;q++){sumsq-=KYtraces2[q]*KYtraces[q];}
like=-.5*smax*(1+log(2*M_PI*sumsq/smax));
lrtstat=2*(like-likenull);
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);

////////

//get estimates excluding block p

for(p=0;p<num_blocks;p++)
{
if(p%10==0)
{
printf("Processing traces for Jackknife Block %d of %d\n", p+1, num_blocks);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Processing traces for Jackknife Block %d of %d\n", p+1, num_blocks);
fclose(output);
}

//read mkinsD2 and mkinsR2 - will already be at correct point in file
for(q=0;q<total;q++)
{
if(fread(mkinsD2+(size_t)q*num_samples_use, sizeof(double), num_samples_use, output2)!=num_samples_use)
{printf("Error reading values from %s\n\n", filename2);exit(1);}
}
for(k=0;k<total3;k++)
{
if(fread(mkinsR2+(size_t)k*num_samples_use, sizeof(double), num_samples_use, output2)!=num_samples_use)
{printf("Error reading values from %s\n\n", filename2);exit(1);}
}

//subtract mkinsD2 and mkinsR2 from mkinsD and mkinsR
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++){mkinsD[(size_t)q*num_samples_use+i]-=mkinsD2[(size_t)q*num_samples_use+i];}
}
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++){mkinsR[(size_t)k*num_samples_use+i]-=mkinsR2[(size_t)k*num_samples_use+i];}
}

alpha=-1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, mkinsD, &num_samples_use, mkinsD, &num_samples_use, &beta, KKtraces, &total);

for(g=0;g<num_vects;g++)
{
token=num_samples_use*total2;
alpha=1.0/num_vects;beta=1.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, mkinsR+(size_t)g*num_samples_use, &token, mkinsR+(size_t)g*num_samples_use, &token, &beta, KKtraces, &total);
}

token=num_samples_use*total2;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &total, &alpha, mkinsR+(size_t)num_vects*num_samples_use, &token, Y2, &one, &beta, KYtraces, &one);
for(q=0;q<total;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=mkinsD[(size_t)q*num_samples_use+i]*pow(Y2[i],2);}
KYtraces[q]-=sum;
}

for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){KKtraces[q+q2*total]=KKtraces[q+q2*total]/kinsums[q]/kinsums[q2]/2;}
KYtraces[q]=KYtraces[q]/kinsums[q]/2;
}

//save, then get estimates
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){fprintf(output3,"%.6f ", KKtraces[q+q2*total]);}
fprintf(output3,"%.6f\n", KYtraces[q]);
}

(void)eigen_invert(KKtraces, total, KKtraces2, 1, KYtraces, 1);
sum=0;for(q=0;q<total;q++){sum+=KYtraces[q];}
for(q=0;q<total;q++){jacks[q+p*(total+1+total)]=KYtraces[q];}
jacks[total+p*(total+1+total)]=sum;
for(q=0;q<total;q++){jacks[total+1+q+p*(total+1+total)]=KYtraces[q]/sum;}

//add mkinsD2 and mkinsR2 back onto mkinsD and mkinsR
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++){mkinsD[(size_t)q*num_samples_use+i]+=mkinsD2[(size_t)q*num_samples_use+i];}
}
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++){mkinsR[(size_t)k*num_samples_use+i]+=mkinsR2[(size_t)k*num_samples_use+i];}
}
}	//end of p loop
printf("\n");

fclose(output2);
fclose(output3);

//can delete temp file
sprintf(cmd,"rm -rf %s",filename2);
system(cmd);

for(q=0;q<total+1+total;q++)	//get sds
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=jacks[q+p*(total+1+total)];sumsq+=pow(jacks[q+p*(total+1+total)],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
if(q<total+1){hersds[q]=pow(var,.5);}
else{sharesds[q-total-1]=pow(var,.5);}
}

//get covariances
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[q+p*(total+1+total)];
sum2+=jacks[q2+p*(total+1+total)];
sumsq+=jacks[q+p*(total+1+total)]*jacks[q2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers[q+q2*total]=var;
}}

printf("Total heritability (SD): %.4f (%.4f)\n\n", hers[total], hersds[total]);

////////

//adjust for tops
for(q=0;q<total;q++){hers[q]*=(1-topher);hersds[q]*=(1-topher);}

//save stuff

sum=0;for(q=0;q<total;q++){sum+=kinsums[q];}

if(mode==126){sprintf(filename2,"%s.fasthe", outfile);}
else{sprintf(filename2,"%s.fastpcgc", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Top_Predictors %d\nNum_Covariates %d\n", total, num_tops, num_covars);
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", num_samples_use, respcounts[0]);
fprintf(output2, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability Her_SD Size Mega_Intensity Int_SD\n");
for(q=0;q<total;q++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", q+1, hers[q], hersds[q], kinsums[q], hers[q]/kinsums[q]*1000000, hersds[q]/kinsums[q]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);
fclose(output2);

sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(mode==126){fprintf(output3, "Component Effect SD P\n");}
else{fprintf(output3, "Component Log_Odds SD P\n");}
fprintf(output3, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output3);

sprintf(filename4,"%s.share", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SD Expected Enrichment SD\n");
for(q=0;q<total;q++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", q+1, shares[q], sharesds[q], kinsums[q]/sum, shares[q]/kinsums[q]*sum, sharesds[q]/kinsums[q]*sum);}
fclose(output4);

if(mode==126&&prev!=-9999)	//binary he
{
factor=get_factor(Y, num_samples_use, prev, -9999, outfile);

sprintf(filename5,"%s.fasthe.liab", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Top_Predictors %d\nNum_Covariates %d\n", total, num_tops, num_covars);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", num_samples_use, respcounts[0]);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability Her_SD Size Mega_Intensity Int_SD\n");
for(q=0;q<total;q++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", q+1, hers[q]*factor, hersds[q]*factor, kinsums[q], hers[q]/kinsums[q]*1000000*factor, hersds[q]/kinsums[q]*1000000*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);
fclose(output5);
}

if(mode==127)	//save marginals for pcgc
{
sprintf(filename5,"%s.fastpcgc.marginal", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Top_Predictors %d\nNum_Covariates %d\n", total, num_tops, num_covars);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", num_samples_use, respcounts[0]);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability Her_SD Size Mega_Intensity Int_SD\n");
for(q=0;q<total;q++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", q+1, hers[q]*(1-covher), hersds[q]*(1-covher), kinsums[q], hers[q]/kinsums[q]*1000000*(1-covher), hersds[q]/kinsums[q]*1000000*(1-covher));}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher), sum, hers[total]/sum*1000000*(1-covher), hersds[total]/sum*1000000*(1-covher));
fclose(output5);
}

sprintf(filename6,"%s.cross", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(q=0;q<total;q++){fprintf(output6, "Her_K%d\t", q+1);}
fprintf(output6, "\n");
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++)
{fprintf(output6, "%.6f\t", cohers[q+q2*total]);}
fprintf(output6, "\n");
}
fclose(output6);

////////

/*
//save heritability model (first indexes of predictors, then scaled pweights for each partition)
sprintf(filename4,"%s.traces.model", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor Centre Scaling Weight");
for(q=0;q<total;q++){fprintf(output4, " Kinship%d", q+1);}
fprintf(output4, "\n");
for(j=0;j<data_length;j++)
{
fprintf(output4, "%s %.6f %.6f %.6f", preds[j], centres[j], mults[j], weights[j]);
for(q=0;q<total;q++){fprintf(output4, " %.6f", pweights[q][j]);}
fprintf(output4, "\n");
}
fclose(output4);

//finish by making root file
sprintf(filename5,"%s.traces.root", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5,"Datafile %s\n", datafile);
fprintf(output5,"Num_Samples %d\nNum_Predictors %d\n", num_samples, num_preds);
fprintf(output5,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, data_length);
fprintf(output5,"Num_Kinships %d\nNum_Blocks %d\nNum_Vectors %d\n", total, num_blocks, num_vects);
fclose(output5);
*/

for(q=0;q<num_parts+1;q++){free(pindexes[q]);free(pweights[q]);}free(pindexes);free(pweights);
free(bstarts);free(bends);
free(data);
free(kinsums);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);free(Y2);
free(R);free(RTdata);free(RTdata2);
free(mkinsD);free(mkinsD2);free(mkinsR);free(mkinsR2);
free(KKtraces);free(KKtraces2);free(KYtraces);free(KYtraces2);
free(hers);free(hersds);free(shares);free(sharesds);free(cohers);free(jacks);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

