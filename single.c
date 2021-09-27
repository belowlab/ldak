/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor testing

///////////////////////////

if(bitsize>data_length){bitsize=data_length;}

//allocate variables
if(num_kins==1){data_warn2(bitsize,2*num_samples_use);}
else{data_warn2(bitsize,num_samples_use);}
data=malloc(sizeof(double)*num_samples_use*bitsize);

order=malloc(sizeof(int)*num_samples_use);
sweights=malloc(sizeof(double)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*(num_fixed+1));

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4);

if(num_fixed>1)
{
ZTZ=malloc(sizeof(double)*(num_fixed-1));
ZTdata=malloc(sizeof(double)*(num_fixed-1)*bitsize);
}

if(num_kins==1)
{
UTY=malloc(sizeof(double)*num_samples_use);
UTZ=malloc(sizeof(double)*num_samples_use*(num_fixed+1));
UTdata=malloc(sizeof(double)*num_samples_use*bitsize);
vstarts=malloc(sizeof(double)*2);
}

if(mode==131&&num_kins==0&&num_fixed==1){YTdata=malloc(sizeof(double)*bitsize);}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

if(strcmp(sampwfile,"blank")!=0)	//get sample weights
{read_sampwfile(sampwfile, sweights, num_samples_use, ids3, ids1, ids2);}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

//fill Y and start of Z, get varphen, and if using, get factor
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]]*pow(sweights[order[i]],.5);
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use]*sweights[order[i]];}
}

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Y[i];sumsq+=pow(Y[i],2);}
varphen=sumsq/num_samples_use-pow(sum/num_samples_use,2);

if(prev!=-9999){factor=get_factor(Y, num_samples_use, prev, -9999, outfile);}

//set tindex - indicates whether each predictor is a top preds (and if so, which)
for(j=0;j<data_length;j++){tindex[j]=-9999;}
if(num_tops>0)	//find the tops
{
mark=0;
for(j=0;j<num_tops;j++)
{
while(tkeeppreds[j]>keeppreds_use[mark]){mark++;}
tindex[mark]=j;
}
}

if(num_fixed>1)	//get sumsqs for fixed effects (except intercept)
{
for(j=0;j<num_fixed-1;j++)
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++)
{sum+=Z[i+(j+1)*num_samples_use];sumsq+=pow(Z[i+(j+1)*num_samples_use],2);}
ZTZ[j]=sumsq-pow(sum,2)/num_samples_use;
}
}

if(num_kins==1)	//get UTY and start of UTZ
{
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_samples_use, &alpha, U, &num_samples_use, Y, &one, &beta, UTY, &one);
dgemm_("T", "N", &num_samples_use, &num_fixed, &num_samples_use, &alpha, U, &num_samples_use, Z, &num_samples_use, &beta, UTZ, &num_samples_use);
}

//solve null model - get thetas for covariates, and when using a kinship, starting variances
if(mode==131)
{
if(num_kins==0)
{reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, NULL, -9999, NULL, NULL);}
else
{
printf("Solving Null Model\n");
linear_reml(num_samples_use, num_fixed, Y, Z, U, E, UTY, UTZ, kintraces, NULL, vstarts, thetas, thetasds, thetapvas, constrain, tol, maxiter, 1);
printf("\n");
}
}
else{reg_covar_log(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, NULL, NULL, NULL, -9999, tol, maxiter);}

//save
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Effect SD P\n");
fprintf(output, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

////////

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.assoc",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
if(mode==131){fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Effect SD Effect_Liability SD_Liability A1_Mean MAF\n");}
else{fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Log_OR SD A1_Mean MAF\n");}

sprintf(filename3,"%s.summaries",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Direction Stat n\n");

sprintf(filename4,"%s.pvalues",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor P\n");

sprintf(filename5,"%s.score",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor A1 A2 Centre ALL P<0.1 P<0.01 P<0.001 P<0.0001 P<0.00001 P<5e-8\n");

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

fclose(output);
printf("Performing single-SNP analysis for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing single-SNP analysis for Chunk %d of %d\n", bit+1, bittotal);
fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output5);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//read data for chunk, and just centre
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, 0, NULL, 1, preds+bitstart);

if(strcmp(sampwfile,"blank")!=0)	//scale by sample weights
{
for(j=0;j<bitlength;j++)
{
value=pow(sweights[j],.5);
for(i=0;i<num_samples_use;i++){data[i+j*num_samples_use]*=value;}
}
}

if(num_fixed>1)	//get ZTdata (excluding intercept), and see if correlations exceed mincor
{
token=num_fixed-1;
alpha=1.0;beta=0.0;
dgemm_("T", "N", &token, &bitlength, &num_samples_use, &alpha, Z+num_samples_use, &num_samples_use, data, &num_samples_use, &beta, ZTdata, &token);

wcount=0;
for(j=0;j<bitlength;j++)
{
if(tindex[bitstart+j]==-9999&&mults[bitstart+j]!=-9999)	//not a top predictor, nor trivial
{
value=0;
for(k=0;k<num_fixed-1;k++)
{
value2=pow(ZTdata[k+j*(num_fixed-1)],2)/ZTZ[k]/sqdevs[bitstart+j]/num_samples_use;
if(value2>value){value=value2;}
}
if(value>mincor)
{
if(wcount<5){printf("Warning, Predictor %s has correlation squared %.2f with a fixed effect, so will be ignored\n", preds[bitstart+j], value);}
mults[bitstart+j]=-9999;
wcount++;
}
}}
if(wcount>5){printf("In total, %d predictors have correlation squared above %.4f with a fixed effect, so will be ignored\n", wcount, mincor);}
if(wcount>0){printf("\n");}
}

if(num_kins==1)	//get UTdata
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_samples_use, &bitlength, &num_samples_use, &alpha, U, &num_samples_use, data, &num_samples_use, &beta, UTdata, &num_samples_use);
}

if(mode==131&&num_kins==0&&num_fixed==1)	//get YTdata
{
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, Y, &one, &beta, YTdata, &one);
}

//ready to test
for(j=bitstart;j<bitend;j++)
{
if(mults[j]==-9999)	//trivial
{stats[0]=-9999;stats[1]=-9999;stats[2]=-9999;stats[3]=-9999;}
else	//will be testing
{
if(mode==131&&num_kins==0&&num_fixed==1)	//have simple case
{
value=YTdata[j-bitstart]/sqdevs[j]/num_samples_use;
value2=(varphen/sqdevs[j]-pow(value,2))/(num_samples_use-2);

stats[0]=value;
stats[1]=pow(value2,.5);
stats[2]=stats[0]/stats[1];
stats[3]=erfc(fabs(stats[0]/stats[1])*M_SQRT1_2);
}
else	//complex case - might be logistic, might have kinships, might have tops
{
if(tindex[j]!=-9999)	//a top predictor, so have already tested
{
stats[0]=thetas[num_covars+tindex[j]];
stats[1]=thetasds[num_covars+tindex[j]];
if(stats[1]!=-9999){stats[2]=stats[0]/stats[1];}
stats[3]=thetapvas[num_covars+tindex[j]];
}
else	//not a top, so must test
{
//add test predictor to end of Z
for(i=0;i<num_samples_use;i++){Z[i+num_fixed*num_samples_use]=data[(size_t)(j-bitstart)*num_samples_use+i];}

if(mode==131)
{
if(num_kins==0){reg_single_lin(Y, Z, num_samples_use, num_fixed+1, stats);}
else	//mixed model
{
//add UTdata for test predictor to end of UTZ, then test
for(i=0;i<num_samples_use;i++){UTZ[i+num_fixed*num_samples_use]=UTdata[(size_t)(j-bitstart)*num_samples_use+i];}
linear_reml(num_samples_use, num_fixed+1, Y, Z, U, E, UTY, UTZ, kintraces, stats, vstarts, NULL, NULL, NULL, constrain, tol, maxiter, exact);
}
}
else{reg_single_log(Y, Z, num_samples_use, num_fixed+1, stats, tol, maxiter);}
}
}	//end of complex case
}	//end of testing

//save results
if(stats[1]!=-9999)	//tested predictor - include in all results
{
//print assoc
fprintf(output2, "%d %s %.0f %c %c ", chr[j], preds[j], cmbp[j], al1[j], al2[j]);
fprintf(output2, "%.4f %.4e %.4e %.4e ", stats[2], stats[3], stats[0], stats[1]);
if(mode==131)
{
if(prev!=-9999){fprintf(output2, "%.4e %.4e ", stats[0]*factor, stats[1]*factor);}
else{fprintf(output2, "NA NA ");}
}
fprintf(output2, "%.6f ", centres[j]);
if(nonsnp==0){fprintf(output2, "%.6f\n", centres[j]/2+(centres[j]>1)*(1-centres[j]));}
else{fprintf(output2, "NA\n");}

//print summaries, pvalues and scores
fprintf(output3, "%s %c %c %d %.4f %d\n", preds[j], al1[j], al2[j], (stats[0]>=0)-(stats[0]<0), pow(stats[2],2), num_samples_use);
fprintf(output4, "%s %.4e\n", preds[j], stats[3]);
fprintf(output5, "%s %c %c %.6f %.6f", preds[j], al1[j], al2[j],  centres[j], stats[0]);
for(k=0;k<6;k++)
{
if(stats[3]<cuts[k]){fprintf(output5, " %.6f", stats[0]);}
else{fprintf(output5, " 0");}
}
fprintf(output5, "\n");
}
else	//trivial predictor - so include only in assoc
{
fprintf(output2, "%d %s %.0f %c %c ", chr[j], preds[j], cmbp[j], al1[j], al2[j]);
if(mode==131){fprintf(output2, "NA NA NA NA NA NA NA NA\n");}
else{fprintf(output2, "NA NA NA NA NA NA\n");}
}
}	//end of j loop
}	//end of bit loop
printf("\n");

fclose(output);
fclose(output2);
fclose(output3);
fclose(output4);
fclose(output5);

printf("Main results saved in %s, with a summary version in %s, p-values in %s and score file in %s\n\n", filename2, filename3, filename4, filename5);

free(data);
free(order);free(sweights);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);
free(tindex);free(stats);
if(num_fixed>1){free(ZTZ);free(ZTdata);}
if(num_kins==1){free(UTY);free(UTZ);free(UTdata);free(vstarts);}
if(mode==131&&num_kins==0&&num_fixed==1){free(YTdata);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

