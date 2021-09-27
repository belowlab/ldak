/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Construct ridge, bolt, bayesr and lasso models from summary statistics - assumes phenotype has variance one

///////////////////////////

//set neff, the assumed number of individuals
sum=0;for(j=0;j<data_length;j++){sum+=nss[j];}
neff=sum/data_length;

if(strcmp(sums2file,"blank")!=0)	//also using neff2
{
sum=0;for(j=0;j<data_length;j++){sum+=nss2[j];}
neff2=sum/data_length;
}

num_try=0;
if(ptype==1)	//lasso-sparse
{num_try+=80;}
if(ptype==0||ptype==2)	//lasso-shrink
{num_try+=1+10*multiher;}
if(ptype==0||ptype==3)	//ridge
{num_try+=1+10*multiher;}
if(ptype==0||ptype==4)
{num_try+=132-114*basicbolt-110*ldpred-(ptype==0);}
if(ptype==0||ptype==5)	//bayesr-sparse
{num_try+=83+(ptype==5);}
if(ptype==6)	//bayesr-shrink
{num_try+=84-(ptype==0);}

trytypes=malloc(sizeof(int)*num_try);
tryhers=malloc(sizeof(double)*num_try);
trylams=malloc(sizeof(double)*num_try);
tryscales=malloc(sizeof(double)*num_try);
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

//set tryscales to 1 and hers, lams, ps & f2s to NA, then will change when required
for(p=0;p<num_try;p++)
{tryscales[p]=1.0;tryhers[p]=-9999;trylams[p]=-9999;tryps[p]=-9999;tryp2s[p]=-9999;tryp3s[p]=-9999;tryp4s[p]=-9999;tryf2s[p]=-9999;}

count=0;

if(ptype==1)	//lasso-sparse
{
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.9;trylams[count]=0.001*pow(100,(float)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.5;trylams[count]=0.001*pow(100,(float)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.2;trylams[count]=0.001*pow(100,(float)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.1;trylams[count]=0.001*pow(100,(float)p/19);count++;}
}

if(ptype==0||ptype==2)	//lasso-shrink
{
if(multiher==1)
{
for(p=5;p<16;p++){trytypes[count]=2;tryhers[count]=p*0.1*her;count++;}
}
else{trytypes[count]=3;tryhers[count]=her;count++;}
}

if(ptype==0||ptype==3)	//ridge
{
if(multiher==1)
{
for(p=5;p<16;p++){trytypes[count]=3;tryhers[count]=p*0.1*her;count++;}
}
else{trytypes[count]=3;tryhers[count]=her;count++;}
}

if(ptype==0||ptype==4)	//bolt - skip first model if doing ridge
{
if(basicbolt==0&&ldpred==0)	//132 pairs (maybe minus ridge and f2=0)
{
loads=malloc(sizeof(double)*12);
loads2=malloc(sizeof(double)*11);
loads[0]=.5;loads[1]=.45;loads[2]=.4;loads[3]=.35;loads[4]=.3;loads[5]=.25;
loads[6]=.2;loads[7]=.15;loads[8]=.1;loads[9]=.05;loads[10]=0.02;loads[11]=0.01;
loads2[0]=.5;loads2[1]=.45;loads2[2]=.4;loads2[3]=.35;loads2[4]=.3;loads2[5]=.25;
loads2[6]=.2;loads2[7]=.15;loads2[8]=.1;loads2[9]=.05;loads2[10]=0;
for(j=0;j<132;j++)
{
if(j>0||ptype==4)
{trytypes[count]=4;tryhers[count]=her;tryps[count]=loads[j%12];tryp2s[count]=1-loads[j%12];tryf2s[count]=loads2[j/12];count++;}
}
free(loads);free(loads2);
}
else
{
if(basicbolt==1)	//just 18 pairs (maybe minus ridge)
{
loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
for(j=0;j<18;j++)
{
if(j>0||ptype==4)
{trytypes[count]=4;tryhers[count]=her;tryps[count]=loads[j%6];tryp2s[count]=1-loads[j%6];tryf2s[count]=loads2[j/6];count++;}
}
free(loads);free(loads2);
}
if(ldpred==1)	//just 22 pairs (maybe minus ridge)
{
loads=malloc(sizeof(double)*21);
loads[0]=.95;loads[1]=.9;loads[2]=.85;loads[3]=.8;loads[4]=.75;loads[5]=.7;loads[6]=.65;
loads[7]=.6;loads[8]=.55;loads[9]=.5;loads[10]=.45;loads[11]=.4;loads[12]=.35;loads[13]=.3;
loads[14]=.25;loads[15]=.2;loads[16]=.15;loads[17]=.1;loads[18]=.05;loads[19]=.02;loads[20]=.01;
if(ptype==4){trytypes[count]=4;tryhers[count]=her;tryps[count]=0.5;tryp2s[count]=0.5;tryf2s[count]=0.5;count++;}
for(j=0;j<21;j++)
{trytypes[count]=4;tryhers[count]=her;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=0;count++;}
}
free(loads);
}
}

if(ptype==0||ptype==5)	//bayesr with point mass (can not have all ps zero) - add in ridge model, unless also doing ridge
{
loads=malloc(sizeof(double)*7);
loads[0]=0;loads[1]=.001;loads[2]=.005;loads[3]=.01;loads[4]=.02;loads[5]=.05;loads[6]=.1;
for(j=0;j<7;j++){
for(j2=0;j2<7;j2++){
for(j3=0;j3<7;j3++){
if(j+j2+j3>0&&j2>=j&&j3>=j2)
{
trytypes[count]=5;tryhers[count]=her;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
if(ptype==5){trytypes[count]=5;tryhers[count]=her;tryps[count]=0;tryp2s[count]=0;tryp3s[count]=0;tryp4s[count]=1;count++;}
free(loads);
}

if(ptype==6)	//bayesr with four gaussians - skip all ps zero if doing ridge
{
loads=malloc(sizeof(double)*7);
loads[0]=0;loads[1]=.001;loads[2]=.005;loads[3]=.01;loads[4]=.02;loads[5]=.05;loads[6]=.1;
for(j=0;j<7;j++){
for(j2=0;j2<7;j2++){
for(j3=0;j3<7;j3++){
if(j2>=j&&j3>=j2&&(j+j2+j3>0||ptype==6))
{
trytypes[count]=6;tryhers[count]=her;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}

if(count!=num_try){printf("doug error %d and %d\n", num_try, count);exit(1);}

////////

//read headers from cors
usedpreds=malloc(sizeof(int)*num_preds);
maxnums=malloc(sizeof(int)*num_preds);
scumsums=malloc(sizeof(size_t)*num_preds);

sprintf(filename,"%s.cors.bin", corsfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read indicators and check compatible with keeppreds_use
fseeko(input, 0, SEEK_SET);
if(fread(usedpreds, sizeof(int), num_preds, input)!=num_preds)
{printf("Error reading predictor indicators from %s\n\n", filename);exit(1);}
for(j=0;j<data_length;j++)
{
if(usedpreds[keeppreds_use[j]]!=1){printf("Error, Predictor %s was not included when calculating predictor-predictor correlations; either remake the correlations or exclude predictors not included\n\n", preds[j]);exit(1);}
}

//read maxnums, and use to set cumsums
fseeko(input, (off_t)sizeof(int)*num_preds, SEEK_SET);
if(fread(maxnums, sizeof(int), num_preds, input)!=num_preds)
{printf("Error reading counts from %s\n\n", filename);exit(1);}
scumsums[0]=0;for(j=1;j<num_preds;j++){scumsums[j]=scumsums[j-1]+maxnums[j-1];}

fclose(input);

////////

//allocate variables
scount=0;for(j=0;j<data_length;j++){scount+=maxnums[keeppreds[j]];}
value=(double)scount/1024/1024/1024*8;
if(value>1){printf("Warning, to read the correlations requires %.1f Gb\n\n", value);}

bigs=malloc(sizeof(int*)*data_length);
rjks=malloc(sizeof(float*)*data_length);
for(j=0;j<data_length;j++)
{
bigs[j]=malloc(sizeof(int)*maxnums[keeppreds_use[j]]);
rjks[j]=malloc(sizeof(float)*maxnums[keeppreds_use[j]]);
}

actnums=malloc(sizeof(int)*data_length);
rjksums=malloc(sizeof(float)*data_length);

datasqs=malloc(sizeof(double)*data_length);
YTdata=malloc(sizeof(double)*data_length);

lambdas=malloc(sizeof(double)*data_length*num_try);
lambdas2=malloc(sizeof(double)*data_length*num_try);
lambdas3=malloc(sizeof(double)*data_length*num_try);
lambdas4=malloc(sizeof(double)*data_length*num_try);

effs=malloc(sizeof(double)*data_length*num_try);
effs2=malloc(sizeof(double)*data_length*num_try);
convs=malloc(sizeof(int)*data_length*num_try);

exps=malloc(sizeof(double)*data_length);

////////

//reopen cors
sprintf(filename,"%s.cors.bin", corsfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}

//read means, allowing for filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(centres+j, sizeof(double), 1, input)!=1)
{printf("Error reading mean for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}

//now scalings, allowing for filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(mults+j, sizeof(double), 1, input)!=1)
{printf("Error reading scaling for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}

//and variances, allowing for filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*2, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*2+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(sqdevs+j, sizeof(double), 1, input)!=1)
{printf("Error reading variance for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}

//finally read correlations - first with crude filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*3, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(maxnums[keeppreds_use[j]]>0)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*3+(sizeof(int)+sizeof(float))*scumsums[keeppreds_use[j]], SEEK_SET);}
if(fread(bigs[j], sizeof(int), maxnums[keeppreds_use[j]], input)!=maxnums[keeppreds_use[j]])
{printf("Error reading indexes for Predictor %d from %s\n\n", j+1, filename);exit(1);}
if(fread(rjks[j], sizeof(float), maxnums[keeppreds_use[j]], input)!=maxnums[keeppreds_use[j]])
{printf("Error reading correlations for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}
}
fclose(input);

//set actnums and rjksums, and replace rjks with strip * covjk (this also finishes any filtering)
for(j=0;j<num_preds;j++){usedpreds[j]=-1;}
for(j=0;j<data_length;j++){usedpreds[keeppreds_use[j]]=j;}

for(j=0;j<data_length;j++)
{
count=0;
sum=1;
for(j2=0;j2<maxnums[keeppreds_use[j]];j2++)
{
j3=bigs[j][j2];
if(usedpreds[j3]!=-1)	//predictor remains (and has been mapped to usedpreds[j3])
{
sum+=pow(rjks[j][j2],2);
bigs[j][count]=usedpreds[j3];
rjks[j][count]=shrink*rjks[j][j2]*pow(sqdevs[j],.5)*pow(sqdevs[usedpreds[j3]],.5);
count++;
}
}
actnums[j]=count;
rjksums[j]=sum;
}

////////

//get expected squared effect size of each predictor, divided by hers

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
set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, 1.0, tryhers[p], trylams, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, trytypes[p]);
}

//blank progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

////////

if(strcmp(sums2file,"blank")!=0)	//solve using training statistics - var(Y)=1, so YTY=neff2
{
//screen and file print
if(num_try==1){printf("Estimating effect sizes using secondary summary statistics\n\n");}
else{printf("Estimating effect sizes for %d models using secondary summary statistics (if using multiple cores, models will finish in a random order)\n\n", num_try);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Solving %d models using secondary summary statistics\n", num_try);
fprintf(output, "Model\tType\tHer\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tNum_Predictors_Failed\n");
fclose(output);

//get XTX and XTY for each predictor
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*neff2;YTdata[j]=rhos2[j]*pow(sqdevs[j],.5)*neff2;}

#pragma omp parallel for private(p,j,j2,count,sum,window_kb2,bitstart,bitstart2,bitend,sumsq,sumsq2) schedule(dynamic)
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[j+p*data_length]=0;}
else	//effect sizes begin at marginal postmeans under own model, weighted by tagging
{effs[j+p*data_length]=get_postmean(YTdata[j], lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p])/rjksums[j];}
}
else{effs[j+p*data_length]=0;}
}

//set convs to zero
for(j=0;j<data_length;j++){convs[j+p*data_length]=0;}

window_kb2=window_kb/segments;
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb2||chr[bitend]!=chr[bitstart]){break;}
}

//calculate ess for bit
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[j+p*data_length]-datasqs[j]*pow(effs[j+p*data_length],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend){sumsq-=rjks[j][j2]*neff2*effs[j+p*data_length]*effs[bigs[j][j2]+p*data_length];}
}
}

//save effect sizes for bit, in case fail to converge
for(j=bitstart;j<bitend;j++){effs2[j+p*data_length]=effs[j+p*data_length];}

count=0;
while(1)
{
count++;

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j];for(j2=0;j2<actnums[j];j2++){sum-=tryscales[p]*rjks[j][j2]*neff2*effs[bigs[j][j2]+p*data_length];}

//update effect size
if(trytypes[p]==1){effs[j+p*data_length]=get_postmean(sum, lambdas[j+p*data_length]*neff2, -9999, -9999, -9999, datasqs[j], 1.0, -9999, -9999, -9999, -9999, NULL, 1);}	
else{effs[j+p*data_length]=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p]);}
}
}	//end of j loop

//save old then calculate new ess for bit
sumsq2=sumsq;
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[j+p*data_length]-datasqs[j]*pow(effs[j+p*data_length],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend){sumsq-=rjks[j][j2]*neff2*effs[j+p*data_length]*effs[bigs[j][j2]+p*data_length];}
}
}

if(fabs(sumsq-sumsq2)/neff2<tol)	//converged
{
for(j=bitstart;j<bitend;j++){convs[j+p*data_length]=1;}

if(window_kb2==window_kb)	//did a full window, so move start forwards window_kb/segments
{
for(bitstart2=bitstart+1;bitstart2<data_length;bitstart2++)
{
if(cmbp[bitstart2]-cmbp[bitstart]>1000*window_kb/segments||chr[bitstart2]!=chr[bitstart]){break;}
}
if(bitstart2==bitend)	//must have been at end of chromosome
{window_kb2=window_kb/segments;}
bitstart=bitstart2;
}
else	//did a small window, so start stays same, but increase window wise
{window_kb2*=2;}

break;
}

if(count==maxiter)	//failed (almost certainly due to end of bit) - undo iterations, skip rest of bit, and reduce window size 
{
for(j=bitstart;j<bitend;j++){effs[j+p*data_length]=effs2[j+p*data_length];}
bitstart=bitend;
window_kb2=window_kb/segments;
break;
}
}	//end of inside while loop

}	//end of outside while loop

//count how many failed to converged
count=0;for(j=0;j<data_length;j++){count+=(exps[j]>0&&convs[j+p*data_length]==0);}

//print update
#pragma omp critical
{
printf("Solved Model %d: ", p+1);
if(trytypes[p]==1){printf("lasso-sparse, lambda %.4f, scale %.2f ", trylams[p], tryscales[p]);}
if(trytypes[p]==2){printf("lasso-shrink, her %.4f ", tryhers[p]);}
if(trytypes[p]==3){printf("ridge, her %.4f ", tryhers[p]);}
if(trytypes[p]==4){printf("bolt, her %.4f, p %.4f, f2 %.4f ", tryhers[p], tryps[p], tryf2s[p]);}
if(trytypes[p]==5){printf("BayesR, her %.4f, p %.4f, p2 %.4f, p3 %.4f ", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p]);}
if(trytypes[p]==6){printf("BayesR-shrink, her %.4f, p %.4f, p2 %.4f, p3 %.4f ", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p]);}
if(count==0){printf(" - effect sizes converged for all predictors\n");}
else{printf(" - effect sizes failed to converge for %d predictors\n", count);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "%d\t", p+1);
if(trytypes[p]==1){fprintf(output, "lasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\t%d\n", trylams[p], tryscales[p], count);}
if(trytypes[p]==2){fprintf(output, "lasso-shrink\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p], count);}
if(trytypes[p]==3){fprintf(output, "ridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p], count);}
if(trytypes[p]==4){fprintf(output, "bolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\t%d\n", tryhers[p], tryps[p], tryf2s[p], count);}
if(trytypes[p]==5){fprintf(output, "BayesR\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p], count);}
if(trytypes[p]==6){fprintf(output, "BayesR-shrink\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p], count);}
fclose(output);
}

}	//end of p loop

//save models
sprintf(filename2,"%s.effects.train",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(num_try==1)
{
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<num_tops;j++)
{fprintf(output2, "%s %c %c %.6f %.4e\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j]);}
for(j=0;j<data_length;j++)
{
if(effs[j]!=0){fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[j]);}
}
printf("\nModel saved in %s\n\n", filename2);
}
else
{
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
printf("\nModels saved in %s\n\n", filename2);
}

fclose(output2);
}	//end of using training samples

////////

//now solve for main statistics

//screen and file print
if(strcmp(sums2file,"blank")!=0)
{
if(num_try==1){printf("Estimating effect sizes using main summary statistics\n\n");}
else{printf("Estimating effect sizes for %d models using main summary statistics (if using multiple cores, models will finish in a random order)\n\n", num_try);}
}
else
{
if(num_try==1){printf("Estimating effect sizes\n");}
else{printf("Estimating effect sizes for %d models (if using multiple cores, models will finish in a random order)\n", num_try);}
}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(strcmp(sums2file,"blank")!=0)
{
if(num_try==1){fprintf(output,"Solving using main summary statistics\n");}
else{fprintf(output,"Solving %d models using main summary statistics\n", num_try);}
}
fprintf(output, "Model\tType\tHer\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tProp_Var_Exp\tLikelihood\tNum_Predictors_Failed\n");
fclose(output);

//get XTX and XTY for each predictor
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*neff;YTdata[j]=rhos[j]*pow(sqdevs[j],.5)*neff;}

#pragma omp parallel for private(p,j,j2,count,sum,window_kb2,bitstart,bitstart2,bitend,sumsq,sumsq2) schedule(dynamic)
for(p=0;p<num_try;p++)
{
//initialize values
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[j+p*data_length]=0;}
else	//effect sizes begin at marginal postmeans under own model, weighted by tagging
{effs[j+p*data_length]=get_postmean(YTdata[j], lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p])/rjksums[j];}	
}
else{effs[j+p*data_length]=0;}
}

//set convs to zero
for(j=0;j<data_length;j++){convs[j+p*data_length]=0;}

window_kb2=window_kb/segments;
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb2||chr[bitend]!=chr[bitstart]){break;}
}

//calculate ess for bit
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[j+p*data_length]-datasqs[j]*pow(effs[j+p*data_length],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend){sumsq-=rjks[j][j2]*neff*effs[j+p*data_length]*effs[bigs[j][j2]+p*data_length];}
}
}

//save effect sizes for bit, in case fail to converge
for(j=bitstart;j<bitend;j++){effs2[j+p*data_length]=effs[j+p*data_length];}

count=0;
while(1)
{
count++;

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j];for(j2=0;j2<actnums[j];j2++){sum-=tryscales[p]*rjks[j][j2]*neff*effs[bigs[j][j2]+p*data_length];}


//update effect size
if(trytypes[p]==1){effs[j+p*data_length]=get_postmean(sum, lambdas[j+p*data_length]*neff, -9999, -9999, -9999, datasqs[j], 1.0, -9999, -9999, -9999, -9999, NULL, 1);}
else{effs[j+p*data_length]=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p]);}
}
}	//end of j loop

//save old then calculate new ess for bit
sumsq2=sumsq;
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[j+p*data_length]-datasqs[j]*pow(effs[j+p*data_length],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend){sumsq-=rjks[j][j2]*neff*effs[j+p*data_length]*effs[bigs[j][j2]+p*data_length];}
}
}

if(fabs(sumsq-sumsq2)/neff<tol)	//converged
{
for(j=bitstart;j<bitend;j++){convs[j+p*data_length]=1;}

if(window_kb2==window_kb)	//did a full window, so move start forwards window_kb/segments
{
for(bitstart2=bitstart+1;bitstart2<data_length;bitstart2++)
{
if(cmbp[bitstart2]-cmbp[bitstart]>1000*window_kb/segments||chr[bitstart2]!=chr[bitstart]){break;}
}
if(bitstart2==bitend)	//must have been at end of chromosome
{window_kb2=window_kb/segments;}
bitstart=bitstart2;
}
else	//did a small window, so start stays same, but increase window size
{window_kb2*=2;}

break;
}

if(count==maxiter)	//failed (almost certainly due to end of bit) - undo iterations, skip rest of bit, and reduce window size 
{
for(j=bitstart;j<bitend;j++){effs[j+p*data_length]=effs2[j+p*data_length];}
bitstart=bitend;
window_kb2=window_kb/segments;
break;
}
}	//end of inside while loop

}	//end of outside while loop

//count how many failed to converged
count=0;for(j=0;j<data_length;j++){count+=(exps[j]>0&&convs[j+p*data_length]==0);}

//print update
#pragma omp critical
{
if(num_try==1)
{
if(count==0){printf("Effect sizes converged for all predictors\n");}
else{printf("Effect sizes failed to converge for %d predictors\n", count);}
}
else
{
printf("Solved Model %d: ", p+1);
if(trytypes[p]==1){printf("lasso-sparse, lambda %.4f, scale %.4f ", trylams[p], tryscales[p]);}
if(trytypes[p]==2){printf("lasso-shrink, her %.4f ", tryhers[p]);}
if(trytypes[p]==3){printf("ridge, her %.4f ", tryhers[p]);}
if(trytypes[p]==4){printf("bolt, her %.4f, p %.4f, f2 %.4f ", tryhers[p], tryps[p], tryf2s[p]);}
if(trytypes[p]==5){printf("BayesR, her %.4f, p %.4f, p2 %.4f, p3 %.4f ", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p]);}
if(trytypes[p]==6){printf("BayesR-shrink, her %.4f, p %.4f, p2 %.4f, p3 %.4f ", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p]);}
if(count==0){printf("- effect sizes converged for all predictors\n");}
else{printf("- effect sizes failed to converge for %d predictors\n", count);}
}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "%d\t", p+1);
if(trytypes[p]==1){fprintf(output, "lasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\t%d\n", trylams[p], tryscales[p], count);}
if(trytypes[p]==2){fprintf(output, "lasso-shrink\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p], count);}
if(trytypes[p]==3){fprintf(output, "ridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p], count);}
if(trytypes[p]==4){fprintf(output, "bolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\t%d\n", tryhers[p], tryps[p], tryf2s[p], count);}
if(trytypes[p]==5){fprintf(output, "BayesR\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p], count);}
if(trytypes[p]==6){fprintf(output, "BayesR-shrink\t%.4f\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p], tryp4s[p], tryp3s[p], tryp2s[p], count);}
fclose(output);
}
}	//end of p loop

//save models
if(strcmp(sums2file,"blank")!=0){sprintf(filename2,"%s.effects.final",outfile);}
else{sprintf(filename2,"%s.effects",outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(num_try==1)
{
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<num_tops;j++)
{fprintf(output2, "%s %c %c %.6f %.4e\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j]);}
for(j=0;j<data_length;j++)
{
if(effs[j]!=0){fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[j]);}
}
printf("\nModel saved in %s\n\n", filename2);
}
else
{
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
printf("\nModels saved in %s\n\n", filename2);
}

fclose(output2);

if(strcmp(sums2file,"blank")!=0)
{printf("You should now decide the best single model; run \"--calc-scores\" with options \"--scorefile %s.effects.train\" and \"--final-effects %s.effects.final\"\n\n", outfile, outfile);}

////////

free(trytypes);free(tryhers);free(trylams);free(tryscales);free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);free(tryf2s);
free(usedpreds);free(maxnums);free(scumsums);
for(j=0;j<data_length;j++){free(bigs[j]);free(rjks[j]);}free(bigs);free(rjks);
free(actnums);free(rjksums);
free(datasqs);free(YTdata);
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(effs);free(effs2);free(convs);
free(exps);

///////////////////////////
