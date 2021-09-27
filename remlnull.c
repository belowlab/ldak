/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You shouldt have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Solves NULL model prior to complex mixed model regression, then makes combined matrix and decomposes

///////////////////////////

void null_reml(int ns, int num_fixed, int num_kins, int num_regs, double *Y, double *Z, double **mkins, double *kintraces, double *kinsums, double *X, int Xtotal, int *Xstarts, int *Xends, double *Xsums, double prev, int np, int constrain, double tol, int maxiter, int memsave, char **kinstems, char **ids3, char *outfile, char **ids1, char **ids2)
{
size_t scount, stotal;
int i, j, j2, k, k2, r, count, token, flag, one=1, info;
double value, value2, sum, sum2, sumsq, max, alpha, beta;

int num_fixed, total, cflag, rflag, *fixed, *fixedsave, nlost, nlost2;
double nfree, varnull, varmin, gam, gam2, gam3, relax;
double likenull, like, likeold, diff, lrtstat, lrtpva, covher, topher, factor;
double *scales, *vars, *vardiffs, *varsds, *hers, *hersds, *shares, *sharesds;

double *AI, *AI2, *AI3, *AIsave, *BI, *BIsave, *J, *JAI, *JAIJT;
double *ZTY, *ZTZ, *ZTZ2, *ZTZZTY, detZTZ;
double detV, *ZTVZ, *ZTVZ2, *ZTVZ3, detZTVZ, *PY, **KPY, **PKPY, *traces;

double *V, *V2, *V3, *VZ, *VZZTVZ, *P;
double *PX, *XTPY;

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500], filename10[500];
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7, *output8, *output9, *output10;

//variables to make compatible with remllikes.c and remlderiv.c
int shortcut=0;
double *UTY, *UTZ, *D, detD, *BUTZ, *H, *HUTY, *HKPY;
double *UTX, *DUTX, detC, *XTVX, *XTVX2, *XTVX3, detXTVX, *F, *FUTZ, *FUTY, *FUTX, *HUTX, *FKPY;


//set nfree, total and stotal
nfree=ns-num_fixed;
total=1+num_kins+num_regs;
stotal=(size_t)ns*ns;

//allocate variables

fixed=malloc(sizeof(int)*total);
fixedsave=malloc(sizeof(int)*total);
vars=malloc(sizeof(double)*total);
vardiffs=malloc(sizeof(double)*total);
varsds=malloc(sizeof(double)*total);
scales=malloc(sizeof(double)*total);
hers=malloc(sizeof(double)*total);
hersds=malloc(sizeof(double)*total);
shares=malloc(sizeof(double)*total);
sharesds=malloc(sizeof(double)*total);

AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
AIsave=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);
BIsave=malloc(sizeof(double)*total);
J=malloc(sizeof(double)*total*total);
JAI=malloc(sizeof(double)*total*total);
JAIJT=malloc(sizeof(double)*total*total);

ZTY=malloc(sizeof(double)*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZZTY=malloc(sizeof(double)*num_fixed);

ZTVZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTVZ2=malloc(sizeof(double)*num_fixed);
ZTVZ3=malloc(sizeof(double)*num_fixed*num_fixed);
PY=malloc(sizeof(double)*ns);
KPY=malloc(sizeof(double*)*total);
for(k=0;k<total;k++){KPY[k]=malloc(sizeof(double)*ns);}
PKPY=malloc(sizeof(double*)*total);
for(k=0;k<total;k++){PKPY[k]=malloc(sizeof(double)*ns);}
traces=malloc(sizeof(double)*total);

V=malloc(sizeof(double)*ns*ns);
V2=malloc(sizeof(double)*ns);
VZ=malloc(sizeof(double)*ns*num_fixed);
VZZTVZ=malloc(sizeof(double)*ns*num_fixed);
P=malloc(sizeof(double)*ns*ns);

if(Xtotal>0)
{
PX=malloc(sizeof(double)*ns*Xtotal);
XTPY=malloc(sizeof(double)*Xtotal);
}

////////

//solve model with just fixed effects to get varnull, varmin and null likelihood
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
for(j=0;j<num_fixed;j++){ZTZZTY[j]=ZTY[j];}
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, 1, ZTZZTY, 1);

sumsq=0;
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
for(j=0;j<num_fixed;j++){sumsq-=ZTY[j]*ZTZZTY[j];}
varnull=sumsq/nfree;
varmin=0.0001*varnull;
likenull=-.5*(nfree+nfree*log(2*M_PI*varnull)+detZTZ);

//set scales
scales[0]=1;
for(k=0;k<num_kins;k++){scales[1+k]=kintraces[k];}
for(r=0;r<num_regs;r++){scales[1+num_kins+r]=1;}

//set starting vars, hers, shares and fixed
hers[0]=.5;shares[0]=1;
for(k=1;k<total;k++){hers[k]=.5/(total-1);shares[k]=1.0/(total-1);}
for(k=0;k<total;k++){vars[k]=hers[k]*varnull/scales[k];}
for(k=0;k<total;k++){fixed[k]=0;}

//prepare to screen and file print
printf("Iter\t");
for(k=0;k<num_kins;k++){printf("Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){printf("Her_R%d\t", r+1);}
printf("Her_All\tLikelihood\tDifference\tTarget\tNum_Constrained\n");

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Iteration\t");
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d\t", r+1);}
fprintf(output, "Her_All\tTotal_Variance\tLikelihood\tDifference\tTarget\tNum_Constrained\n");
fclose(output);

////////

//now iterate
count=0;
rflag=0;	//0 if normal moves, 1 if reduced moves, 2 if transitioning from reduced to normal moves
cflag=1;
while(1)
{
//compute invV, detV, invZTVZ, detZTVZ, P, PY and gam (plus other variables depending on shortcut)
#include "remllike.c"

//get likelihood
like=-.5*gam-.5*detZTVZ-.5*detV-.5*nfree*log(2*M_PI);

if(count>0)	//set diff and decide what type of move to do
{
if(like>likeold-tol)	//move was fine, so next move will be normal or transitioning
{
diff=like-likeold;
if(rflag==1){rflag=2;}
else{rflag=0;}
}
else	//move was poor, so return to previous state and next move will be reduced
{
printf("Warning, the last move reduced the likelihood, so have returned to the previous state\n");
for(k=0;k<total;k++){fixed[k]=fixedsave[k];}
for(k=0;k<total;k++){vars[k]-=vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}
like=likeold;
diff=0;
rflag=1;
}
}
likeold=like;
for(k=0;k<total;k++){fixedsave[k]=fixed[k];}


if(rflag==0||rflag==2)	//compute derivatives
{
//compute PX, XTPY, KPY, PKPY, (inverse) AI and BI (plus other variables, depending on shortcut)
#include "remlderiv.c"

//save derivatives in case returning
for(k=0;k<total;k++)
{
BIsave[k]=BI[k];
for(k2=0;k2<total;k2++){AIsave[k+k2*total]=AI[k+k2*total];}
}
}
else	//recover saved values
{
for(k=0;k<total;k++)
{
BI[k]=BIsave[k];
for(k2=0;k2<total;k2++){AI[k+k2*total]=AIsave[k+k2*total];}
}
}

//print update
sum=0;for(k=0;k<num_kins+num_regs;k++){sum+=hers[1+k];}
value=0;for(k=0;k<total;k++){value+=vars[k];}
nlost=0;for(k=0;k<total;k++){nlost+=(fixed[k]>=3);}
nlost2=0;for(k=0;k<total;k++){nlost2+=(fixed[k]==1||fixed[k]==2);}

if(count==0){printf("Start\t");}
else{printf("%d\t", count);}
for(k=0;k<num_kins+num_regs;k++){printf("%.4f\t", hers[1+k]);}
printf("%.4f\t", sum);
if(rflag==0||rflag==2){printf("%.2f\t", like);}
else{printf("%.2f*\t", like);}
if(count==0){printf("n/a\t\t%.6f\t", tol);}
else{printf("%.6f\t%.6f\t", diff, tol);}
if(nlost2==0){printf("%d\n", nlost);}
else{printf("%d*\n", nlost);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "%d\t", count);
for(k=0;k<num_kins+num_regs;k++){fprintf(output, "%.6f\t", hers[1+k]);}
if(count==0){fprintf(output, "%.6f\t%.6f\t%.6f\tNA\t%.6f\t%d\n", sum, value, like, tol, nlost);}
else{fprintf(output, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n", sum, value, like, diff, tol, nlost);}
fclose(output);

//see if breaking (normally can only break if rflag=0 and nlost2=0, unless at iter limit)

if(num_kins+num_regs==0){break;}	//null model
if(nlost>=num_kins+num_regs){printf("All heritabilities are constrained\n");break;}
if(count>0)
{
if(fabs(diff)<tol&&rflag==0&&nlost2==0){break;}
}
if(count==maxiter)
{
printf("\nWarning, REML failed to converge within %d iterations; consider using \"--max-iter\" and/or \"--tolerance\" to change the iteration limit and tolerance", maxiter);
if(constrain==0)
{printf("; additionally, it might help to use \"--constrain YES\" to constrain (region) heritabilities to [0,1]\n");}
else
{printf("; additionally, it might help to use \"--constrain NO\" to allow heritabilities outside [0,1]\n");}
cflag=0;
break;
}

////////

//update variances using NR

//decide how far to move
if(rflag==0||rflag==2){relax=1;}
else{relax*=.5;}

//get proposed moves, ensuring not too large
alpha=relax;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, vardiffs, &one);

max=0;
for(k=0;k<total;k++)
{
if(fabs(vardiffs[k])>max){max=fabs(vardiffs[k]);}
}
if(max>.1*varnull)	//then reduce moves
{
relax*=.1*varnull/max;
for(k=0;k<total;k++){vardiffs[k]*=.1*varnull/max;}
}

if(constrain==1)	//variances can not be negative
{
for(k=0;k<total;k++)
{
if(fixed[k]<3)	//free to update
{
if(vars[k]+vardiffs[k]<varmin){vardiffs[k]=varmin-vars[k];fixed[k]++;}
else{fixed[k]=0;}
if(fixed[k]==3){vardiffs[k]=-vars[k];}
}
}}

//now move
for(k=0;k<total;k++){vars[k]+=vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}

count++;
}	//end of while loop
printf("\n");

//get some stats
lrtstat=2*(like-likenull);
lrtpva=.5*erfc(pow(lrtstat,.5)*M_SQRT1_2);
if(lrtstat<0){lrtpva=.75;}
if(hers[0]<=0){lrtpva=1;}

////////

//get SDs - for transformed variances, must compute J invAI JT, where Jij=dnewi/dvarsj
//values of AI corresponding to fixed components will have been set to zero; do same for J

//load up SDs of vars direct from AI
for(k=0;k<total;k++)
{
if(AI[k+k*total]>=0){varsds[k]=pow(AI[k+k*total],.5);}
else{varsds[k]=-9999;}
}

//get SDs of hers - Jij=delta*scalei/sum-scalei*scalej*vari/sum^2 (where sum across all elements)
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}

for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
if(fixed[k]<3&&fixed[k2]<3){J[k+k2*total]=-scales[k]*vars[k]*scales[k2]*pow(sum,-2);}
else{J[k+k2*total]=0;}
}
if(fixed[k]<3){J[k+k*total]+=scales[k]/sum;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

//load up
for(k=1;k<total;k++)
{
if(JAIJT[k+k*total]>=0){hersds[k]=pow(JAIJT[k+k*total],.5);}
else{hersds[k]=-9999;}
}

//save details for her_all as first element of hers and hersds
hers[0]=0;for(k=1;k<total;k++){hers[0]+=hers[k];}
sum=0;
for(k=1;k<total;k++)
{
for(k2=1;k2<total;k2++){sum+=JAIJT[k+k2*total];}
}
if(sum>=0){hersds[0]=pow(sum,.5);}
else{hersds[0]=-9999;}


//get SDs of shares - Jij=delta*scalei/sum-scalei*scalej*vari/sum^2 (where sum excludes noise term)
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}

//first row and column (easier to simply set new0=var0)
if(fixed[0]<3){J[0]=1;}
else{J[0]=0;}
for(k=1;k<total;k++){J[k]=0;J[k*total]=0;}

//rest
for(k=1;k<total;k++)
{
for(k2=1;k2<total;k2++)
{
if(fixed[k]<3&&fixed[k2]<3){J[k+k2*total]=-scales[k]*vars[k]*scales[k2]*pow(sum,-2);}
else{J[k+k2*total]=0;}
}
if(fixed[k]<3){J[k+k*total]+=scales[k]/sum;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

//load up
for(k=0;k<total;k++)
{
if(JAIJT[k+k*total]>=0){sharesds[k]=pow(JAIJT[k+k*total],.5);}
else{sharesds[k]=-9999;}
}

///////////////////////////

//get variance due to covariates and top predictors - will already have ZTY

sumsq=-ZTY[0]/ns*ZTY[0];
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
value=-ZTY[0]/ns*ZTY[0];for(j=0;j<num_covars;j++){value+=ZTY[j]*thetas[j];}
value2=0;for(j=num_covars;j<num_fixed;j++){value2+=ZTY[j]*thetas[j];}
covher=value/sumsq;
topher=value2/(sumsq-value);

if(num_covars>1){printf("Proportion of variance explained by the %d covariates: %.4f\n", num_covars, covher);}
if(num_tops==1){printf("Proportion of variance explained by the top predictor: %.4f\n", topher);}
if(num_tops>1){printf("Proportion of variance explained by the %d top predictors: %.4f\n", num_tops, topher);}
if(num_covars>1||num_tops>1){printf("\n");}

//adjust for tops
for(k=0;k<total;k++)
{
hers[k]*=(1-topher);
if(hersds[k]!=-9999){hersds[k]*=(1-topher);}
}

//save .reml, .coeff, .share, .vars, indi.res, .indi.blp, .reg.blup, .reg.score and .cross, then liability versions

flag=0;
for(k=0;k<num_kins;k++){flag+=(kinsums[k]==-9999);}
if(flag==0)
{
sum=0;
for(k=0;k<num_kins;k++){sum+=kinsums[k];}
for(r=0;r<num_regs;r++){sum+=Xsums[r];}
}

sprintf(filename2,"%s.reml", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
if(num_kins+num_regs>0){fprintf(output2, "Blupfile %s.indi.blp\n", outfile);}
else{fprintf(output2, "Blupfile none\n");}
if(num_regs>0||num_tops>0){fprintf(output2, "Regfile %s.reg.blup\n", outfile);}
else{fprintf(output2, "Regfile none\n");}
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
if(cflag==1){fprintf(output2,"Converged YES\n");}
else{fprintf(output2,"Converged NO\n");}
fprintf(output2, "Null_Likelihood %.4f\nAlt_Likelihood %.4f\n", likenull, like);
if(num_kins+num_regs==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability Her_SD Size Mega_Intensity Int_SD\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[1+k], hersds[1+k], kinsums[k], hers[1+k]/kinsums[k]*1000000, hersds[1+k]/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[1+num_kins+r], hersds[1+num_kins+r], Xsums[r], hers[1+num_kins+r]/Xsums[r]*1000000, hersds[1+num_kins+r]/Xsums[r]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
if(num_kins+num_regs==0){fprintf(output2, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[0]+topher, hersds[0], sum, hers[0]/sum*1000000, hersds[0]/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[1+k], hersds[1+k]);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[1+num_kins+r], hersds[1+num_kins+r]);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f NA NA NA\n", hers[0]+topher, hersds[0]);
}
fclose(output2);

if(compact==0)
{
sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Component Effect SD P\n");
fprintf(output3, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);} 
for(j=0;j<num_envs;j++){fprintf(output3, "Enviromental_%d %.6f %.6f %.4e\n",j, thetas[num_covars+j], thetasds[num_covars+j], thetapvas[num_covars+j]);}
fclose(output3);

sprintf(filename4,"%s.share", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SD Expected Enrichment SD\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, shares[1+k], sharesds[1+k], kinsums[k]/sum, shares[1+k]/kinsums[k]*sum, sharesds[1+k]/kinsums[k]*sum);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f %.6f %.6f %.6f\n", r+1, shares[1+num_kins+r], sharesds[1+num_kins+r], Xsums[r]/sum, shares[1+num_kins+r]/Xsums[r]*sum, sharesds[1+num_kins+r]/Xsums[r]*sum);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f NA NA NA\n", k+1, shares[1+k], sharesds[1+k]);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f NA NA NA\n", r+1, shares[1+num_kins+r], sharesds[1+num_kins+r]);}
}
fclose(output4);

sprintf(filename5,"%s.vars", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Component Variance SD\n");
for(k=0;k<num_kins;k++){fprintf(output5, "Var_K%d %.6f %.6f\n", k+1, vars[1+k], varsds[1+k]);}
for(r=0;r<num_regs;r++){fprintf(output5, "Var_R%d %.6f %.6f\n", r+1, vars[1+num_kins+r], varsds[1+num_kins+r]);}
fprintf(output5, "Var_E %.6f %.6f\n", vars[0], varsds[0]);
fclose(output5);

sprintf(filename6,"%s.indi.res", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "ID1\tID2\tPhenotype\tFitted\tResidual\n");
for(i=0;i<ns;i++){fprintf(output6, "%s\t%s\t%f\t%f\t%f\n", ids1[i], ids2[i], Y[i], Yadj[i], Y[i]-Yadj[i]);}
fclose(output6);

if(num_kins+num_regs>0)
{
sprintf(filename7,"%s.indi.blp", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(i=0;i<ns;i++)
{
fprintf(output7, "%s\t%s\t", ids1[i], ids2[i]);
for(k=0;k<num_kins;k++){fprintf(output7, "%.6f\t%.6f\t", mg2[k][i], mg[k][i]);}
for(r=0;r<num_regs;r++){fprintf(output7, "0\t%.6f\t", mg[num_kins+r][i]);}
fprintf(output7, "\n");
}
fclose(output7);
}

if(num_regs>0||num_tops>0)
{
sprintf(filename8,"%s.reg.blup", outfile);
if((output8=fopen(filename8,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename8);exit(1);}
fprintf(output8, "Predictor\tA1\tA2\tCentre\t");
for(r=0;r<num_regs;r++){fprintf(output8, "Region%d\t", r+1);}
if(num_tops>0){fprintf(output8, "Top_Preds\n");}
else{fprintf(output8, "\n");}
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{
fprintf(output8,"%s\t%c\t%c\t%.6f\t", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]);
for(r=0;r<num_regs;r++){fprintf(output8, "%.6f\t", effects[r][j]);}
if(num_tops>0){fprintf(output8, "%.6f\n", effects[num_regs][j]);}
else{fprintf(output8,"\n");}
}
}
fclose(output8);

sprintf(filename9,"%s.reg.score", outfile);
if((output9=fopen(filename9,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename9);exit(1);}
fprintf(output9, "Predictor\tA1\tA2\tCentre\tEffect\n");
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output9,"%s\t%c\t%c\t%.6f\t%.6f\n", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j], effects[num_regs+1][j]);}
}
fclose(output9);
}

sprintf(filename10,"%s.cross", outfile);
if((output10=fopen(filename10,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename10);exit(1);}
for(k=0;k<num_kins;k++){fprintf(output10, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output10, "Her_R%d\t", r+1);}
fprintf(output10, "\n");
for(k=1;k<total;k++)
{
for(k2=1;k2<total;k2++)
{fprintf(output10, "%.6f\t", JAIJTsave[k+k2*total]);}
fprintf(output10, "\n");
}
fclose(output10);
}	//end of compact=0

////////

if(prev!=-9999)
{
factor=get_factor(Y, ns, prev, -9999, outfile);

sprintf(filename,"%s.liab", filename2);	//.reml
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
if(num_kins+num_regs>0){fprintf(output, "Blupfile %s.indi.blp.liab\n", outfile);}
else{fprintf(output, "Blupfile none\n");}
if(num_regs>0||num_tops>0){fprintf(output, "Regfile %s.reg.blup.liab\n", outfile);}
else{fprintf(output, "Regfile none\n");}
fprintf(output, "Coeffsfile %s.coeff.liab\nCovar_Heritability %.4f\n", outfile, covher*factor);
fprintf(output, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
if(cflag==1){fprintf(output2,"Converged YES\n");}
else{fprintf(output2,"Converged NO\n");}
fprintf(output, "Null_Likelihood %.4f\nAlt_Likelihood %.4f\n", likenull, like);
if(num_kins+num_regs==1){fprintf(output, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output, "Component Heritability Her_SD Size Mega_Intensity Int_SD\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[1+k]*factor, hersds[1+k]*factor, kinsums[k], hers[1+k]/kinsums[k]*1000000*factor, hersds[1+k]/kinsums[k]*1000000*factor);}
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[1+num_kins+r]*factor, hersds[1+num_kins+r]*factor, Xsums[r], hers[1+num_kins+r]/Xsums[r]*1000000*factor, hersds[1+num_kins+r]/Xsums[r]*1000000*factor);}
fprintf(output, "Her_Top %.6f NA NA NA NA\n", topher*factor);
if(num_kins+num_regs==0){fprintf(output, "Her_All %.6f NA NA NA NA\n", topher*factor);}
else{fprintf(output, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[0]*factor+topher*factor, hersds[0]*factor, sum, hers[0]/sum*1000000*factor, hersds[0]/sum*1000000*factor);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[1+k]*factor, hersds[1+k]*factor);}
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[1+num_kins+r]*factor, hersds[1+num_kins+r]*factor);}
fprintf(output, "Her_Top %.6f NA NA NA NA\n", topher*factor);
fprintf(output, "Her_All %.6f %.6f NA NA NA\n", hers[0]*factor+topher*factor, hersds[0]*factor);
}
fclose(output);

if(compact==0)
{
sprintf(filename,"%s.liab", filename3);	//.coeff
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Effect SD P\n");
fprintf(output, "Intercept %.6f %.6f %.4e\n", thetas[0]*pow(factor,.5), thetasds[0]*pow(factor,.5), thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j]*pow(factor,.5), thetasds[j]*pow(factor,.5), thetapvas[j]);} 
for(j=0;j<num_envs;j++){fprintf(output, "Enviromental_%d %.6f %.6f %.4e\n",j, thetas[num_covars+j]*pow(factor,.5), thetasds[num_covars+j]*pow(factor,.5), thetapvas[num_covars+j]);}
fclose(output);

sprintf(filename,"%s.liab", filename7);	//.indi.blp
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(i=0;i<ns;i++)
{
fprintf(output, "%s\t%s\t", ids1[i], ids2[i]);
for(k=0;k<num_kins;k++){fprintf(output, "%.6f\t%.6f\t", mg2[k][i]*pow(factor,.5), mg[k][i]*pow(factor,.5));}
for(r=0;r<num_regs;r++){fprintf(output, "0\t%.6f\t", mg[num_kins+r][i]*pow(factor,.5));}
fprintf(output, "\n");
}
fclose(output);

if(num_regs>0||num_tops>0)
{
sprintf(filename,"%s.liab", filename8);	//.reg.blup
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\t");
for(r=0;r<num_regs;r++){fprintf(output, "Region%d\t", r+1);}
if(num_tops>0){fprintf(output, "Top_Preds\n");}
else{fprintf(output, "\n");}
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output,"%s\t%c\t%c\t%.6f\t", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]);}
for(r=0;r<num_regs;r++){fprintf(output, "%.6f\t", effects[r][j]*pow(factor,.5));}
if(num_tops>0){fprintf(output, "%.6f\n", effects[num_regs][j]*pow(factor,.5));}
else{fprintf(output,"\n");}
}
fclose(output);

sprintf(filename,"%s.liab", filename9);	//.reg.score
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\tEffect\n");
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output9,"%s\t%c\t%c\t%.6f\t%.6f\n", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]*pow(factor,.5), effects[num_regs+1][j]*pow(factor,.5));}
}
fclose(output);
}
}	//end of compact=0
}	//end of binary

if(compact==1){printf("Results saved in %s", filename2);}
else{printf("Main results saved in %s", filename2);}
if(prev!=-9999){printf(", with a liability version saved in %s.liab", filename2);}
printf("\n\n");

////////

if(discenv==1)	//there are 1+num_envs+num_regs genetic kinships - get proportion of variance these explain
//there is no need to compute the entire jacobian, as we only need the top line (gen sums / total sums)
{
printf("Computing heritabilities for %d subgroups\n", num_envs);

//first do for all individuals - can use existing scales
sum=0;
for(k=0;k<1+num_envs;k++){sum+=scales[1+k]*vars[1+k];}
for(r=0;r<num_regs;r++){sum+=scales[1+num_kins+r]*vars[1+num_kins+r];}
sum2=0;for(k=0;k<total;k++){sum2+=scales[k]*vars[k];}

//new1 = sum/sum2 (equals sum scale*var for gens / sum scale*var for all)
for(k=0;k<total;k++){SJ[k]=-scales[k]*sum*pow(sum2,-2);}
for(k=0;k<1+num_envs;k++){SJ[1+k]+=scales[1+k]/sum2;}
for(r=0;r<num_regs;r++){SJ[1+num_kins+r]+=scales[1+num_kins+r]/sum2;}

Snums[num_envs]=ns;
Shers[num_envs]=sum/sum2;
value=0;for(k=0;k<total;k++){for(k2=0;k2<total;k2++){value+=AI[k+k2*total]*SJ[k]*SJ[k2];}}
if(value>=0){Shersds[num_envs]=pow(value,.5);}
else{Shersds[num_envs]=-9999;}

//now for each subgroup
for(j=0;j<num_envs;j++)
{
//work out individuals in subgroup
count=0;
for(i=0;i<ns;i++)
{
if(Z[i+(num_covars+num_tops+j)*ns]==1){Sindexer[count]=i;count++;}
}

if(memsave==1)	//load these into wantids
{
for(i=0;i<count;i++){wantids[i]=ids3[Sindexer[i]];}
}

//get traces across these individuals
Sscales[0]=1;
for(k=0;k<num_kins;k++)
{
if(memsave==0)
{
sum=0;for(i=0;i<count;i++){sum+=mkins[k][(size_t)Sindexer[i]*ns+Sindexer[i]];}
Sscales[1+k]=sum/count;
}
else{Sscales[1+k]=read_kin_trace(kinstems[k], count, wantids, 1);}
}
for(r=0;r<num_regs;r++)
{
sum=0;
for(i=0;i<count;i++)
{
for(k=Xstarts[r];k<Xends[r];k++){sum+=pow(X[Sindexer[i]+k*ns],2);}
}
Sscales[1+num_kins+r]=sum/Xsums[r]/count;
}

sum=0;
for(k=0;k<1+num_envs;k++){sum+=Sscales[1+k]*vars[1+k];}
for(r=0;r<num_regs;r++){sum+=Sscales[1+num_kins+r]*vars[1+num_kins+r];}
sum2=0;for(k=0;k<total;k++){sum2+=Sscales[k]*vars[k];}

//new1 = sum/sum2 (equals sum scale*var for gens / sum scale*var for all)
for(k=0;k<total;k++){SJ[k]=-scales[k]*sum*pow(sum2,-2);}
for(k=0;k<1+num_envs;k++){SJ[1+k]+=scales[1+k]/sum2;}
for(r=0;r<num_regs;r++){SJ[1+num_kins+r]+=scales[1+num_kins+r]/sum2;}

Snums[j]=count;
Shers[j]=sum/sum2;
value=0;for(k=0;k<total;k++){for(k2=0;k2<total;k2++){value+=AI[k+k2*total]*SJ[k]*SJ[k2];}}
Shersds[j]=pow(value,.5);
if(value<0){Shersds[j]=-9999;}
}	//end of j loop

sprintf(filename,"%s.subgroups", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Num_Samples Heritability Her_SD\n");
for(j=0;j<num_envs;j++){fprintf(output, "Her_Sub%d %d %.6f %.6f\n", j+1, Snums[j], Shers[j], Shersds[j]);}
fprintf(output, "Her_All %d %.6f %.6f\n", Snums[num_envs], Shers[num_envs], Shersds[num_envs]);
fclose(output);

printf("Estimates saved in %s\n\n", filename);
free(Sindexer);free(Snums);free(Sscales);free(Shers);free(Shersds);free(SJ);
if(memsave==1){free(wantids);}
}	//end of discenv==1

////////

free(fixed);free(fixedsave);free(vars);free(vardiffs);free(varsds);free(scales);free(hers);free(hersds);free(shares);free(sharesds);
free(ZTY);free(ZTZ);free(ZTZ2);free(ZTZZTY);
free(AI);free(AI2);free(AI3);free(AIsave);free(BI);free(BIsave);free(J);free(JAI);free(JAIJT);free(JAIJTsave);
free(ZTVZ);free(ZTVZ2);free(ZTVZ3);free(PY);
for(k=0;k<total;k++){free(KPY[k]);free(PKPY[k]);}free(KPY);free(PKPY);free(traces);
free(ZTVY);free(thetas);free(thetasds);free(thetapvas);free(Yadj);
if(num_kins+num_regs>0)
{
for(k=0;k<num_kins+num_regs;k++){free(mg[k]);}free(mg);
for(k=0;k<num_kins+num_regs;k++){free(mg2[k]);}free(mg2);
}
if(num_regs>0||num_tops>0)
{
for(r=0;r<num_regs+4;r++){free(effects[r]);}free(effects);
}

if(shortcut==0){free(V);free(V2);free(VZ);free(VZZTVZ);free(P);}
if(shortcut==1){free(UTY);free(UTZ);free(D);free(BUTZ);free(H);free(HUTY);free(HKPY);}

if(Xtotal>0)
{
free(PX);free(XTPY);
if(shortcut==1)
{
free(UTX);free(DUTX);free(XTVX);free(XTVX2);free(XTVX3);free(F);free(FUTZ);free(FUTY);free(FUTX);free(HUTX);free(FKPY);}
}

if(discenv==1)
{
free(Sindexer);free(Snums);free(Sscales);free(SJ);free(Shers);free(Shersds);
if(memsave==1){free(wantids);}
}

}	//end of multi_reml

///////////////////////////

