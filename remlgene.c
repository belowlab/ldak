/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//REML with one gene - have not added Bayesian version

///////////////////////////

double gene_reml(int ns, int num_fixed, double *Y, double *Z, double *X, int Xtotal, double Xsum, double *Xnss, double *Xrhos, double *stats, double tol, int maxiter, double shrink, double strip)
{
int i, j, j2, k, count, one=1, lwork, info;
double nfree, scale, sum, sum2, value, alpha, beta, wkopt, *work;
double lambda, lambda2, lambdadiff, her, likenull, like, like2, likeold, diff;

double S1, S2, S3, T1, T2, T3, gam, deriv, dderiv, dderiv2;
double *ZTY, *ZTX, *ZTZ, *ZTZ2, *ZTZ3, detZTZ, *ZTZZTY, *ZTZZTX;
double YTCY, *XTCY, *XTCX, *U, *E, *D;

double begin[19]={.00005,.0001,.0002,.0005,.001,.002,.005,.01,.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9};


//set nfree and scale (how much to scale up)
if(Xnss==NULL)
{
nfree=ns-num_fixed;
scale=1;
}
else
{
sum=0;for(j=0;j<Xtotal;j++){sum+=Xnss[j];}
nfree=sum/Xtotal;
scale=nfree/ns;
}

//allocate variables - will always have Xtotal>0

ZTY=malloc(sizeof(double)*num_fixed);
ZTX=malloc(sizeof(double)*num_fixed*Xtotal);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZZTY=malloc(sizeof(double)*num_fixed);
ZTZZTX=malloc(sizeof(double)*num_fixed*Xtotal);

XTCY=malloc(sizeof(double)*Xtotal);
XTCX=malloc(sizeof(double)*Xtotal*Xtotal);
U=malloc(sizeof(double)*Xtotal*Xtotal);
E=malloc(sizeof(double)*Xtotal);
D=malloc(sizeof(double)*Xtotal);

//fill some variables

//get ZTY, ZTX, ZTZ
alpha=scale;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_fixed, &Xtotal, &ns, &alpha, Z, &ns, X, &ns, &beta, ZTX, &num_fixed);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);

//invert ZTZ and get ZTZZTY and ZTZZTX
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, ZTZZTY, &one);
dgemm_("N", "N", &num_fixed, &Xtotal, &num_fixed, &alpha, ZTZ, &num_fixed, ZTX, &num_fixed, &beta, ZTZZTX, &num_fixed);

//get YTCY (when using summary statistics, this equals var(Y)*nfree) and likenull
YTCY=0;
for(i=0;i<ns;i++){YTCY+=scale*pow(Y[i],2);}
for(j=0;j<num_fixed;j++){YTCY-=ZTY[j]*ZTZZTY[j];}
likenull=-.5*(nfree+nfree*log(2*M_PI*YTCY/nfree)+detZTZ);

//get XTCX (when using summary statistics, this equals covar(X)*nfree)
alpha=scale;beta=0.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &ns, &alpha, X, &ns, X, &ns, &beta, XTCX, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &num_fixed, &alpha, ZTX, &num_fixed, ZTZZTX, &num_fixed, &beta, XTCX, &Xtotal);

//get XTCY - only use summary statistics for real analysis (not when permuting)
if(Xnss==NULL||stats==NULL)
{
alpha=scale;beta=0.0;
dgemv_("T", &ns, &Xtotal, &alpha, X, &ns, Y, &one, &beta, XTCY, &one);
alpha=-1.0;beta=1.0;
dgemv_("T", &num_fixed, &Xtotal, &alpha, ZTX, &num_fixed, ZTZZTY, &one, &beta, XTCY, &one);
}
else
{
for(j=0;j<Xtotal;j++){XTCY[j]=Xrhos[j]*pow(XTCX[j+j*Xtotal]*YTCY,.5);}
}

if(shrink!=1)	//alter XTCX (must be using summary statistics, so XTCX=XTX)
{
if(shrink<1)	//deflate off diagonal terms
{
for(j=0;j<Xtotal;j++)
{
for(j2=j+1;j2<Xtotal;j2++){XTCX[j+j2*Xtotal]*=shrink;XTCX[j2+j*Xtotal]*=shrink;}
}
}
if(shrink>1)	//inflate diagonal terms
{
for(j=0;j<Xtotal;j++){XTCX[j+j*Xtotal]*=shrink;}
}
}

//decomp XTCX
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){U[j+j2*Xtotal]=XTCX[j+j2*Xtotal];}
}
lwork=-1;
dsyev_("V", "U", &Xtotal, U, &Xtotal, E, &wkopt, &lwork, &info);
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &Xtotal, U, &Xtotal, E, work, &lwork, &info);
free(work);

if(strip>0)	//set end proportion of E to zero (must be using summary statistics, so XTCX=XTX)
{
sum=0;for(j=0;j<Xtotal;j++){sum+=E[j];}
sum2=0;
for(j=Xtotal-1;j>=0;j--)
{
sum2+=E[j];
if(sum2/sum>1-strip){E[j]=0;}
}
}

//get D = UTXTCY
alpha=1.0;beta=0.0;
dgemv_("T", &Xtotal, &Xtotal, &alpha, U, &Xtotal, XTCY, &one, &beta, D, &one);

//////////////////////////

//likelihood -.5 [nfree (1 + log(2pi gamma/nfree)) + log |E+lambda| - Xtotal log(lambda) + detZTZ]
//where gamma is YTCY - DT (E + lambda)^-1 D

//test some heritabilities and start with best fitting
for(k=0;k<19;k++)
{
her=begin[k];
lambda2=Xsum*(1-her)/her;

S1=0;T1=0;
for(j=0;j<Xtotal;j++)
{
S1+=pow(D[j],2)/(E[j]+lambda2);
if(E[j]+lambda2>0){T1+=log(E[j]+lambda2);}
}
gam=YTCY-S1;
like2=-.5*(nfree+nfree*log(2*M_PI*gam/nfree)+T1-Xtotal*log(lambda2)+detZTZ);

if(k==0){lambda=lambda2;like=like2;}
if(like2>like){lambda=lambda2;like=like2;}
}

//now iterate
count=0;
while(1)
{
//get derivatives (already have likelihood)
S1=0;S2=0;S3=0;T1=0;T2=0;T3=0;
for(j=0;j<Xtotal;j++)
{
S1+=pow(D[j],2)/(E[j]+lambda);
S2+=pow(D[j],2)*pow(E[j]+lambda,-2);
S3+=pow(D[j],2)*pow(E[j]+lambda,-3);
if(E[j]+lambda>0){T1+=log(E[j]+lambda);}
T2+=1/(E[j]+lambda);
T3+=pow(E[j]+lambda,-2);
}
gam=YTCY-S1;
if(gam<0){printf("Error, negative gamma, please tell Doug %f %f %f | %f %f %f | %f\n\n", S1, S2, S3, T1, T2, T3, gam);exit(1);}

deriv=-.5*nfree/gam*S2-.5*T2+.5*Xtotal/lambda;
dderiv=-.5*nfree*pow(gam,-2)*pow(S2,2)+nfree/gam*S3+.5*T3-.5*Xtotal*pow(lambda,-2);

//see if breaking
if(count>0)	//test for convergence
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;
if(count==maxiter){printf("Warning, REML failed to converge within %d iterations\n", maxiter);break;}

//get proposed move, ensuring remains positive
lambdadiff=-deriv/dderiv;
if(lambda+lambdadiff<=0){lambdadiff=1e-6-lambda;}

value=1;
while(value>0.0001)
{
//get likelihood based on moving value*lambdadiff
lambda2=lambda+value*lambdadiff;
S1=0;T1=0;
for(j=0;j<Xtotal;j++)
{
S1+=pow(D[j],2)/(E[j]+lambda2);
if(E[j]+lambda2>0){T1+=log(E[j]+lambda2);}
}
gam=YTCY-S1;
like2=-.5*(nfree+nfree*log(2*M_PI*gam/nfree)+T1-Xtotal*log(lambda2)+detZTZ);

if(like2>like-tol)	//accept move
{lambda=lambda2;like=like2;break;}
else	//next turn try smaller move
{value*=.5;}
}

count++;
}	//end of while loop

if(stats!=NULL)	//get many statistics
{
her=Xsum/(Xsum+lambda);
dderiv2=pow(Xsum,2)*pow(her,-4)*dderiv+2*Xsum*pow(her,-3)*deriv;

stats[0]=her;
if(dderiv2<0){stats[1]=pow(-dderiv2,-.5);}
else{stats[1]=-9999;}
stats[2]=likenull;
stats[3]=like;
stats[4]=2*(like-likenull);
if(stats[4]>0){stats[5]=.5*erfc(pow(stats[4],.5)*M_SQRT1_2);}
else{stats[5]=.75;}
}

free(ZTY);free(ZTX);free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTZZTY);free(ZTZZTX);
free(XTCY);free(XTCX);free(U);free(E);free(D);

return(2*(like-likenull));
}	//end of reml_lite

///////////////////////////

