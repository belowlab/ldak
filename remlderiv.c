/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You shouldt have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Compute REML derivatives

///////////////////////////

//want AI (minus second derivatives) and BI (first derivatives)
//require YTPKPY, YTPKPKPY and tr(PK) for each kinship or pair of kinships

if(Xtotal>0)	//will need PX and XTPY (for shortcut=1, actually get UTPX and XTPY)
{
if(shortcut==0)	//get both directly
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &Xtotal, &ns, &alpha, P, &ns, X, &ns, &beta, PX, &ns);
dgemv_("T", &ns, &Xtotal, &alpha, PX, &ns, Y, &one, &beta, XTPY, &one);
}

if(shortcut==1)	//PX = (invD - DUTXF - BUTZH) UTX, can get XTPY directly
{
for(i=0;i<ns;i++)
{
for(j=0;j<Xtotal;j++){PX[i+j*ns]=UTX[i+j*ns]/D[i];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "N", &Xtotal, &Xtotal, &ns, &alpha, F, &Xtotal, UTX, &ns, &beta, FUTX, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &Xtotal, &Xtotal, &alpha, DUTX, &ns, FUTX, &Xtotal, &beta, PX, &ns);
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_fixed, &Xtotal, &ns, &alpha, H, &num_fixed, UTX, &ns, &beta, HUTX, &num_fixed);
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &Xtotal, &num_fixed, &alpha, BUTZ, &ns, HUTX, &num_fixed, &beta, PX, &ns);
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &Xtotal, &alpha, PX, &ns, UTY, &one, &beta, XTPY, &one);
}
}

////////

//get KPYs

//get KPY for noise
for(i=0;i<ns;i++){KPY[0][i]=PY[i];}

for(k=0;k<num_kins;k++)	//get KPYs for kinships
{
if(shortcut==0)
{
if(memsave==1)
{
mkins[k]=malloc(sizeof(double)*ns*ns);
read_kins(kinstems[k], mkins[k], NULL, 1.0, ns, ids3, 2);
}
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, mkins[k], &ns, PY, &one, &beta, KPY[1+k], &one);
if(memsave==1){free(mkins[k]);}
}

if(shortcut==1)	//must have k=0
{
for(i=0;i<ns;i++){KPY[1][i]=E[i]*PY[i];}
}
}

for(r=0;r<num_regs;r++)	//get KPYs for regions
{
token=Xends[r]-Xstarts[r];
if(shortcut==0)
{
alpha=1.0/Xsums[r];beta=0.0;
dgemv_("N", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, XTPY+Xstarts[r], &one, &beta, KPY[1+num_kins+r], &one);
}

if(shortcut==1)
{
alpha=1.0/Xsums[r];beta=0.0;
dgemv_("N", &ns, &token, &alpha, UTX+Xstarts[r]*ns, &ns, XTPY+Xstarts[r], &one, &beta, KPY[1+num_kins+r], &one);
}
}	//end of r loop

////////

//get PKPYs

for(k=0;k<total;k++)
{
if(shortcut==0)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, KPY[k], &one, &beta, PKPY[k], &one);
}

if(shortcut==1)	//get (invD - DUTXF - BUTZH) KPY
{
//set PKPY = invD KPY
for(i=0;i<ns;i++){PKPY[k][i]=KPY[k][i]/D[i];}
if(Xtotal>0)	//using X, so subtract DUTXF KPY
{
alpha=1.0;beta=0.0;
dgemv_("N", &Xtotal, &ns, &alpha, F, &Xtotal, KPY[k], &one, &beta, FKPY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &Xtotal, &alpha, DUTX, &ns, FKPY, &one, &beta, PKPY[k], &one);
}
//subtract BUTZH KPY
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &ns, &alpha, H, &num_fixed, KPY[k], &one, &beta, HKPY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, BUTZ, &ns, HKPY, &one, &beta, PKPY[k], &one);
}
}	//end of k loop

////////

//get trace PKs

//get trace PK for noise
traces[0]=0;
if(shortcut==0)
{
for(i=0;i<ns;i++){traces[0]+=P[(size_t)i*ns+i];}
}

if(shortcut==1)	//trace P = trace (invD - DUTXF - BUTZH)
{
for(i=0;i<ns;i++)
{
traces[0]+=1.0/D[i];
for(j=0;j<Xtotal;j++){traces[0]-=DUTX[i+j*ns]*F[j+i*Xtotal];}
for(j=0;j<num_fixed;j++){traces[0]-=BUTZ[i+j*ns]*H[j+i*num_fixed];}
}
}

for(k=0;k<num_kins;k++)	//get trace PK for kinships
{
if(shortcut==0)
{
if(memsave==1)
{
mkins[k]=malloc(sizeof(double)*ns*ns);
read_kins(kinstems[k], mkins[k], NULL, 1.0, ns, ids3, 2);
}
traces[1+k]=0;
for(scount=0;scount<stotal;scount++){traces[1+k]+=P[scount]*mkins[k][scount];}
if(memsave==1){free(mkins[k]);}
}

if(shortcut==1)	//must have k=0 and tr(PK) = tr(U(...)UTUEUT) = tr(E(invD - DUTXF - BUTZH))
{
traces[1+k]=0;
for(i=0;i<ns;i++)
{
traces[1]+=E[i]/D[i];
for(j=0;j<Xtotal;j++){traces[1]-=E[i]*DUTX[i+j*ns]*F[j+i*Xtotal];}
for(j=0;j<num_fixed;j++){traces[1]-=E[i]*BUTZ[i+j*ns]*H[j+i*num_fixed];}
}
}
}	//end of k loop

for(r=0;r<num_regs;r++)	//get trace PK for regions
{
if(shortcut==0)
{
traces[1+num_kins+r]=0;
for(i=0;i<ns;i++)
{
for(j=Xstarts[r];j<Xends[r];j++){traces[1+num_kins+r]+=PX[i+j*ns]*X[i+j*ns]/Xsums[r];}
}
}

if(shortcut==1)
{
traces[1+num_kins+r]=0;
for(i=0;i<ns;i++)
{
for(j=Xstarts[r];j<Xends[r];j++){traces[1+num_kins+r]+=PX[i+j*ns]*UTX[i+j*ns]/Xsums[r];}
}
}
}

////////

//fill AI = -2nd deriv = .5*YTPKPKPY and BI = 1st deriv = .5*(YTPKPKY -trace(PK))

//start assuming all terms fixed (AI diag, BI=0)
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){AI[k+k2*total]=0;}
AI[k+k*total]=1;BI[k]=0;
}

//now fill up
for(k=0;k<total;k++)
{
if(fixed[k]<3)
{
//diagonals and B
gam2=0;for(i=0;i<ns;i++){gam2+=KPY[k][i]*PY[i];}
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[k][i]*KPY[k][i];}
AI[k+k*total]=.5*gam3;
BI[k]=.5*(gam2-traces[k]);

//off-diagonals
for(k2=0;k2<k;k2++)
{
if(fixed[k2]<3)
{
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[k][i]*KPY[k2][i];}
AI[k+k2*total]=.5*gam3;
AI[k2+k*total]=.5*gam3;
}}	//end of k2 free and k2 loop
}}	//end of k free and k loop

//for stability, scale AI so has trace one
sum=0;for(k=0;k<total;k++){sum+=AI[k+k*total];}
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){AI[k+k2*total]*=pow(sum,-1);}
}

//invert then set fixed to zero
(void)eigen_invert(AI, total, AI2, -1, AI3, 1);
for(k=0;k<total;k++)
{
if(fixed[k]>=3){AI[k+k*total]=0;}
}

//undo scaling of AI
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){AI[k+k2*total]*=pow(sum,-1);}
}

///////////////////////////

