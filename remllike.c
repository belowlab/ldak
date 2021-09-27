/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You shouldt have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Compute REML likelihood

///////////////////////////

//V has form var0 I + var1 K1 + var2 K2 + ... + var(1+nk) X1X1T + var(1+nk+1) X2X2T + ... 
//when shortcut=1, will have nk=0 or nk=1 and write V as var0 I + var1 UEUT + XCXT

if(shortcut==0)	//get invV directly
{
//load up V
for(scount=0;scount<stotal;scount++){V[scount]=0;}
for(i=0;i<ns;i++){V[(size_t)i*ns+i]=vars[0];}

for(k=0;k<num_kins;k++)	//add on varK for kinships
{
if(vars[1+k]!=0)
{
if(memsave==0)	//have kinships saved in mkins
{
for(scount=0;scount<stotal;scount++){V[scount]+=vars[1+k]*mkins[k][scount];}
}
else{read_kins(kinstems[k], V, NULL, vars[1+k], ns, ids3, 3);}
}
}

for(r=0;r<num_regs;r++)	//and for regions
{
if(vars[1+num_kins+r]!=0)
{
token=Xends[r]-Xstarts[r];
alpha=vars[1+num_kins+r]/Xsums[r];beta=1.0;
dgemm_("N", "T", &ns, &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, V, &ns);
}
}

//invert V - save V in case ldlt fails
V3=malloc(sizeof(double)*ns*ns);
for(scount=0;scount<stotal;scount++){V3[scount]=V[scount];}
detV=ldlt_invert(V, ns, -1, NULL, &info, 1);
if(info!=0)	//eigen will definitely work, so can use V3 as work space
{
for(scount=0;scount<stotal;scount++){V[scount]=V3[scount];}
detV=eigen_invert(V, ns, V2, -1, V3, 1);
}
free(V3);
}	//end of shortcut=0

if(shortcut==1)	//get invV = U (invD - invD UTX F) UT, where F= (t(UTX) invDUTX + invC)^-1 t(invDUTX)
{
//find D and detD where var0 I + var1 K1 = UDUT (D = vars0 I or D = vars0 I + vars1 E)
detD=0;
for(i=0;i<ns;i++)
{
if(num_kins==1){D[i]=vars[0]+vars[1]*E[i];}
else{D[i]=vars[0];}
if(D[i]!=0){detD+=log(fabs(D[i]));}
}

if(Xtotal>0)	//using X, so find F
{
//get DUTX = invD UTX
for(r=0;r<num_regs;r++)
{
for(j=Xstarts[r];j<Xends[r];j++)
{
if(vars[1+num_kins+r]!=0)	//fill column
{
for(i=0;i<ns;i++){DUTX[i+j*ns]=UTX[i+j*ns]/D[i];}
}
else	//blank column
{
for(i=0;i<ns;i++){DUTX[i+j*ns]=0;}
}
}}

//get XTVX = t(UTX) DUTX + invC, where C = vars/Xsums (recall for regions with var=0, DUTX is blank)
//detV = |D + UTX invC XTU| = |D| |C| |XUT invD UTX + invC| = |D| |C| |XTVX|
alpha=1.0;beta=0.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &ns, &alpha, UTX, &ns, DUTX, &ns, &beta, XTVX, &Xtotal);

detC=0;
for(r=0;r<num_regs;r++)
{
for(j=Xstarts[r];j<Xends[r];j++)
{
if(vars[1+num_kins+r]!=0)	//add invC to diagonal and contribution to detC
{
XTVX[j+j*Xtotal]+=Xsums[r]/vars[1+num_kins+r];
detC+=log(fabs(vars[1+num_kins+r]/Xsums[r]));
}
else	//blank row except for 1 on diagonal (no contribution to detC)
{
for(j2=0;j2<Xtotal;j2++){XTVX[j+j2*Xtotal]=0;}
XTVX[j+j*Xtotal]=1;
}
}}	//end of j and r loops

//get F = invXTVX t(DUTX) - save XTVX in case ldlt fails
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){XTVX3[j+j2*Xtotal]=XTVX[j+j2*Xtotal];}
}

//first set F = t(DUTX)
for(j=0;j<Xtotal;j++)
{
for(i=0;i<ns;i++){F[j+i*Xtotal]=DUTX[i+j*ns];}
}
detXTVX=ldlt_invert(XTVX, Xtotal, ns, F, &info, 1);
if(info!=0)	//inversion failed, so use eigen (XTVX saved in XTVX3, F will still equal DUTX)
{detXTVX=eigen_invert(XTVX3, Xtotal, XTVX2, ns, F, 1);}
detV=detD+detC+detXTVX;
}	//end of using X
else	//not using X, so F not used and V = UDUT
{detV=detD;}
}	//end of shortcut=1

////////

//now want P = invV - invVZ invZTVZ ZTinvV or U (invD-DUTXF-BUTZH) UT, H = invZTVZ t(BUTZ)

if(shortcut==0)	//get invVZ then ZTVZ = ZT invVZ
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_fixed, &ns, &alpha, V, &ns, Z, &ns, &beta, VZ, &ns);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, VZ, &ns, &beta, ZTVZ, &num_fixed);
}

if(shortcut==1)	//get BUTZ = UTinvV UTZ, then ZTVZ = ZTU BUTZ
{
//set BUTZ = invD UTZ
for(i=0;i<ns;i++)
{
for(j=0;j<num_fixed;j++){BUTZ[i+j*ns]=UTZ[i+j*ns]/D[i];}
}
if(Xtotal>0)	//using X, so subtract DUTZ F UTZ
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &Xtotal, &num_fixed, &ns, &alpha, F, &Xtotal, UTZ, &ns, &beta, FUTZ, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &num_fixed, &Xtotal, &alpha, DUTX, &ns, FUTZ, &Xtotal, &beta, BUTZ, &ns);
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, UTZ, &ns, BUTZ, &ns, &beta, ZTVZ, &num_fixed);
}

//invert ZTVZ
detZTVZ=eigen_invert(ZTVZ, num_fixed, ZTVZ2, -1, ZTVZ3, 1);

if(shortcut==0)	//P = invV - invVZ invZTVZ ZTinvV
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_fixed, &num_fixed, &alpha, VZ, &ns, ZTVZ, &num_fixed, &beta, VZZTVZ, &ns);
for(scount=0;scount<stotal;scount++){P[scount]=V[scount];}
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &ns, &ns, &num_fixed, &alpha, VZZTVZ, &ns, VZ, &ns, &beta, P, &ns);
}

if(shortcut==1)	//P = U (invD-DUTXF-BUTZH) UT, where H = invZTVZ t(BUTZ)
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_fixed, &ns, &num_fixed, &alpha, ZTVZ, &num_fixed, BUTZ, &ns, &beta, H, &num_fixed);
}

////////

//get PY (for shortcut=1, actually get UTPY) and gam = YTPY

if(shortcut==0)	//get direct
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, Y, &one, &beta, PY, &one);
gam=0;for(i=0;i<ns;i++){gam+=PY[i]*Y[i];}
}

if(shortcut==1) //get UTPY = (invD - DUTXF - BUTZH) UTY
{
//set PY = invD UTY
for(i=0;i<ns;i++){PY[i]=UTY[i]/D[i];}
if(Xtotal>0)	//using X, so subtract DUTXF UTY
{
alpha=1.0;beta=0.0;
dgemv_("N", &Xtotal, &ns, &alpha, F, &Xtotal, UTY, &one, &beta, FUTY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &Xtotal, &alpha, DUTX, &ns, FUTY, &one, &beta, PY, &one);
}
//subtract BUTZH UTY
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &ns, &alpha, H, &num_fixed, UTY, &one, &beta, HUTY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, BUTZ, &ns, HUTY, &one, &beta, PY, &one);

gam=0;for(i=0;i<ns;i++){gam+=PY[i]*UTY[i];}
}

///////////////////////////

