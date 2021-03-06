/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You shouldt have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Decompositions

///////////////////////////

double eigen_invert(double *mat, int length, double *mat2, int ncol, double *mat3, int type)
{
//ncol=0 - get "cholesky", ncol=-1 - get inverse, ncol>0 - solve mat X = mat3
//ncol=0 - mat3 NULL, ncol=-1 - mat3 workspace size(mat), ncol>0 - mat3 = RHS
//type=0 - quiet, type=1 - complain
int i, j, count, lwork, info;
double det, value, alpha, beta, wkopt, *work, *mat4;


if(length==0){return(0);}

lwork=-1;
dsyev_("V", "U", &length, mat, &length, mat2, &wkopt, &lwork, &info);
if(info!=0)
{printf("Error, eigen priming failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);

dsyev_("V", "U", &length, mat, &length, mat2, work, &lwork, &info);
if(info!=0)
{printf("Error, eigen decomp failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
free(work);

//get determinant
det=0;count=0;
for(i=0;i<length;i++)
{
if(fabs(mat2[i])>=0.000001){det+=log(fabs(mat2[i]));}
else{det+=log(0.000001);}
if(mat2[i]<=-0.000001){count++;}
}
if(count>1&&type==1)
{printf("Warning, %d eigenvalues were negative\n", count);}

////////

if(ncol==0)	//solve XXT = mat, return in mat - X = UE^.5
{
for(j=0;j<length;j++)
{
if(mat2[j]>0.000001)
{
value=pow(mat2[j],.5);
for(i=0;i<length;i++){mat[(size_t)j*length+i]=mat[(size_t)j*length+i]*value;}
}
else
{
for(i=0;i<length;i++){mat[(size_t)j*length+i]=0;}
}
}
}

if(ncol==-1)	//solve X mat = I, return in mat - X = (UE^-.5) t(UE^-.5) or UE^-1UT
{
if(count==0)	//no negative values, so can use square root
{
for(j=0;j<length;j++)
{
if(mat2[j]>0.000001)
{
value=pow(mat2[j],.5);
for(i=0;i<length;i++){mat3[(size_t)j*length+i]=mat[(size_t)j*length+i]/value;}
}
else
{
for(i=0;i<length;i++){mat3[(size_t)j*length+i]=0;}
}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &length, &length, &length, &alpha, mat3, &length, mat3, &length, &beta, mat, &length);
}	//end of count=0
else	//there are negative values, so need extra matrix
{
mat4=malloc(sizeof(double)*length*length);
for(i=0;i<length;i++)
{
for(j=0;j<length;j++)
{
mat3[(size_t)j*length+i]=mat[(size_t)j*length+i];
if(fabs(mat2[i])>=0.000001){mat4[(size_t)j*length+i]=mat[(size_t)j*length+i]/mat2[j];}
else{mat4[(size_t)j*length+i]=0;}
}}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &length, &length, &length, &alpha, mat4, &length, mat3, &length, &beta, mat, &length);
free(mat4);
}
}	//end of ncol=-1

if(ncol>0)	//solve mat X = mat3 - X = UE^-1UT mat3
{
mat4=malloc(sizeof(double)*length*ncol);
alpha=1.0;beta=0.0;
dgemm_("T", "N", &length, &ncol, &length, &alpha, mat, &length, mat3, &length, &beta, mat4, &length); 
for(i=0;i<length;i++)
{
for(j=0;j<ncol;j++)
{
if(fabs(mat2[i])>=0.000001){mat4[(size_t)j*length+i]=mat4[(size_t)j*length+i]/mat2[i];}
else{mat4[(size_t)j*length+i]=0;}
}
}
dgemm_("N", "N", &length, &ncol, &length, &alpha, mat, &length, mat4, &length, &beta, mat3, &length);
free(mat4);
}

return(det);
}	//end of eigen_invert

////////

void eigen_strip(double *mat, int start, int end, int total, double strip)
{
int i, j, length, lwork, info;
double *U, *E, sum, sum2, value, alpha, beta, wkopt, *work;


length=end-start;
U=malloc(sizeof(double)*length*length);
E=malloc(sizeof(double)*length);
for(i=0;i<length;i++)
{
for(j=0;j<length;j++){U[(size_t)j*length+i]=mat[(start+i)+(start+j)*total];}
}

lwork=-1;
dsyev_("V", "U", &length, U, &length, E, &wkopt, &lwork, &info);
if(info!=0)
{printf("Error, eigen priming failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &length, U, &length, E, work, &lwork, &info);
if(info!=0)
{printf("Error, eigen decomp failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
free(work);

sum=0;for(i=0;i<length;i++){sum+=E[i];}
sum2=0;
for(i=length-1;i>=0;i--)
{
sum2+=E[i];
if(sum2/sum>1-strip){E[i]=0;}
}

for(j=0;j<length;j++)
{
value=pow(E[j],.5);
for(i=0;i<length;i++){U[(size_t)j*length+i]*=value;}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &length, &length, &length, &alpha, U, &length, U, &length, &beta, mat+start+start*total, &total);

free(U);free(E);
}

///////////////////////////

double ldlt_invert(double *mat, int length, int ncol, double *mat2, int *info, int type)
{
//ncol=-1 - get inverse, ncol>0 - solve mat X = mat2
//type=0 - quiet, type=1 - complain
int i, j, *ipiv, lwork;
double wkopt, *work, det;


if(length>0)
{
ipiv=malloc(sizeof(int)*length);
lwork=-1;
dsytrf_("U", &length, mat, &length, ipiv, &wkopt, &lwork, info);
if(*info!=0)
{printf("Error, LDLT priming failed; please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsytrf_("U", &length, mat, &length, ipiv, work, &lwork, info);
free(work);
if(*info!=0)
{
if(type==1){printf("Warning, LDLT failed (info %d, length %d), resorting to eigen-decomposition\n", *info, length);}
free(ipiv);
}
}

//get determinant
det=0;
for(i=0;i<length;i++)
{
if(fabs(mat[(size_t)i*length+i])>=0.000001){det+=log(fabs(mat[(size_t)i*length+i]));}
else{det+=log(0.000001);}
}

////////

if(ncol==-1)	//put inverse into mat
{
work=malloc(sizeof(double)*length);
dsytri_("U", &length, mat, &length, ipiv, work, info);
if(*info!=0)
{printf("Error LDLT invert failed, please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
for(i=0;i<length;i++)
{
for(j=0;j<i;j++){mat[(size_t)j*length+i]=mat[j+i*length];}
}
free(work);
}

if(ncol>0)	//solve
{
dsytrs_("U", &length, &ncol, mat, &length, ipiv, mat2, &length, info);
if(*info!=0)
{printf("Error LDLT solve failed, please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
}

free(ipiv);

return(det);
}

///////////////////////////

