/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You shouldt have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for reading and writing kinships

///////////////////////////

int read_kins(char *kinstem, double *kins, float *kins_single, double denom, int ns, char **ids3, int type)
//type=0 - just read kins, type=1 - add kins, type=2 - read kins quiet, type=3 - add kins quiet
//in all cases, actually using kins*denom
{
int i, i2, count, count2, current;
int *indexer, *indexer2;

float *datatemp;
char **wantids;

char filename[500];
FILE *input;


if(type==0||type==1){printf("Reading kinship matrix with stem %s\n",kinstem);}

//first get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*ns);
indexer2=malloc(sizeof(int)*ns);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);

if(count2!=ns){printf("Error finding indexes read_kin, please tell doug\n\n");exit(1);}

//open kinship and check size
sprintf(filename, "%s.grm.bin", kinstem);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*count*(count+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*count*(count+1)/2, ftello(input));exit(1);}

fseeko(input, 0, SEEK_SET);
current=0;

//read one "row" at a time
datatemp=malloc(sizeof(float)*count);
for(i=0;i<ns;i++)
{
if((type==0||type==1)&&count>40000&&i%20000==0)
{printf("Reading kinships for Sample %d out of %d\n", i+1, ns);}

if(indexer[i]!=current)	//get to start of Row indexer[i]
{fseeko(input, (off_t)sizeof(float)*indexer[i]*(indexer[i]+1)/2, SEEK_SET);}
if(fread(datatemp, sizeof(float), indexer[i]+1, input)!=indexer[i]+1)
{printf("Error reading Row %d of %s\n\n", indexer[i]+1, filename);exit(1);};
current=indexer[i]+1;

if(kins!=NULL)	//saving as doubles
{
if(type==0||type==2)	//just read kins
{
//read off diagonals, then diagonal - this works because indexer monotonic
for(i2=0;i2<i;i2++)
{
kins[(size_t)indexer2[i2]*ns+indexer2[i]]=datatemp[indexer[i2]]*denom;
kins[(size_t)indexer2[i]*ns+indexer2[i2]]=datatemp[indexer[i2]]*denom;
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]=datatemp[indexer[i]]*denom;
}
if(type==1||type==3)	//add kins
{
//read off diagonals, then diagonal - this works because indexer monotonic
for(i2=0;i2<i;i2++)
{
kins[(size_t)indexer2[i2]*ns+indexer2[i]]+=datatemp[indexer[i2]]*denom;
kins[(size_t)indexer2[i]*ns+indexer2[i2]]+=datatemp[indexer[i2]]*denom;
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]+=datatemp[indexer[i]]*denom;
}
}
else	//saving as floats
{
if(type==0||type==2)	//just read kins
{
//read off diagonals, then diagonal - this works because indexer monotonic
for(i2=0;i2<i;i2++)
{
kins_single[(size_t)indexer2[i2]*ns+indexer2[i]]=datatemp[indexer[i2]]*denom;
kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]=datatemp[indexer[i2]]*denom;
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]=datatemp[indexer[i]]*denom;
}
if(type==1||type==3)	//add kins
{
//read off diagonals, then diagonal - this works because indexer monotonic
for(i2=0;i2<i;i2++)
{
kins_single[(size_t)indexer2[i2]*ns+indexer2[i]]+=datatemp[indexer[i2]]*denom;
kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]+=datatemp[indexer[i2]]*denom;
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]+=datatemp[indexer[i]]*denom;
}
}
}	//end of i loop
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(datatemp);

return(0);
}	//end of read_kins

////////

double read_kin_trace(char *kinstem, int ns, char **ids3, int type)
//type=0 - normal, type=1 - quiet
{
int i, count, count2;
int *indexer, *indexer2;
double sum;

char **wantids;
float readfloat;

char filename[500];
FILE *input;


if(type==0){printf("Reading trace of kinship matrix with stem %s\n",kinstem);}

//first get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count2!=ns){printf("Error find indexes read_kin_trace, please tell doug\n\n");exit(1);}

//open kinship and check size
sprintf(filename, "%s.grm.bin", kinstem);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*count*(count+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*count*(count+1)/2, ftello(input));exit(1);}

//now tally trace
sum=0;
for(i=0;i<ns;i++)
{
//get to ith diagonal element
fseeko(input, (off_t)sizeof(float)*indexer[i]*(indexer[i]+3)/2, SEEK_SET);
if(fread(&readfloat, sizeof(float), 1, input)!=1)
{printf("Error reading Row %d of %s\n\n", i+1, filename);exit(1);};
sum+=readfloat;
}
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);

return(sum/ns);
}

///////////////////////////

int read_details(char *filename, char **kpreds, int *kindex, double *kcentres, double *kmults, double *kweights, char *kal1, char *kal2, int length)
{
int j;

int readint;
char readchar, *rs;

FILE *input;

rs=malloc(sizeof(char)*10000000);


if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}

for(j=0;j<length;j++)
{
if(fscanf(input, "%s %d %lf %lf %lf %c %c ", rs, &readint, kcentres+j, kmults+j, kweights+j, kal1+j, kal2+j)!=7)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation\n\n", j+2, filename);exit(1);}
copy_string(kpreds,j,rs);
kindex[j]=readint-1;
}
fclose(input);

free(rs);
return(0);
}	//end of read_details

///////////////////////////

int read_kinsgz(char *kinstem, float *kins_single, int ns, char **ids3)
{
int i, i2, readi, readi2, count, count2, found, size, offset;
int *indexer, *indexer2;

double denom, *datatemp;
char **wantids, *gzbuffer;

char filename[500];
gzFile inputgz;


//size of buffer
size=1000;

printf("Reading %s.grm.id and %s.grm.gz\n\n", kinstem, kinstem);

//get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, filename, NULL, NULL, NULL, 3);

//get ready to read
sprintf(filename, "%s.grm.gz", kinstem);
if((inputgz=gzopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

gzbuffer=malloc(sizeof(char)*size);

//read set of i rows at a time
datatemp=malloc(sizeof(double)*count);
found=0;count2=0;
for(i=0;i<count;i++)
{
for(i2=0;i2<=i;i2++)
{
if(gzgets(inputgz,gzbuffer,size)==NULL)
{printf("Error reading %s; there is no Row %d\n\n", filename, count2);exit(1);}
if(strlen(gzbuffer)==size-1)
{printf("Error reading %s; Row %d is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", filename, count2+1, (int)strlen(gzbuffer), size);exit(1);}

//check starting elements
if(sscanf(gzbuffer, "%d %d %n", &readi, &readi2, &offset)!=2)
{printf("Error reading start of Row %d of %s\n\n", count2+1, filename);exit(1);}
if(readi!=i+1||readi2!=i2+1)
{printf("Error reading %s; Row %d starts %d %d (not %d %d)\n\n", filename, count2+1, readi, readi2, i+1, i2+1);exit(1);}

//then save counts and kinships
if(sscanf(gzbuffer+offset, "%lf %lf ", &denom, datatemp+i2)!=2)
{printf("Error reading end of Row %d of %s\n\n", count2+1, filename);exit(1);}
count2++;
}

if(i==indexer[found])	//using i, so save
{
for(i2=0;i2<found;i2++)
{
kins_single[(size_t)indexer2[i2]*ns+indexer2[found]]=datatemp[indexer[i2]];
kins_single[(size_t)indexer2[found]*ns+indexer2[i2]]=datatemp[indexer[i2]];
}
kins_single[(size_t)indexer2[found]*ns+indexer2[found]]=datatemp[indexer[found]];
found++;
if(found==ns){break;}
}
}
gzclose(inputgz);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(datatemp);

return(0);
}	//end of read_kinfilegz

////////

int read_kinsraw(char *kinstem, float *kins_single, int ns, char **ids3)
{
int i, i2, count, count2, found;
int *indexer, *indexer2;

double *datatemp; 
char **wantids;

char filename[500];
FILE *input;


printf("Reading %s.grm.id and %s.grm.raw\n\n", kinstem, kinstem);

//get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, filename, NULL, NULL, NULL, 3);
if(count2<ns)
{printf("Error 83D, please tell Doug %d %d\n", count2, ns);exit(1);}

//get ready to read
sprintf(filename, "%s.grm.raw", kinstem);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

datatemp=malloc(sizeof(double)*count);
found=0;
for(i=0;i<count;i++)
{
for(i2=0;i2<count;i2++)
{
if(fscanf(input, "%lf\n", datatemp+i2)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", i2+1, i+1, filename);exit(1);}
}

if(i==indexer[found])	//using i, so save
{
for(i2=0;i2<ns;i2++){kins_single[(size_t)indexer2[found]*ns+indexer2[i2]]=datatemp[indexer[i2]];}
found++;
if(found==ns){break;}
}
}
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(datatemp);

return(0);
}	//end of read_kinfileraw

///////////////////////////

int write_kins(char *outfile, double *kins, float *kins_single, int ns, char **ids1, char **ids2, int kindetails, char **preds, int *keeppreds_use, double *centres, double *mults, double *weights, char *al1, char *al2, int length, char *datafile, double power, int kingz, int kinraw, int type)
//type=0 - quiet, type=1 - normal, type=2 - null for association testing, type=3 - trun/pca/square/gxemm kinships
{
int i, i2, j, count;
float *datatemp;
double sum, sumsq, value, denom;

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], cmd[500];
FILE *output, *output2, *output3, *output4, *output5;


if(type!=3)	//get denom
{
sum=0;for(i=0;i<ns;i++)
{
if(kins!=NULL){sum+=kins[(size_t)i*ns+i];}
else{sum+=kins_single[(size_t)i*ns+i];}
}
denom=sum/ns;
}
else{denom=1;}

//write kins to grm.bin
sprintf(filename,"%s.grm.bin", outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

datatemp=malloc(sizeof(float)*ns);
sum=0;sumsq=0;
for(i2=0;i2<ns;i2++)
{
for(i=0;i<=i2;i++)
{
if(kins!=NULL){datatemp[i]=kins[(size_t)i2*ns+i]/denom;}
else{datatemp[i]=kins_single[(size_t)i2*ns+i]/denom;}
}
fwrite(datatemp, sizeof(float), i2+1, output);
for(i=0;i<i2;i++){sum+=datatemp[i];sumsq+=pow(datatemp[i],2);}
}
fclose(output);
free(datatemp);

if(kindetails==1)	//write details and adjustments
{
sprintf(filename2,"%s.grm.details", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor Index Centre Scaling Weight A1 A2\n");
count=0;value=0;
for(j=0;j<length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output2, "%s %d %.6f %.6f %.6f %c %c\n", preds[j], keeppreds_use[j]+1, centres[j], mults[j], weights[j], al1[j], al2[j]);
count++;
value+=weights[j];
}
}
fclose(output2);

sprintf(filename3,"%s.grm.adjust", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

if(datafile!=NULL){fprintf(output3, "Datafile %s\n", datafile);}
else{fprintf(output3, "Multiple Datafiles\n");}
if(power!=-9999){fprintf(output3, "Power %.6f\n", power);}
else{fprintf(output3, "Multiple Powers\n");}

fprintf(output3, "Num_Preds %d\nSum_Weights %.4f\nDenominator %.4f\nOff_Diag_Variance %e\n", count, value, denom, 2*sumsq/ns/(ns-1)-pow(2*sum/ns/(ns-1),2));
fclose(output3);
}

//write IDs
sprintf(filename4,"%s.grm.id", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(i=0;i<ns;i++){
fprintf(output4, "%s %s\n", ids1[i], ids2[i]);}
fclose(output4);

if(type==1)
{
printf("Kinship matrix saved in %s and %s", filename, filename4);
if(kindetails==1)
{printf(", with details saved in %s and %s", filename2, filename3);}
printf("\n");
}

////////

if(kingz==1)
{
sprintf(filename5, "%s.grm", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
value=denom;
if(value<1){value=1;}
for(i=0;i<ns;i++)
{
for(i2=0;i2<=i;i2++)
{
if(kins!=NULL){fprintf(output5, "%d %d %d %.6f\n", i+1, i2+1, (int)value, kins[(size_t)i*ns+i2]/denom);}
else{fprintf(output5, "%d %d %d %.6f\n", i+1, i2+1, (int)value, kins_single[(size_t)i*ns+i2]/denom);}
}
}
fclose(output5);

//try to compress the kinship file
sprintf(cmd, "gzip -f %s", filename5);
system(cmd);
if(type==1){printf("Gzipped version saved in %s.gz\n", filename5);}
}

if(kinraw==1)
{
sprintf(filename5, "%s.grm.raw", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++)
{
if(kins!=NULL){fprintf(output5,"%.6f ", kins[(size_t)i*ns+i2]/denom);}
else{fprintf(output5,"%.6f ", kins_single[(size_t)i*ns+i2]/denom);}
}
fprintf(output5,"\n");
}
fclose(output5);
if(type==1){printf("Plain text version saved in %s\n", filename5);}
}

if(type==1){printf("\n");}

return(0);
}	//end of write_kins

///////////////////////////

int write_eigen(char *outfile, double *U, double *E, int ns, char **ids1, char **ids2, char *kinstem, int eigenraw)
{
int i, i2;
float *datatemp;

char filename[500], filename2[500], filename3[500], filename4[500];
FILE *output, *output2, *output3, *output4;


sprintf(filename,"%s.eigen.bin",outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}

datatemp=malloc(sizeof(float)*ns);
for(i=0;i<ns;i++){datatemp[i]=E[i];}
fwrite(datatemp, sizeof(float), ns, output);
for(i2=0;i2<ns;i2++)
{
for(i=0;i<ns;i++){datatemp[i]=U[(size_t)i2*ns+i];}
fwrite(datatemp, sizeof(float), ns, output);
}
fclose(output);
free(datatemp);

sprintf(filename2,"%s.eigen.id",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(i=0;i<ns;i++){fprintf(output2,"%s %s\n", ids1[i], ids2[i]);}
fclose(output2);

sprintf(filename3,"%s.eigen.root",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
fprintf(output3,"Kinship %s\n", kinstem);
fclose(output3);

if(eigenraw==1)
{
sprintf(filename4,"%s.eigen.raw",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
for(i=ns-1;i>=0;i--){fprintf(output4, "%.6f ", E[i]);}
fprintf(output4, "\n");
for(i2=ns-1;i2>=0;i2--)
{
for(i=0;i<ns;i++){fprintf(output4, "%.6f ", U[(size_t)i2*ns+i]);}
fprintf(output4, "\n");
}
fclose(output4);
}

printf("Eigen-decomposition saved in %s, %s and %s", filename, filename2, filename3);
if(eigenraw==1)
{printf(", with plain text version saved in %s", filename4);}
printf("\n\n");

return(0);
}	//end of write_eigen

////////

int read_eigens(char *eigenfile, double *U, double *E, int ns, char **ids3, int axes)
{
int i, i2, count, count2;
int *indexer;

float *datatemp;
char **wantids;

char filename[500];
FILE *input;


printf("Reading eigen-decomposition with stem %s\n", eigenfile);

//first get indexes of individuals we want
sprintf(filename, "%s.eigen.id", eigenfile);
count=countrows(filename);
if(count!=ns)
{printf("Error, the decomposition is not suitable (%s contains %d samples, not %d)\nWhen using \"--decompose\" make sure you provide the same sample filterings as now (and if there are missing phenotypes, also use \"--pheno\")\n\n", filename, count, ns);exit(1);}
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(ids3, ns, wantids, count, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count2==0)
{printf("Error, the decomposition is not suitable (%s contains none of the %d samples)\nWhen using \"--decompose\" make sure you provide the same sample filterings as now (and if there are missing phenotypes, also use \"--pheno\")\n\n", filename, ns);exit(1);}
if(count2<ns)
{printf("Error, the decomposition is not suitable (%s contains only %d of the %d samples)\nWhen using \"--decompose\" make sure you provide the same sample filterings as now (and if there are missing phenotypes, also use \"--pheno\")\n\n", filename, count2, ns);exit(1);}

//open decomp and check size
sprintf(filename,"%s.eigen.bin",eigenfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*(count+1)*count)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*(count+1)*count, ftello(input));exit(1);}

datatemp=malloc(sizeof(float)*count);

if(axes==-9999)	//read and save all eigenvalues and (probably) eigenvectors
{
fseeko(input, 0, SEEK_SET);
if(fread(datatemp, sizeof(float), count, input)!=count)
{printf("Error reading eigenvalues from %s\n\n", filename);exit(1);}
for(i=0;i<count;i++){E[i]=datatemp[i];}

if(U!=NULL)	//now each eigenvector in turn
{
for(i=0;i<count;i++)
{
if(fread(datatemp, sizeof(float), count, input)!=count)
{printf("Error reading values for Eigenvector %d from %s\n\n", i+1, filename);exit(1);}
for(i2=0;i2<count;i2++){U[(size_t)i*count+indexer[i2]]=datatemp[i2];}
}
}
}
else	//read and save last axes eigenvalues
{
fseeko(input, (off_t)sizeof(float)*(count-axes), SEEK_SET);
if(fread(datatemp, sizeof(float), axes, input)!=axes)
{printf("Error reading eigenvalues from %s\n\n", filename);exit(1);}
for(i=0;i<axes;i++){E[i]=datatemp[i];}

fseeko(input, (off_t)sizeof(float)*(1+count-axes)*count, SEEK_SET);
for(i=0;i<axes;i++)
{
if(fread(datatemp, sizeof(float), count, input)!=count)
{printf("Error reading values for Eigenvector %d from %s\n\n", count-axes+i+1, filename);exit(1);}
for(i2=0;i2<count;i2++){U[(size_t)i*count+indexer[i2]]=datatemp[i2];}
}
}
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);};free(wantids);
free(indexer);free(datatemp);

printf("\n");
return(0);
}	//end of read_eigenfile

///////////////////////////

