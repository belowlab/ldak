//I previously allowed for projecting the data, to speed up single-snp lmm
//Here are the related scripts

//some declarations in declare.c
char projfile2[500]="blank", projfile[500]="";

//command line stuff from readargs.c
if(strcmp(argv[count],"--project")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=144;strcpy(outfile2,argv[count+1]);found=1;
}
if(strcmp(argv[count],"--projection")==0)
{strcpy(projfile2,argv[count+1]);found=1;}

//requirements from required.c
if(mode==144)	//project
{
if(strcmp(eigenfile,"blank")==0)
{printf("Error, you must use \"--eigen\" to provide an eigen-decomposition (onto which to project the data)\n\n");exit(1);} 
}

//told in param.c
if(mode==144)
{
printf("Projecting a dataset onto an eigen-decomposition of a kinship matrix\n\n");
}

//appended in append.c
(void)append_check(projfile,projfile2,workdir);

//some things in getnums.c

if(mode==144)	//reduce to samples in eigen-decomposition
{
sprintf(filename,"%s.eigen.id", eigenfile);
count=countrows(filename);
printf("Reading list of %d samples to use from %s\n", count, filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, idsorder, NULL, 3);
if(count2==0){printf("Error, can not find any of these\n\n");exit(1);}
if(count2<count){printf("Error, can only find %d of these\n", count2);exit(1);}
for(i=0;i<count2;i++){usedids[indexer[i]]++;}
for(i=0;i<num_samples;i++){usedids[i]=(usedids[i]==2);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

if(mode==144&&strcmp(eigenfile,"blank")!=0)	//check samples align with projection
{
sprintf(filename,"%s.project.id", eigenfile);
count=countrows(filename);
if(count!=num_samples_use)
{printf("Error, the number of samples being used does not match the number in %s (%d), indicating different sample filterings were used when making the eigen-decomposition\n\n", filename, count);exit(1);}
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

for(i=0;i<num_samples_use;i++)
{
if(strcmp(ids3[i],wantids[i])!=0)
{printf("Error, ID %d in %s should be %s (not %s), indicating different sample filterings were used when making the eigen-decomposition\n\n", i+1, filename, ids3[i], wantids[i]);exit(1);}
}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}


if(mode==151&&strcmp(projfile,"blank")!=0)	//check projfile predictors and size
{
sprintf(filename,"%s.project.predictors", projfile);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=num_preds_use)
{printf("Error, %s contains %d predictors, indicating different sample filterings are being used\n\n", filename, count);exit(1);} 

printf("Reading list of %d predictors from %s\n", count, filename);
wantpreds=malloc(sizeof(char*)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);
indexer=malloc(sizeof(int)*count);
count2=find_strings(wantpreds, count, allpreds, num_preds, NULL, indexer, NULL, NULL, NULL, predorder, 3);
if(count2<count)
{printf("Error, only %d of these are in the data, indicating a different dataset was used with \"--project\"\n\n", count2);exit(1);}
 
count2=0;
for(j=0;j<count;j++){count2+=(usedpreds[indexer[j]]==1);}
if(count2<count)	//this covers the unlikely case that the count=num_preds_used but not same predictors
{printf("Error, after filtering samples, only %d of these remain, indicating different sample filterings have been used\n\n", count2);exit(1);}

for(j=0;j<num_preds_use;j++)
{
if(strcmp(allpreds[keeppreds[j]],wantpreds[j])!=0){printf("Error 811E, please tell Doug %d of %d %s %s\n", j+1, num_preds_use, allpreds[keeppreds[j]], wantpreds[j]);exit(1);}
}

sprintf(filename,"%s.project.bin",projfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*count*num_samples_use)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*count*num_samples_use, ftello(input));exit(1);}
fclose(input);

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}

////////

//some checks in consistent.c

if((strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0)&&mode==144)
{printf("Error, you can not use \"--keep\" or \"--remove\" with \"--project\" (the samples will be read from the eigen-decomposition)\n\n");exit(1);}

if(strcmp(projfile,"blank")!=0)
{printf("Error, you must use \"--grm\" to provide the kinship matrix to which the projection corresponds\n\n");exit(1);}
if(strcmp(projfile,"blank")!=0&&mode!=151)
{printf("Error, you can only use \"--projection\" with \"--linear\"\n\n");exit(1);}

if(checkroot!=-9999&&mode!=132&&mode!=133&&mode!=134&&mode!=142&&strcmp(tracefile,"blank")==0&&strcmp(eigenfile,"blank")==0&&strcmp(projfile,"blank")==0)
{printf("Error, you can only use \"--check-root\" with \"--calc-blups\", \"--he\", \"--pcgc\" or \"--calc-pca-loads\", or when using \"--traces\", \"--eigen\" or \"--projection\"\n\n");exit(1);}

//and some more in parsefiles.c

if(strcmp(projfile,"blank")!=0)
{
if(strcmp(eigenfile,"blank")==0)
{printf("Error, if using \"--project\", you must also use \"--eigen\" to provide the eigen-decomposition used to make the projection\n\n");exit(1);}

sprintf(filename,"%s.project.id", projfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--project\"\n\n", filename);exit(1);}
sprintf(filename,"%s.project.predictors", projfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--project\"\n\n", filename);exit(1);}
sprintf(filename,"%s.project.bin", projfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--project\"\n\n", filename);exit(1);}

sprintf(filename,"%s.project.root", projfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--project\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Datafile %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Datafile\"), suggesting the file has been changed since creation with \"--decompose\"\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches datafile and decomposition
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(datastems[0],readstring)!=0&&strcmp(datastems[0],readstring2)!=0)
{printf("Error, the projection %s corresponds to the predictor file %s, which appears to be different to that provided now (%s); if you are sure the projection is correct, use \"--check-root NO\"\n\n", projfile, readstring, datastems[0]);exit(1);}
if(fscanf(input, "Decomposition %s ", readstring)!=1)
{printf("Error reading Row 2 of %s (should begin \"Decomposition\"), suggesting the file has been changed since creation with \"--decompose\"\"\n\n",filename);exit(1);}
(void)append_check(readstring2,readstring,workdir);
if(strcmp(eigenfile,readstring)!=0&&strcmp(eigenfile,readstring2)!=0)
{printf("Error, the projection %s corresponds to the eigen-decomposition %s, which appears to be different to that provided now (%s); if you are sure the projection is correct, use \"--check-root NO\"\n\n", projfile, readstring, kinstems[0]);exit(1);}
}
fclose(input);
}

////////

//mode=144 created the projection

if(mode==144)	//project
{
if(bitsize>data_length){bitsize=data_length;}
data_warn2(bitsize,2*num_samples_use);
eigen_warn(num_samples_use);

data=malloc(sizeof(double)*num_samples_use*bitsize);

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);
UTdata=malloc(sizeof(double)*num_samples_use*bitsize);

read_eigens(eigenfile, U, E, num_samples_use, ids3, -9999);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.project.bin",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Projecting for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Projecting for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data, centres+bitstart, mults+bitstart, NULL, num_samples_use, bitlength, missingvalue, 0, 0, 0, NULL, 1, preds+bitstart);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_samples_use, &bitlength, &num_samples_use, &alpha, U, &num_samples_use, data, &num_samples_use, &beta, UTdata, &num_samples_use);

for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++)
{
writefloat=UTdata[(size_t)j*num_samples_use+i];
if(fwrite(&writefloat, sizeof(float), 1, output2)!=1)
{printf("Error writing projected values for Predictor %s to %s\n\n", preds[j], filename2);exit(1);}
}}
}	//end of bit loop
printf("\n");
fclose(output2);

sprintf(filename3,"%s.project.id",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
for(i=0;i<num_samples_use;i++){fprintf(output3,"%s %s\n", ids1[i], ids2[i]);}
fclose(output3);

sprintf(filename4,"%s.project.predictors",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
for(j=0;j<num_preds_use;j++){fprintf(output4,"%s\n", preds[j]);}
fclose(output4);

sprintf(filename5,"%s.project.root",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename5);exit(1);}
fprintf(output5,"Datafile %s\nKinship %s", datafile, kinstems[0]);
fclose(output5);

printf("Projection saved in %s, %s, %s and %s", filename2, filename3, filename4, filename5);

free(data);
if(binary==0){gzclose(datainputgz);}
free(U);free(E);free(UTdata);
}	//end of mode=144

///////////////////////////

