/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Joining correlations

///////////////////////////

//read first root file, which gets datafile, num_samples, num_preds, num_samples_use, cutoff, window_kb/window_cm
sprintf(filename, "%s.cors.root", corstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "Datafile %s ", datafile)!=1)
{printf("Error reading Row 1 of %s (should begin \"Datafile\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Samples %d ", &num_samples)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Samples\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors %d ", &num_preds)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Samples_Used %d ", &num_samples_use)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Predictors_Used %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}

if(fscanf(input, "Threshold %lf ", &mincor)!=1)
{printf("Error reading Row 7 of %s (should begin \"Threshold\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 8 of %s, suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(strcmp(readstring,"Window_kb")!=0&&strcmp(readstring,"Window_cM")!=0)
{printf("Error, Row 8 of %s should begin with \"Window_kb\" or \"Window_cM\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n", filename, readstring);exit(1);}
if(strcmp(readstring,"Window_kb")==0){window_kb=atof(readstring2);}
else{window_cm=atof(readstring2);}
fclose(input);

//check size of first bin file
sprintf(filename,"%s.cors.bin", corstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)(sizeof(int)*2+sizeof(double)*3)*num_preds+(sizeof(int)+sizeof(float))*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)(sizeof(int)*2+sizeof(double)*3)*num_preds+(sizeof(int)+sizeof(float))*scount, ftello(input));exit(1);}
fclose(input);

//now read remaining roots and bins, checking consistent and sizes
for(k=1;k<num_cors;k++)
{
sprintf(filename, "%s.cors.root", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "Datafile %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Datafile\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches datafile
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(datafile,readstring)!=0&&strcmp(datafile,readstring2)!=0)
{printf("Error, the correlations %s and %s appear to correspond different data files (%s and %s); if you are sure both are correct, use \"--check-root NO\"\n\n", corstems[0], corstems[k], datafile, readstring);exit(1);}
}

if(fscanf(input, "Num_Samples %d ", &readint)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Samples\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(readint!=num_samples)
{printf("Error, the correlations %s and %s correspond to data files with different numbers of samples (%d and %d)\n\n", corstems[0], corstems[k], num_samples, readint);exit(1);}

if(fscanf(input, "Num_Predictors %d ", &readint)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(readint!=num_preds)
{printf("Error, the correlations %s and %s correspond to data files with different numbers of predictors (%d and %d)\n\n", corstems[0], corstems[k], num_preds, readint);exit(1);}

if(fscanf(input, "Num_Samples_Used %d ", &readint)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(readint!=num_samples_use)
{printf("Error, the correlations %s and %s were made using different numbers of samples (%d and %d)\n\n", corstems[0], corstems[k], num_samples_use, readint);exit(1);}

if(fscanf(input, "Num_Predictors_Used %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}

if(fscanf(input, "Threshold %lf ", &readdouble)!=1)
{printf("Error reading Row 7 of %s (should begin \"Threshold\"), suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(readdouble!=mincor)
{printf("Error, the correlations %s and %s were constructed using different thresholds (%f and %f)\n\n", corstems[0], corstems[k], mincor, readdouble);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 8 of %s, suggesting the file has been changed since creation with \"--calc-cors\"\"\n\n",filename);exit(1);}
if(strcmp(readstring,"Window_kb")!=0&&strcmp(readstring,"Window_cM")!=0)
{printf("Error, Row 8 of %s should begin with \"Window_kb\" or \"Window_cM\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n", filename, readstring);exit(1);}

if(atof(readstring2)!=window_kb&&atof(readstring2)!=window_cm)
{printf("Error, the correlations %s and %s were constructed using different window sizes\n\n", corstems[0], corstems[k]);exit(1);}
fclose(input);

sprintf(filename,"%s.cors.bin", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)(sizeof(int)*2+sizeof(double)*3)*num_preds+(sizeof(int)+sizeof(float))*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)(sizeof(int)*2+sizeof(double)*3)*num_preds+(sizeof(int)+sizeof(float))*scount, ftello(input));exit(1);}
fclose(input);
}	//end of k loop

////////

//allocate variables
inputs=malloc(sizeof(FILE*)*num_cors);

usedpreds=malloc(sizeof(int)*num_preds);
actnums=malloc(sizeof(int)*num_preds);
centres=malloc(sizeof(double)*num_preds);
mults=malloc(sizeof(double)*num_preds);
sqdevs=malloc(sizeof(double)*num_preds);

writeints=malloc(sizeof(int)*num_preds);
writedoubles=malloc(sizeof(double)*num_preds);

bigs=malloc(sizeof(int*)*num_preds);
rjks=malloc(sizeof(float*)*num_preds);

//open all bin files
for(k=0;k<num_cors;k++)
{
sprintf(filename,"%s.cors.bin", corstems[k]);
if((inputs[k]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
}

//get predictors used for each correlation, counts, means, scalings and variances
for(j=0;j<num_preds;j++){usedpreds[j]=-9999;actnums[j]=0;centres[j]=-9999;mults[j]=-9999;sqdevs[j]=-9999;}

for(k=0;k<num_cors;k++)
{
//check predictors are new
fseeko(inputs[k], 0, SEEK_SET);
if(fread(writeints, sizeof(int), num_preds, inputs[k])!=num_preds)
{printf("Error reading predictor indicators from %s.cors.bin\n\n", corstems[k]);exit(1);}
for(j=0;j<num_preds;j++)
{
if(writeints[j]==1)
{
if(usedpreds[j]!=-9999){printf("Error, a predictor was used when making both %s and %s; each predictor can be used only to make one of the correlations\n\n", corstems[usedpreds[j]], corstems[k]);exit(1);}
usedpreds[j]=k;
}
}

//update counts
if(fread(writeints, sizeof(int), num_preds, inputs[k])!=num_preds)
{printf("Error reading predictor counts from %s.cors.bin\n\n", corstems[k]);exit(1);}
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==k){actnums[j]=writeints[j];}
}

//update means
if(fread(writedoubles, sizeof(double), num_preds, inputs[k])!=num_preds)
{printf("Error reading predictor means from %s.cors.bin\n\n", corstems[k]);exit(1);}
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==k){centres[j]=writedoubles[j];}
}

//update scalings
if(fread(writedoubles, sizeof(double), num_preds, inputs[k])!=num_preds)
{printf("Error reading predictor scalings from %s.cors.bin\n\n", corstems[k]);exit(1);}
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==k){mults[j]=writedoubles[j];}
}

//update means
if(fread(writedoubles, sizeof(double), num_preds, inputs[k])!=num_preds)
{printf("Error reading predictor variances from %s.cors.bin\n\n", corstems[k]);exit(1);}
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==k){sqdevs[j]=writedoubles[j];}
}
}

////////

//open output bin and save headers
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

for(j=0;j<num_preds;j++){writeints[j]=(usedpreds[j]!=-9999);}
fwrite(writeints, sizeof(int), num_preds, output2);
fwrite(actnums, sizeof(int), num_preds, output2);
fwrite(centres, sizeof(double), num_preds, output2);
fwrite(mults, sizeof(double), num_preds, output2);
fwrite(sqdevs, sizeof(double), num_preds, output2);

//now copy values from input bins to output bin
count=0;
for(j=0;j<num_preds;j++)
{
if(j%100000==0){printf("Processing Predictor %d out of %d\n", j+1, num_preds);}

if(usedpreds[j]!=-9999)
{
if(actnums[j]>0)
{
k=usedpreds[j];

bigs[j]=malloc(sizeof(int)*actnums[j]);
rjks[j]=malloc(sizeof(float)*actnums[j]);

if(fread(bigs[j], sizeof(int), actnums[j], inputs[k])!=actnums[j])
{printf("Error reading predictor indexes from %s.cors.bin\n\n", corstems[k]);exit(1);}
if(fread(rjks[j], sizeof(float), actnums[j], inputs[k])!=actnums[j])
{printf("Error reading predictor correlations from %s.cors.bin\n\n", corstems[k]);exit(1);}

fwrite(bigs[j], sizeof(int), actnums[j], output2);
fwrite(rjks[j], sizeof(float), actnums[j], output2);

free(bigs[j]);
free(rjks[j]);
}

count++;
}
}
printf("\n");

for(k=0;k<num_cors;k++){fclose(inputs[k]);}
fclose(output2);

//print root
scount=0;for(j=0;j<num_preds;j++){scount+=actnums[j];}

sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Num_Samples %d\nNum_Predictors %d\n", num_samples, num_preds);
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, count);
fprintf(output3,"Num_Pairs %jd\n", scount);
fprintf(output3,"Threshold %.4e\n", mincor);
if(window_cm!=-9999){fprintf(output3,"Window_cM %.4f\n", window_cm);}
else{fprintf(output3,"Window_kb %.4f\n", window_kb);}
fclose(output3);

printf("The joined correlations are saved in %s and %s\n\n", filename2, filename3);


free(inputs);
free(usedpreds);free(actnums);free(centres);free(mults);free(sqdevs);
free(writeints);free(writedoubles);
free(bigs);free(rjks);

num_samples=-9999;
num_preds=-9999;
num_samples_use=-9999;

///////////////////////////

