/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculating correlations, saving those above a threshold - ignore self correlations

///////////////////////////

if(mincor==-9999)	//find threshold corresponding to P=0.01
{
mincor=1-exp(-6.634897/num_samples_use);
printf("Will record pairs of predictors with correlation squared at least %.4e (this corresponds to P<0.01)\n", mincor);
}

if(bitsize==-9999)	//get optimal bitsize (average number of neighbours) - always using window_kb
{
scount=0;
k=1;
for(j=0;j<data_length;j++)
{
while(1)	//move along until out of range, tallying each data_length predictors encountered
{
if(k==data_length){break;}
if(chr[k]!=chr[j]||cmbp[k]-cmbp[j]>1000*window_kb){break;}
k++;
}
scount+=k-j-1;
}
bitsize=(int)((scount-1)/data_length)+1;

if(bitsize<20){bitsize=20;}
if(bitsize>8000){bitsize=8000;}
printf("The bit-size will be set to %d (you can change this using \"--bit-size\")\n\n", bitsize);
}
else
{
if(bitsize>data_length){bitsize=data_length;}
}

step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//work out bitmax
bitmax=bitsize;
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
if(bitend2-bitstart>bitmax){bitmax=bitend2-bitstart;}
}

//work out max number of correlations for each predictor
maxnums=malloc(sizeof(int)*data_length);
for(j=0;j<data_length;j++){maxnums[j]=0;}

//first go to left
k=data_length-2;
for(j=data_length-1;j>=0;j--)
{
while(1)
{
if(k==-1){break;}
if(chr[k]!=chr[j]||cmbp[j]-cmbp[k]>1000*window_kb){break;}
k--;
}
maxnums[j]+=j-k-1;
}

//now right
k=1;
for(j=0;j<data_length;j++)
{
while(1)
{
if(k==data_length){break;}
if(chr[k]!=chr[j]||cmbp[k]-cmbp[j]>1000*window_kb){break;}
k++;
}
maxnums[j]+=k-j-1;
}

//work out max allocations within any bit (as multiple of bitmax)
bitmax=bitsize;
smax=0;
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
if(bitend2-bitstart>bitmax){bitmax=bitend2-bitstart;}

scount=0;for(j=bitstart;j<bitend2;j++){scount+=maxnums[j];}
if(scount/bitmax>smax){smax=scount/bitmax;}
}

////////

//allocate variables
data_warn3(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);

anal_warn(bitmax, bitmax+smax);
cors=malloc(sizeof(double)*bitmax*bitmax);
bigs=malloc(sizeof(int*)*data_length);
rjks=malloc(sizeof(float*)*data_length);

actnums=malloc(sizeof(int)*data_length);

writeints=malloc(sizeof(int)*num_preds);
writedoubles=malloc(sizeof(double)*num_preds);

//set actnums to zero
for(j=0;j<data_length;j++){actnums[j]=0;}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//open bin file
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//write indicators of predictors used
for(j=0;j<num_preds;j++){writeints[j]=0;}
for(j=0;j<data_length;j++){writeints[keeppreds_use[j]]=1;}
fwrite(writeints, sizeof(int), num_preds, output2);

//write placeholders for counts, means, scalings and variances
for(j=0;j<num_preds;j++){writeints[j]=0;writedoubles[j]=0;}
fwrite(writeints, sizeof(int), num_preds, output2);
fwrite(writedoubles, sizeof(double), num_preds, output2);
fwrite(writedoubles, sizeof(double), num_preds, output2);
fwrite(writedoubles, sizeof(double), num_preds, output2);

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
bitlength=bitend2-bitstart;

if(bit%step==0)
{
printf("Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

shuffle=0;
for(j=0;j<end-bitstart;j++)	//using values already in data, so shuffle back
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(bitstart-start+j)*num_samples_use+i];}
shuffle++;
}

//allocate for predictors not in previous chunk
for(j=bitstart+shuffle;j<bitend2;j++)
{
if(maxnums[j]>0)
{
bigs[j]=malloc(sizeof(int)*maxnums[j]);
rjks[j]=malloc(sizeof(float)*maxnums[j]);
}
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, bitstart+shuffle, bitend2, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data+(size_t)shuffle*num_samples_use, centres+bitstart+shuffle, mults+bitstart+shuffle, sqdevs+bitstart+shuffle,  num_samples_use, bitlength-shuffle, missingvalue, -1, 0, 0, NULL, 1, preds+bitstart+shuffle);

//get correlation
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors, &bitlength);

//now loop through predictors in chunk (have already set actnums to zero)
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)	//non trivial
{
for(k=j+1;k<bitend2;k++)
{
if(chr[k]!=chr[j]){break;}
if(cmbp[k]-cmbp[j]>1000*window_kb){break;}

if(mults[k]!=-9999)	//non-trivial and within range
{
value=cors[(size_t)(k-bitstart)*bitlength+(j-bitstart)];

if(pow(value,2)>=mincor)	//above threshold
{
bigs[j][actnums[j]]=keeppreds_use[k];
rjks[j][actnums[j]]=value;
actnums[j]++;

bigs[k][actnums[k]]=keeppreds_use[j];
rjks[k][actnums[k]]=value;
actnums[k]++;
}
}}	//end of using k and k loop

if(actnums[j]>0)	//print correlations for this predictor
{
fwrite(bigs[j], sizeof(int), actnums[j], output2);
fwrite(rjks[j], sizeof(float), actnums[j], output2);
}
}}	//end of using j and j loop

//free predictors that will not be used again
for(j=bitstart;j<bitend;j++)
{
if(maxnums[j]>0){free(bigs[j]);free(rjks[j]);}
}

start=bitstart;if(bitend2>end){end=bitend2;}
}	//end of bit loop
printf("\n");

fclose(output2);

//re-open cors file
if((output2=fopen(filename2,"rb+"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fseeko(output2, sizeof(int)*num_preds, SEEK_SET);

//update counts
for(j=0;j<num_preds;j++){writeints[j]=0;}
for(j=0;j<data_length;j++){writeints[keeppreds_use[j]]=actnums[j];}
fwrite(writeints, sizeof(int), num_preds, output2);

//update means
for(j=0;j<num_preds;j++){writedoubles[j]=-9999;}
for(j=0;j<data_length;j++){writedoubles[keeppreds_use[j]]=centres[j];}
fwrite(writedoubles, sizeof(double), num_preds, output2);

//update scalings
for(j=0;j<num_preds;j++){writedoubles[j]=-9999;}
for(j=0;j<data_length;j++){writedoubles[keeppreds_use[j]]=mults[j];}
fwrite(writedoubles, sizeof(double), num_preds, output2);

//update variances
for(j=0;j<num_preds;j++){writedoubles[j]=-9999;}
for(j=0;j<data_length;j++){writedoubles[keeppreds_use[j]]=sqdevs[j];}
fwrite(writedoubles, sizeof(double), num_preds, output2);

fseeko(output2, 0, SEEK_END);
fclose(output2);

//print root
scount=0;for(j=0;j<data_length;j++){scount+=actnums[j];}

sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Num_Samples %d\nNum_Predictors %d\n", num_samples, num_preds);
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, data_length);
fprintf(output3,"Num_Pairs %jd\n", scount);
fprintf(output3,"Threshold %.4e\n", mincor);
if(window_cm!=-9999){fprintf(output3,"Window_cM %.4f\n", window_cm);}
else{fprintf(output3,"Window_kb %.4f\n", window_kb);}
fclose(output3);

printf("For each predictor, there are on average %.2f other predictors with correlation squared at least %.6f\n\n", (double)scount/data_length, mincor);

printf("The correlations are saved in %s and %s\n\n", filename2, filename3);

free(maxnums);
free(data);
free(cors);free(bigs);free(rjks);
free(actnums);
free(writeints);free(writedoubles);
if(binary==0){gzclose(datainputgz);}

///////////////////////////
