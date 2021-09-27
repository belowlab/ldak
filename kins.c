/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Computing kinships

///////////////////////////

if(bitsize>data_length){bitsize=data_length;}

//allocate variables (and set kins to zero)

if(single==0)
{
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
}
else
{
data_warn2(bitsize,num_samples_use*3/2);
data=malloc(sizeof(double)*num_samples_use*bitsize);
data_single=malloc(sizeof(float)*num_samples_use*bitsize);
}

if(single==0)
{
anal_warn(num_samples_use, num_samples_use);
kins=calloc((size_t)num_samples_use*num_samples_use,sizeof(double));
}
else
{
anal_warn(num_samples_use, num_samples_use/2);
kins_single=calloc((size_t)num_samples_use*num_samples_use,sizeof(single));
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
if(mode==112){sprintf(filename,"%s/progress.%d",folder,partition);}
else{sprintf(filename,"%s.progress",outfile);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
count=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating kinships for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating kinships for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, power, strcmp(centresfile,"blank")!=0, hwestand, weights+bitstart, 1, preds+bitstart);

if(single==0)
{
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &bitlength, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);
}
else
{
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++){data_single[(size_t)j*num_samples_use+i]=(float)data[(size_t)j*num_samples_use+i];}
}

alpha_single=1.0;beta_single=1.0;
sgemm_("N", "T", &num_samples_use, &num_samples_use, &bitlength, &alpha_single, data_single, &num_samples_use, data_single, &num_samples_use, &beta_single, kins_single, &num_samples_use);
}
}	//end of bit loop
printf("\n");

count=0;for(j=0;j<data_length;j++){count+=(mults[j]!=-9999);}
if(count==0){printf("Error, all %d predictors trivial\n\n", data_length);exit(1);}

//save kins
if(mode==112){sprintf(outfile,"%skinships.%d", folder, partition);}
if(single==0){write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, preds, keeppreds_use, centres, mults, weights, al1, al2, data_length, datafile, power, kingz, kinraw, 1);}
else{write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 1, preds, keeppreds_use, centres, mults, weights, al1, al2, data_length, datafile, power, kingz, kinraw, 1);}

free(data);if(single==1){free(data_single);}
if(single==0){free(kins);}else{free(kins_single);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

