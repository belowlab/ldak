/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//do work for splitting summary statistics

///////////////////////////

if(bitsize>data_length){bitsize=data_length;}

//allocate variables
if(strcmp(proffile,"blank")==0)
{
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
}

randnorms=malloc(sizeof(double)*num_samples_use);
subrhos=malloc(sizeof(double)*data_length);
restrhos=malloc(sizeof(double)*data_length);
datarands=malloc(sizeof(double)*bitsize);

//set randnorms
for(i=0;i<num_samples_use;i++){randnorms[i]=RNOR;}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

////////

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating summaries for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating summaries for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

//read data and standardize
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1, preds+bitstart);

//divide predictor j by root(nss[j])
for(j=bitstart;j<bitend;j++)
{
value=pow(nss[j],-.5);
for(i=0;i<num_samples_use;i++){data[i+(j-bitstart)*num_samples_use]*=value;}
}

//get t(data) * randnorms * root(1/num_samples_use*(1-subprop)/subprop)
alpha=pow(1.0/num_samples_use*(1-subprop)/subprop,.5);beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, randnorms, &one, &beta, datarands, &one);

//subrhos is rhos plus datarands
for(j=bitstart;j<bitend;j++){subrhos[j]=rhos[j]+datarands[j-bitstart];}

//restrhos is complement
for(j=bitstart;j<bitend;j++){restrhos[j]=(rhos[j]-subprop*subrhos[j])/(1-subprop);}
}	//end of bit loop
printf("\n");

////////

//save summaries
sprintf(filename2,"%s.train.summaries",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor A1 A2 Direction Stat n\n");

sprintf(filename3,"%s.test.summaries",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Direction Stat n\n");

for(j=0;j<data_length;j++)
{
fprintf(output2, "%s %c %c %d %.4f %.1f\n", preds[j], al1[j], al2[j], (subrhos[j]>=0)-(subrhos[j]<0), subprop*nss[j]*pow(subrhos[j],2)/(1-pow(subrhos[j],2)), subprop*nss[j]);
fprintf(output3, "%s %c %c %d %.4f %.1f\n", preds[j], al1[j], al2[j], (restrhos[j]>=0)-(restrhos[j]<0), (1-subprop)*nss[j]*pow(restrhos[j],2)/(1-pow(restrhos[j],2)), (1-subprop)*nss[j]);
}

fclose(output2);
fclose(output3);

printf("New summary statistics saved in %s and %s\n\n", filename2, filename3);

free(data);
free(randnorms);free(subrhos);free(restrhos);free(datarands);

///////////////////////////

