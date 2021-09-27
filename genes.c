/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Gene-based analyses

///////////////////////////

//count number of gene in this partition and set genemax to size of longest gene
count=0;genemax=0;
for(g=0;g<num_genes;g++)
{
if(gparts[g]==partition)
{
count++;
if(gends[g]-gstarts[g]>genemax){genemax=gends[g]-gstarts[g];}
}
}
printf("The longest gene/chunk contains %d predictors\n\n", genemax);

//allocate variables
data_warn(num_samples_use, genemax);
data=malloc(sizeof(double)*num_samples_use*genemax);

if(mode==137)
{
anal_warn(num_samples_use, num_samples_use);
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
}
else
{
if(num_kins==0){anal_warn(num_samples_use+2*genemax, genemax);}
else{anal_warn(5*num_samples_use+3*genemax, genemax);}

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
X=malloc(sizeof(double)*num_samples_use*genemax);
if(strcmp(sumsfile,"blank")!=0)
{
Xnss=malloc(sizeof(double)*genemax);
Xrhos=malloc(sizeof(double)*genemax);
}

retain=malloc(sizeof(int)*genemax);
stats=malloc(sizeof(double)*8);
if(gene_perms>0){permlike=malloc(sizeof(double)*gene_perms);}
if(num_kins==1){vstarts=malloc(sizeof(double)*3);}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

if(mode==138)	//fill Y and Z, and if using, get factor - also stuff for nk=1
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

if(prev!=-9999){factor=get_factor(Y, num_samples_use, prev, ascer, NULL);}

if(num_kins==1)	//get starting heritabilities, varnull and likenull (saved in vstarts)
{
printf("Solving Null Model\n");
(void)adv_reml(num_samples_use, num_fixed, 0, Y, Z, U, E, kintraces, NULL, -9999, -9999, NULL, vstarts, tol, maxiter);
printf("\n");
}
}

//deal with progress and on-the-fly files
sprintf(filename,"%sprogress.%d", folder, partition);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

if(mode==138)
{
sprintf(filename2,"%sremls.%d", folder, partition);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Gene_Number Gene_Name Num_Predictors Heritability SD Null_Likelihood Alt_Likelihood LRT_Stat LRT_P_Raw Her_Liability SD_Liability");
for(p=0;p<gene_perms;p++){fprintf(output2, " Perm_%d", p+1);}
fprintf(output2, "\n");
}

////////

count2=0;
for(g=0;g<num_genes;g++)
{
total=gends[g]-gstarts[g];

if(gparts[g]==partition)
{
if(mode==137)
{
if(count2%200==0)
{
printf("Calculating kinships for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating kinships for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
}
}
else
{
if((num_kins==0&&count2%200==0)||(num_kins==1&&count2%20==0))
{
printf("Performing REML for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Performing REML for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
}
}

shuffle=0;	//this is how many predictors we already have data for (can be more than total)
for(j=gstarts[g];j<end;j++)
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(j-start)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, gstarts[g]+shuffle, gends[g], keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp);
stand_data(data+(size_t)shuffle*num_samples_use, centres+gstarts[g]+shuffle, mults+gstarts[g]+shuffle, sqdevs+gstarts[g]+shuffle, num_samples_use, total-shuffle, missingvalue, power, 0, hwestand, weights+gstarts[g]+shuffle, 1, preds+gstarts[g]+shuffle);

//ready to do work
if(mode==137)	//calc-kins-genes
{
//find number of non-trivial predictors
count3=0;for(j=gstarts[g];j<gends[g];j++){count3+=(mults[j]!=-9999&&weights[j]>0);}

if(count3==0)	//gene trivial
{
printf("Warning, all %d predictors in Gene/Chunk %d are trivial", gends[g]-gstarts[g], g+1);
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n");
}
else	//non-trivial
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &total, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);

sprintf(outfile,"%sgeneships.%d", folder, g+1);
write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, preds+gstarts[g], keeppreds_use+gstarts[g], centres+gstarts[g], mults+gstarts[g], weights+gstarts[g], al1+gstarts[g], al2+gstarts[g], total, datafile, power, kingz, kinraw, 0);
}
}
else	//calc-kins-reml
{
//find number of non-trivial predictors (indexed by retain)
if(strcmp(sumsfile,"blank")==0)
{count3=prune_gene(retain, gprune, data, num_samples_use, total, mults+gstarts[g], weights+gstarts[g], NULL, X);}
else
{count3=prune_gene(retain, gprune, data, num_samples_use, total, mults+gstarts[g], weights+gstarts[g], nss+gstarts[g], X);}

if(count3==0)	//gene trivial
{
printf("Warning, all %d predictors in Gene/Chunk %d are trivial", gends[g]-gstarts[g], g+1);
if(strcmp(sumsfile,"blank")==0)
{
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n");
}
else
{
if(strcmp(weightsfile,"blank")!=0){printf(", have weight zero");}
printf(" or are missing summary statistics\n");
}
fprintf(output2, "%d %s 0 NA NA NA NA NA NA NA NA", g+1, gnames[g]);
for(p=0;p<gene_perms;p++){fprintf(output2, " NA");}
fprintf(output2, "\n");
}
else	//non-trivial
{
//analyse real data
if(strcmp(sumsfile,"blank")==0)
{
sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, NULL, retain, data, NULL, NULL);
if(num_kins==0)
{(void)gene_reml(num_samples_use, num_fixed, Y, Z, X, count3, sumsq, NULL, NULL, stats, tol, maxiter, shrink, strip);}
else
{(void)adv_reml(num_samples_use, num_fixed, 1, Y, Z, U, E, kintraces, X, count3, sumsq, stats, vstarts, tol, maxiter);}
}
else
{
sumsq=fill_gene(X, Xnss, Xrhos, num_samples_use, total, NULL, retain, data, nss+gstarts[g], rhos+gstarts[g]);
(void)gene_reml(num_samples_use, num_fixed, Y, Z, X, count3, sumsq, Xnss, Xrhos, stats, tol, maxiter, shrink, strip);
}

for(p=0;p<gene_perms;p++)
{
permute_int(order,num_samples_use);

if(strcmp(sumsfile,"blank")==0)
{
sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, order, retain, data, NULL, NULL);
if(num_kins==0)
{permlike[p]=gene_reml(num_samples_use, num_fixed, Y, Z, X, count3, sumsq, NULL, NULL, NULL, tol, maxiter, shrink, strip);}
else
{permlike[p]=adv_reml(num_samples_use, num_fixed, 1, Y, Z, U, E, kintraces, X, count3, sumsq, NULL, vstarts, tol, maxiter);}
}
else
{
sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, order, retain, data, nss+gstarts[g], rhos+gstarts[g]);
permlike[p]=gene_reml(num_samples_use, num_fixed, Y, Z, X, count3, sumsq, Xnss, Xrhos, NULL, tol, maxiter, shrink, strip);
}
}

//print results
fprintf(output2, "%d %s %d ", g+1, gnames[g], count3);
if(stats[1]!=-9999){fprintf(output2, "%.6f %.6f ", stats[0], stats[1]);}
else{fprintf(output2, "%.6f NA ", stats[0]);}
fprintf(output2, "%.4f %.4f %.4f %.4e ", stats[2], stats[3],  stats[4], stats[5]);
if(prev!=-9999)
{
if(stats[1]!=-9999){fprintf(output2,"%.6f %.6f ", stats[0]*factor, stats[1]*factor);}
else{fprintf(output2, "%.6f NA", stats[0]*factor);}
}
else{fprintf(output2, "NA NA");}
for(p=0;p<gene_perms;p++){fprintf(output2, " %.4f", permlike[p]);}
fprintf(output2, "\n");
}
}	//end of mode=138

count2++;
start=gstarts[g];if(gends[g]>end){end=gends[g];}
}}	//end of in right partition and g loop
printf("\n");

fclose(output);
if(mode==138){fclose(output2);}

if(mode==137){printf("Kinship matrices saved with stem %sgeneships\n\n", folder);}
else{printf("REML estimates saved in %sremls.%d\n\n", folder, partition);}

free(data);
if(mode==137){free(kins);}
else
{
free(order);free(Y);free(Z);free(X);
if(strcmp(sumsfile,"blank")!=0){free(Xnss);free(Xrhos);}
free(retain);free(stats);
if(gene_perms>0){free(permlike);}
if(num_kins==1){free(vstarts);}
}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

