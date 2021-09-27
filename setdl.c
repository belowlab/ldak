/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//set data_length and keeppreds - will only be here if use_data=1

///////////////////////////

//set data_length to num_preds_use (number after filtering), then overwrite (or update) when necessary
data_length=num_preds_use;
keeppreds_use=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds_use;j++){keeppreds_use[j]=keeppreds[j];}

////////

//mode=101 - cut-weights
//nothin=0 - data_length starts as number after filtering (will later be reduced by thinning)
//nothin=1 - data_length is number after filtering
//nothin=2 - data_length is number after thinning

//mode=102/103/104 - calc-weights, join-weights, calc-weights-direct
//nothin=0 - data_length is number after thinning
//nothin=1 - data_length is number after filtering

if((mode==101&&nothin==2)||((mode==102||mode==103||mode==104)&&nothin==0))
{
//read in thinned data, and check all predictors present in data (before and after any predictor filterings)
sprintf(filename,"%sthin.in",folder);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--cut-weights\"\n\n", filename);exit(1);}

data_length=countrows(filename);
printf("Reading list of %d thinned predictors from %s\n", data_length, filename);
wantpreds=malloc(sizeof(char*)*data_length);
read_strings(filename, wantpreds, data_length, NULL, 1, 0);

count=find_strings(allpreds, num_preds, wantpreds, data_length, keeppreds_use, NULL, NULL, filename, predorder, NULL, 3);
if(count<data_length){printf("Error, only %d of these are in the data\n\n", count);exit(1);}

usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds[j]]++;}

count=0;for(j=0;j<data_length;j++){count+=(usedpreds[keeppreds_use[j]]==1);}
if(count<data_length)	//must have failed due to extract
{printf("Error, after filtering predictors, only %d remain, indicating %s was created using different predictor filterings\n\n", count, filename);exit(1);}
for(j=0;j<data_length;j++){free(wantpreds[j]);}free(wantpreds);free(usedpreds);
printf("\n");
}

//mode=105 - adjust-weights - data_length is number after filtering

//mode=106 - thin - data_length is number after filtering
//mode=107 - thin-tops - data_length is number after filtering (already filtered based on pvalues)

if(mode==108)	//find-tags
{
//read in scorefile - keeppreds_use will first be set to give indexes of tags
blupcentres=malloc(sizeof(double*)*num_scores);
blupfactors=malloc(sizeof(double*)*num_scores);
num_seek=read_scores(blupcentres, blupfactors, keeppreds_use,num_scores, scorefile, num_preds, allpreds, allal1, allal2, predorder, bimfile, -9999, NULL, 1);

//reset datalength and keeppreds_use to include both filtered predictors and those in scorefile
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds[j]]=1;}
for(j=0;j<num_seek;j++){usedpreds[keeppreds_use[j]]+=2;}

//will be using usedpreds 1, 2 and 3
data_length=0;for(j=0;j<num_preds;j++){data_length+=(usedpreds[j]>0);}
sindex=malloc(sizeof(int)*num_seek);
stypes=malloc(sizeof(int)*data_length);
count=0;count2=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]>0)	//will be using
{
keeppreds_use[count]=j;stypes[count]=usedpreds[j];
if(usedpreds[j]>1)	//in scorefile
{sindex[count2]=count;count2++;}
count++;
}
}

free(usedpreds);
}

if(mode==109)	//remove-tags
{
//read in targetfile - keeppreds_use will first be set to give indexes of target predictors
count=countrows(targetfile);
printf("Reading list of %d target predictors from %s\n", count, targetfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(targetfile, wantpreds, count, NULL, 1, 0);

num_seek=find_strings(allpreds, num_preds, wantpreds, count, keeppreds_use, NULL, NULL, targetfile, predorder, NULL, 3);
if(num_seek==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(num_seek<count){printf("Warning, only %d of these are in the data\n\n", num_seek);}
else{printf("All of these are in the data\n\n");}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);

//reset datalength and keeppreds_use to include both filtered predictors and those in targetfile
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds[j]]=1;}
for(j=0;j<num_seek;j++){usedpreds[keeppreds_use[j]]+=2;}

//will be using usedpreds 1, 2 and 3
data_length=0;for(j=0;j<num_preds;j++){data_length+=(usedpreds[j]>0);}
sindex=malloc(sizeof(int)*num_seek);
stypes=malloc(sizeof(int)*data_length);
count=0;count2=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]>0)	//will be using
{
keeppreds_use[count]=j;stypes[count]=usedpreds[j];
if(usedpreds[j]>1)	//in targetfile
{sindex[count2]=count;count2++;}
count++;
}
}

free(usedpreds);
}

////////

//mode=111 - cut-kins - data_length is number after filtering

if(mode==112)	//calc-kins - data_length is number predictors in partition
{
if(pstarts[partition-1]!=-9999)
{
data_length=pends[partition-1]-pstarts[partition-1];
for(j=0;j<data_length;j++){keeppreds_use[j]=keeppreds[pstarts[partition-1]+j];}
}
else	//read from file (allow for predictor filterings)
{
//keeppreds_use will first be set to give indexes of partition predictors
sprintf(filename, "%s%d", partpref, partition);
count=countrows(filename);
printf("Reading %d predictors from %s\n", count, filename);
wantpreds=malloc(sizeof(char*)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);

count2=find_strings(allpreds, num_preds, wantpreds, count, keeppreds_use, NULL, NULL, NULL, predorder, NULL, 3);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the data\n", count2);}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);

//reset datalength and keeppreds_use to allow for any predictor filtering
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds[j]]=1;}
for(j=0;j<count2;j++){usedpreds[keeppreds_use[j]]++;}

data_length=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==2){keeppreds_use[data_length]=j;data_length++;}
}
if(data_length==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
if(data_length<count2){printf("Warning, after filtering predictors, only %d remain\n", data_length);}
printf("\n");

free(usedpreds);
}
}

//mode=113 - join-kins - use_data!=1
//mode=114 - calc-kins-direct - data_length is number after filtering

//mode=115 - filter - use_data!=1
//mode=116 - add-grm - use_data!=1
//mode=117 - sub-grm - if(extract=1), data_length is number after filtering, otherwise use_data!=1
//mode=118 - convert-gz - use_data!=1
//mode=119 - convert-raw - use_data!=1
//mode=120 - calc-sim-grm - use_data!=1

///////

//mode=121 - reml - use_data!=1

if(mode==122)	//calc-blups - get effects and rands (will not have predictor filterings)
{
blupcentres=malloc(sizeof(double*)*(num_kins+num_scores+(num_tops>0)));
blupfactors=malloc(sizeof(double*)*(num_kins+num_scores+(num_tops>0)));
data_length=read_effects(blupcentres, blupfactors, keeppreds_use, num_kins, kinstems, num_scores+(num_tops>0), regfile, num_preds, allpreds, allal1, allal2, predorder, bimfile);

if(num_kins>0)	//get random effects
{
bluprands=malloc(sizeof(double)*num_samples_use*num_kins);
read_rand(blupfile, bluprands, num_kins, num_samples_use, ids3);
}
printf("\n");
}

//mode=123 - he - use_data!=1
//mode=124 - pcgc - use_data!=1

//mode=126 - fast-he - data_length is number after filtering
//mode=127 - fast-pcgc - data_length is number after filtering

////////

//mode=131 - linear - data_length is number after filtering
//mode=132 - logistic - data_length is number after filtering
//mode=133 - solve-null - use_data!=1

//mode=136 - cut-genes - data_length is number after filtering
//mode=137 - calc-genes-kins - data_length is number after filtering
//mode=138 - calc-genes-reml - data_length is number after filtering
//mode=139 - join-genes-reml - if(cut1!=-9999), data_length is number after filtering, otherwise use_data!=1

////////

//mode=141 - calc-tagging - data_length is number after filtering
//mode=142 - join-tagging - use_data!=1
//mode=143 - merge-tagging - use_data!=1
//mode=144 - reduce-tagging - use_data!=1

//mode=146 - sum-hers - use_data!=1
//mode=147 - sum-cors - use_data!=1
//mode=149 - calc-exps - use_data!=1
//mode=150 - calc-posts - use_data!=1

////////

//mode=151 - ridge - data_length is number after filtering
//mode=152 - bolt-predict - data_length is number after filtering
//mode=153 - bayesr - data_length is number after filtering

//mode=156 - calc-cors - data_length is number after filtering
//mode=157 - join-cors - use_data!=1
//mode=158 - mega-prs - data_length is number after filtering
//mode=159 - pseudo-summaries - data_length is number after filtering

////////

//mode=161 - pca - use_data!=1

if(mode==162)	//calc-pca-loads - get effects, vectors and values (will not have predictor filterings)
{
blupcentres=malloc(sizeof(double*));
blupfactors=malloc(sizeof(double*));
data_length=read_effects(blupcentres, blupfactors, keeppreds_use, 1, kinstems, 0, NULL, num_preds, allpreds, allal1, allal2, predorder, bimfile);

bluprands=malloc(sizeof(double*)*num_samples_use*axes);
blupvalues=malloc(sizeof(double)*axes);
read_pca(pcastem, bluprands, blupvalues, axes, num_samples_use, ids3);
}

//mode=163 - decompose - use_data!=1
//mode=164 - adjust-grm - use_data!=1

//mode=166 - truncate -grm- use_data!=1
//mode=167 - pca-grm - use_data!=1
//mode=168 - square-grm - use_data!=1
//mode=169 - gxemm-iid - use_data!=1
//mode=170 - gxemm-free - use_data!=1

////////

//mode=171 - calc-stats - data_length is number after filtering

if(mode==172)	//calc-scores - get effects (might have predictor filterings)
{
blupcentres=malloc(sizeof(double*)*num_scores);
blupfactors=malloc(sizeof(double*)*num_scores);
data_length=read_scores(blupcentres, blupfactors, keeppreds_use, num_scores, scorefile, num_preds, allpreds, allal1, allal2, predorder, bimfile, num_preds_use, keeppreds, 0);
}

if(mode==173)	//make-phenos (might have predictor filterings)
{
effects=malloc(sizeof(double*)*num_phenos);
data_length=get_effects(effects, bivar, keeppreds_use, num_phenos, num_causals, causalsfile, effectsfile, num_preds_use, keeppreds, num_preds, allpreds, predorder);
}

//mode=174 - make-snps - use_data!=1

//mode=176 - jackknife - use_data!=1
//mode=177 - cut-folds - use_data!=1
//mode=178 - find-gaussian - use_data!=1

////////

//mode=181/182/183/184/185 - making data - use_data!=1
//mode=186/187/188/189 - condensing data - use_data!=1
//mode=190 - calc-sim-data - use_data!=1

///////////////////////////

//now can reduce annotations and load cmbp (will contain either 1000000xcm or bp)

chr=malloc(sizeof(int)*data_length);
cm=malloc(sizeof(double)*data_length);
bp=malloc(sizeof(double)*data_length);
cmbp=malloc(sizeof(double)*data_length);
preds=malloc(sizeof(char*)*data_length);
al1=malloc(sizeof(char)*data_length);
al2=malloc(sizeof(char)*data_length);
for(j=0;j<data_length;j++)
{
chr[j]=allchr[keeppreds_use[j]];
cm[j]=allcm[keeppreds_use[j]];
bp[j]=allbp[keeppreds_use[j]];
if(window_cm!=-9999){cmbp[j]=1000000*allcm[keeppreds_use[j]];}
else{cmbp[j]=allbp[keeppreds_use[j]];}
copy_string(preds,j,allpreds[keeppreds_use[j]]);
al1[j]=allal1[keeppreds_use[j]];
al2[j]=allal2[keeppreds_use[j]];
}

///////////////////////////

