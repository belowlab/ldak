/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Checks and defaults specific to different modes - plus some general ones at bottom

///////////////////////////

//calculating weights, thinning and finding/removing tags

if(mode==101)	//cut-weights
{
if(nothin==-9999){nothin=0;}
if(wprune==-9999){wprune=0.98;}
if(window_cm!=-9999){window_kb=window_cm*1000;}
if(window_kb==-9999&&window_length==-9999){window_kb=100;}

if(section_cm!=-9999){section_kb=section_cm*1000;}
if(section_kb==-9999&&section_length==-9999){section_length=1000;}
if(buffer_cm!=-9999){buffer_kb=buffer_cm*1000;}
if(buffer_kb==-9999&&buffer_length==-9999){buffer_kb=window_kb;buffer_length=window_length;}
}

if(mode==102)	//calc-weights
{
if(window_cm!=-9999){window_kb=window_cm*1000;}

if(section==-9999)
{printf("Error, you must use \"--section\" to specify which of the %d section to analyse, or use \"--calc-weights-all\" to analyse each section in turn\n\n", num_sections);exit(1);}
section_start=section;

if(lddecay==-9999){lddecay=0;}
if(lddecay==1)
{
if(halflife==-9999)
{
if(window_cm!=-9999){printf("Error, when using \"--decay YES\" you must also use \"--half-life\" to specify the assumed rate of LD decay in cM (to match the units used with \"--cut-weights\"); we suggest 1cM\n\n");}
else{printf("Error, when using \"--decay YES\" you must also use \"--half-life\" to specify the assumed rate of LD decay in kb (to match the units used with \"--cut-weights\"); we suggest 1000kb\n\n");}
exit(1);
}
decay=log(2)/halflife;
}
else{decay=0;}

if(maxtime==-9999){maxtime=360;}
if(fudge==-9999){fudge=0;}
if(simplex==-9999){simplex=0;}

if(mincor==-9999){mincor=0;}
}

if(mode==103)	//join-weights
{
if(spread==-9999){spread=0;}
}

if(mode==104)	//calc-weights-all
{
if(window_cm!=-9999){window_kb=window_cm*1000;}

if(section_start==-9999){section_start=1;}

if(lddecay==-9999){lddecay=0;}
if(lddecay==1)
{
if(halflife==-9999)
{
if(window_cm!=-9999){printf("Error, when using \"--decay YES\" you must also use \"--half-life\" to specify the assumed rate of LD decay in cM (to match the units used with \"--cut-weights\"); we suggest 1cM\n\n");}
else{printf("Error, when using \"--decay YES\" you must also use \"--half-life\" to specify the assumed rate of LD decay in kb (to match the units used with \"--cut-weights\"); we suggest 1000kb\n\n");}
exit(1);
}
decay=log(2)/halflife;
}
else{decay=0;}

if(maxtime==-9999){maxtime=360;}
if(fudge==-9999){fudge=0;}
if(simplex==-9999){simplex=0;}
if(spread==-9999){spread=0;}

if(mincor==-9999){mincor=0;}
}

if(mode==105)	//adjust-weights
{
if(strcmp(partpref,"blank")==0)
{printf("Error, you should use \"--partition-number\" and \"--partition-prefix\" to specify which predictors are in each partition\n\n");exit(1);}

if(reduce==-9999){reduce=1;}
}

if(mode==106)	//thin
{
if(wprune==-9999)
{printf("Error, you must use \"--window-prune\" to specify the correlation squared threshold\n\n");exit(1);}

if(window_kb==-9999&&window_length==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\", \"--window-length\" or \"--window-cm\" to specify the window size in kb, as a fixed number of predictors, or in centiMorgans\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}
}

if(mode==107)	//thin-tops
{
if(strcmp(pvafile,"blank")==0||cutoff==-9999)
{printf("Error, you must use both \"--pvalues\" and \"--cutoff\" (predictors with p-values below the cutoff will be treated as top predictors)\n\n");exit(1);}

if(wprune==-9999)
{printf("Error, you must use \"--window-prune\" to specify the correlation squared threshold\n\n");exit(1);}

if(window_kb==-9999&&window_length==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\", \"--window-length\" or \"--window-cm\" to specify the window size in kb, as a fixed number of predictors, or in centiMorgans\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}
}

if(mode==108)	//find-tags
{
if(strcmp(scorefile,"blank")==0)
{printf("Error, you must use \"--scorefile\" to provide a file with at least five columns: \"Predictor\", \"A1\", \"A2\" and \"Centre\" (the average allele count, set to NA if unknown), followed by a column of effect sizes for each profile\n\n");exit(1);}

if(extract==0)
{printf("Error, you must use \"--extract\" and/or \"--exclude\" to specify which predictors can be used as tags\n\n");exit(1);}

if(window_kb==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\" or \"--window-cm\" to specify the window size in kb or centiMorgans\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}

if(mincor==-9999){mincor=0;}
}

if(mode==109)	//remove-tags
{
if(strcmp(targetfile,"blank")==0)
{printf("Error, you must use \"--targets\" to specify the predictors for which you seek tags\n\n");exit(1);}

if(window_kb==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\" or \"--window-cm\" to specify the window size in kb or centiMorgans\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}

if(mincor==-9999){printf("Error, you must use \"--min-cor\" to specify a minimum correlation squared threshold\n\n");exit(1);}
if(mincor==0){printf("Error, \"--min-cor\", must be followed by a float greater than 0\n\n");exit(1);}
}

///////////////////////////

//calculating and manipulating kinships

if(mode==111)	//cut-kins
{
if(part_length==-9999&&bychr!=1&&strcmp(partpref,"blank")==0)
{printf("Error, you should use \"--partition-number\" and \"--partition-prefix\" to specify which predictors are in each partition, else use \"--partition-length\" or \"--by-chr YES\" to divide predictors evenly or by chromosome (note that if your aim is to compute kinships across all predictors, it is generally easier to use \"--calc-kins-direct\")\n\n");exit(1);}

if(bychr==-9999){bychr=0;}
if(checkpart==-9999){checkpart=1;}
}

if(mode==112)	//calculating kinships for one partition
{
if(partition==-9999&&num_parts==1){partition=1;}
if(partition==-9999){printf("Error, you must use \"--partition\" to specify which of the %d partitions to analyse\n\n", num_parts);exit(1);}

if(single==-9999){single=0;}
}

//mode=113 - nothing to do

if(mode==114)	//calc-kins-direct
{
if(single==-9999){single=0;}
}

if(mode==115)	//filter
{
if(kinstand==-9999){kinstand=1;}
}

//mode=116/117/118/119 - nothing to do

//mode=120 - nothing to do

///////////////////////////

//reml, blup and he/pcgc (shortcut used indirectly by other functions)

if(mode==121)	//reml
{
if(shortcut==-9999)
{
if(num_kins<2){shortcut=1;}
else{shortcut=0;}
}

if(compact==-9999){compact=0;}

if(liab==-9999){liab=0;}
if(liab==1&&prev==-9999)
{printf("Error, you must use \"--prevalence\" to specify the proportion of samples who are cases in the population\n\n");exit(1);}

if(constrain==-9999){constrain=0;}

if(permute==-9999){permute=0;}
if(shrink==-9999){shrink=1;}
if(strip==-9999){strip=0;}
}

if(mode==122)	//calc-blups
{
if(strcmp(remlfile,"blank")==0)
{printf("Error, you must use \"--remlfile\" to provide results from \"--reml\"\n\n");exit(1);}
}

if(mode==123)	//he
{
if(num_blocks==-9999){num_blocks=-1;}
if(permute==-9999){permute=0;}
}

if(mode==124)	//pcgc
{
if(prev==-9999){printf("Error, you must use \"--prevalence\" to specify the proportion of samples who are cases in the population\n\n");exit(1);}

if(num_blocks==-9999){num_blocks=-1;}
if(permute==-9999){permute=0;}
}

if(mode==126)	//fast-he
{
if(strcmp(partpref,"blank")==0&&strcmp(annpref,"blank")==0)	//just base
{parttype=0;num_parts=0;}
else
{
if(strcmp(partpref,"blank")!=0)	//using partitions, so num_parts already set
{parttype=1;}
if(strcmp(annpref,"blank")!=0)	//using annotations - easier to copy num_anns and annpref into num_parts and partpref
{parttype=0;num_parts=num_anns;strcpy(partpref,annpref);}
}

if(num_vects==-9999){printf("Error, you must use \"--repetitions\" to specify the number of random vectors (we suggest 100)\n\n");exit(1);}

if(num_blocks==-9999){num_blocks=200;}
}

if(mode==127)	//fast-pcgc
{
if(prev==-9999){printf("Error, you must use \"--prevalence\" to specify the proportion of samples who are cases in the population\n\n");exit(1);}

if(strcmp(partpref,"blank")==0&&strcmp(annpref,"blank")==0)	//just base
{parttype=0;num_parts=0;}
else
{
if(strcmp(partpref,"blank")!=0)	//using partitions, so num_parts already set
{parttype=1;}
if(strcmp(annpref,"blank")!=0)	//using annotations - easier to copy num_anns and annpref into num_parts and partpref
{parttype=0;num_parts=num_anns;strcpy(partpref,annpref);}
}

if(num_vects==-9999){printf("Error, you must use \"--repetitions\" to specify the number of random vectors (we suggest 100)\n\n");exit(1);}

if(num_blocks==-9999){num_blocks=200;}
}

///////////////////////////

//association analysis

if(mode==131)	//linear
{
if(exact==-9999){exact=1;}

if(mincor==-9999){mincor=0.8;}
if(constrain==-9999){constrain=1;}
if(permute==-9999){permute=0;}
}

if(mode==132)	//logistic
{
if(mincor==-9999){mincor=0.8;}
if(permute==-9999){permute=0;}
}

if(mode==133)	//solve-null
{
if(constrain==-9999){constrain=1;}
}

////////

if(mode==136)	//cut-genes
{
if(bychr==-9999){bychr=0;}

if(strcmp(genefile,"blank")==0&&chunks==-9999&&chunksbp==-9999)
{printf("Error, you should either use \"--genefile\" to divide predictors based on annotations, else use \"--chunks\" or \"--chunks-bp\" to divide into chunks containing a specified number of predictors or basepairs; if dividing based on annotations, the first four columns of the genefile should specify the name, chromosome, start and end basepairs of each gene (using 0-start, half-open positions); if the genefile has a fifth column, this should specify the orientation (+ or -)\n\n");exit(1);}

if(gene_buffer==-9999){gene_buffer=0;}
if(up_buffer==-9999){up_buffer=gene_buffer;down_buffer=gene_buffer;}
if(minweight==-9999){minweight=1e-10;}
if(overlap==-9999){overlap=1;}
}

if(mode==137)	//calc-genes-kins
{
if(partition==-9999&&num_parts==1){partition=1;}
if(partition==-9999){printf("Error, you must use \"--partition\" to specify which of the %d partitions to analyse\n\n", num_parts);exit(1);}
}

if(mode==138)	//calc-gene-reml
{
if(partition==-9999&&num_parts==1){partition=1;}
if(partition==-9999){printf("Error, you must use \"--partition\" to specify which of the %d partitions to analyse\n\n", num_parts);exit(1);}

if(gene_perms==-9999){gene_perms=1;}

if(gprune==-9999)
{
if(strcmp(weightsfile,"blank")!=0)	//with weights, no pruning should be required
{gprune=2;}
else	//without, pruning depends on whether using summaries or not
{
if(strcmp(sumsfile,"blank")!=0){gprune=0.5;}
else{gprune=0.98;}
}
}

if(permute==-9999){permute=0;}
if(shrink==-9999){shrink=1;}
if(strip==-9999){strip=0;}
}

//mode=139 - nothing to do

///////////////////////////

//sumher

if(mode==141)	//calc-tagging
{
if(window_kb==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\" or \"--window-cm\" to specify the window size in which to search for tags; we recommend 1cM (or 1000kb if genetic distances are not available)\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}

if(strcmp(partpref,"blank")==0&&strcmp(annpref,"blank")==0)	//just base
{parttype=0;num_parts=0;}
else
{
if(strcmp(partpref,"blank")!=0)	//using partitions, so num_parts already set
{parttype=1;}
if(strcmp(annpref,"blank")!=0)	//using annotations - easier to copy num_anns and annpref into num_parts and partpref
{parttype=0;num_parts=num_anns;strcpy(partpref,annpref);}
}

if(reduce==-9999){reduce=0;}

if(unbias==-9999){unbias=1;}
if(savemat==-9999){savemat=0;}
if(cover==-9999){cover=0;}

if(mincor==-9999){mincor=0;}
}

if(mode==142)	//join tagging
{
if(strcmp(taglist,"blank")==0)
{printf("Error, you must use \"--taglist\" to provide a list of tagging files\n\n");exit(1);}

if(checkdups==-9999){checkdups=1;}
}

if(mode==143)	//merge tagging
{
if(strcmp(taglist,"blank")==0)
{printf("Error, you must use \"--taglist\" to provide a list of tagging files\n\n");exit(1);}
}

if(mode==144)	//reduce tagging
{
if(strcmp(tagfile,"blank")==0)
{printf("Error, you must use \"--tagfile\" to provide the tagging file created by \"--calc-tagging\"\n\n");exit(1);}

if(strcmp(catfile,"blank")==0)
{printf("Error, you must use \"--categories\" to specify to which partitions/annotations to retain\n\n");exit(1);}
}

if(mode==146)	//sum-hers
{
if(strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--summary\" to provide the results from single-predictor analysis\n\n");exit(1);}

if(strcmp(tagfile,"blank")==0)
{printf("Error, you must use \"--tagfile\" to provide the tagging file created by \"--calc-tagging\"\n\n");exit(1);}

if(strcmp(taufile,"blank")!=0)
{
count=countrows(taufile);
count2=countcols(taufile);
if(count!=num_parts+gcon+cept||count2!=1)
{printf("Error, %s should have %d rows and 1 column (not %d and %d); you can construct this file from the \".taus\" file produced by \"--sum-hers\" by extracting the second column, then appending the scaling factor (if using \"--genomic-control YES\") followed by the intercept (if using \"--intercept YES\")\n\n", taufile, num_parts+gcon+cept, count, count2);exit(1);}
}

if(checksums==-9999){checksums=1;}
if(gcon==-9999){gcon=0;}
if(cept==-9999){cept=0;}

if(ldsc==-9999){ldsc=0;}

if(chisol==-9999)
{
if(ldsc==1){chisol=0;}
else{chisol=1;}
}

if(tagone==-9999){tagone=0;}

if(divide!=-9999)
{
if(parttype==0||num_parts==1)
{printf("Error, you can only use \"--divisions\" if the tagging file contains multiple partitions\n\n");exit(1);}
if(num_parts%divide!=0)
{printf("Error, the argument to \"--divisions\" (%d) is not a factor of the number of partitions (%d)\n\n", divide, num_parts);exit(1);}

if(uptaus==-9999){uptaus=0;}
if(strcmp(powfile,"blank")!=0)
{
num_pows=countrows(powfile);
if(num_pows!=num_parts/divide)
{printf("Error, %s should have %d rows, providing the power for each analyses (not %d)\n\n", powfile, num_parts/divide, num_pows);exit(1);}
}
}

if(num_blocks==-9999)
{
if(ldsc==1||chisol==0){num_blocks=200;}
}
}

if(mode==147)	//sum-cors
{
if(strcmp(sumsfile,"blank")==0||strcmp(sums2file,"blank")==0)
{printf("Error, you must use \"--summary\" and \"--summary2\" to provide two sets of results from single-predictor analysis\n\n");exit(1);}

if(strcmp(tagfile,"blank")==0)
{printf("Error, you must use \"--tagfile\" to provide the tagging file created by \"--calc-tagging\"\n\n");exit(1);}

if(checksums==-9999){checksums=1;}

if(gcon==-9999)
{
if(cept==1){gcon=0;}
else{gcon=1;}
}
if(cept==-9999){cept=0;}

if(tagone==-9999){tagone=0;}

if(plet==-9999){plet=0;}

if(num_blocks==-9999){num_blocks=200;}
}

if(mode==148)	//calculating individual heritabilities
{
if(strcmp(taufile,"blank")==0)
{printf("Error, you must use \"--taus\" to provide the results from \"--sum-hers\"\n\n");exit(1);}

if(strcmp(matfile,"blank")==0)
{printf("Error, you must use \"--matrix\" to provide the heritability matrix created by \"--calc-tagging\" with \"--save-matrix YES\"\n\n");exit(1);}

count=countrows(taufile);
count2=countcols(taufile);
if(count!=1+num_parts||count2!=2+num_parts)
{printf("Error, %s should have %d rows and %d columns (not %d and %d)\n\n", taufile, 1+num_parts, 2+num_parts, count, count2);exit(1);}
}

if(mode==149)	//calculating expectations
{
if(strcmp(tagfile,"blank")==0)
{printf("Error, you must use \"--tagfile\" to provide the tagging file created by \"--calc-tagging\"\n\n");exit(1);}

if(strcmp(taufile,"blank")==0&&her==-9999)
{printf("Error, you must either use \"--taus\", to provide the results from \"--sum-hers\", or use \"--her\" and \"--props\", to specify the total heritabilty explained by all predictors and the proportion of this explained by each category\n\n");exit(1);}

if(strcmp(taufile,"blank")!=0)
{
count=countrows(taufile);
count2=countcols(taufile);
if(count!=1+num_parts||count2!=2+num_parts)
{printf("Error, %s should have %d rows and %d columns (not %d and %d)\n\n", taufile, 1+num_parts, 2+num_parts, count, count2);exit(1);}
}
else	//using her
{
if(num_parts>1&&strcmp(propfile,"blank")==0)
{printf("Error, you must use \"--props\" to provide the proportion of heritability explained by each category in %s\n\n", tagfile);exit(1);}

if(strcmp(propfile,"blank")!=0)
{
count=countrows(propfile);
count2=countcols(propfile);
if(count!=num_parts||count2!=1){printf("Error, %s should have %d rows and 1 column (not %d and %d), providing the proportion of heritability explained by each category in %s\"\n\n", propfile, num_parts, count, count2, tagfile);exit(1);}
}
}
}

if(mode==150)	//calculating bfs
{
if(strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--summary\" to provide the results from single-predictor analysis\n\n");exit(1);}

if(strcmp(expfile,"blank")==0)
{printf("Error, you must use \"--expectations\" to provide the results from \"--calc-exps\"\n\n");exit(1);}
}

///////////////////////////

//individual-level data prediction, then megaprs

if(mode==151)	//ridge
{
if(strcmp(indhers,"blank")==0&&her==-9999)	//can only be here if have provided weights and power, but not her
{printf("Error, you must use \"--her\" to specify the heritability contributed by all predictors\n\n");exit(1);}

if(restart==-9999){restart=0;}
if(skipcv==-9999){skipcv=1;}
}

if(mode==152)	//bolt
{
if(strcmp(indhers,"blank")==0&&her==-9999)	//can only be here if have provided weights and power, but not her
{printf("Error, you must use \"--her\" to specify the heritability contributed by all predictors\n\n");exit(1);}

if(restart==-9999){restart=0;}

if(skipcv==-9999){skipcv=0;}
if(cvprop==-9999&&strcmp(bvsfile,"blank")==0&&skipcv==0)
{printf("Error, you must use \"--cv-proportion\" to specify the fraction of samples used for testing the different parameters (we suggest 0.1); alternatively use \"--cv-samples\" to explicitly specify the test samples, or \"--skip-cv YES\" to turn off the cross-validation\n\n");exit(1);}
}

if(mode==153)	//bayesr
{
if(strcmp(indhers,"blank")==0&&her==-9999)	//can only be here if have provided weights and power, but not her
{printf("Error, you must use \"--her\" to specify the heritability contributed by all predictors\n\n");exit(1);}

if(restart==-9999){restart=0;}

if(skipcv==-9999){skipcv=0;}
if(cvprop==-9999&&strcmp(bvsfile,"blank")==0&&skipcv==0)
{printf("Error, you must use \"--cv-proportion\" to specify the fraction of samples used for testing the different parameters (we suggest 0.1); alternatively use \"--cv-samples\" to explicitly specify the test samples, or \"--skip-cv YES\" to turn off the cross-validation\n\n");exit(1);}

if(pointmass==-9999){pointmass=1;}
if(fullspace==-9999){fullspace=0;}
}

////////

if(mode==156)	//calc-cors
{
if(window_kb==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\" or \"--window-cm\" to specify the window size in kb or centiMorgans; we recommend using \"--window-cm 3\" (or if genetic distances are not available, \"--window-kb 3000\")\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}
}

if(mode==157)	//join-cors
{
if(strcmp(corslist,"blank")==0)
{printf("Error, you must use \"--corslist\" to provide a list of correlations\n\n");exit(1);}
}

if(mode==158)	//mega-prs
{
if(strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--summary\" to provide the results from single-predictor analysis\n\n");exit(1);}

if(strcmp(sums2file,"blank")==0&&onesums!=1)
{printf("Error, you must use \"--summary2\" to provide training summary statistics (or use \"--one-sums YES\")\n\n");exit(1);}

if(ptype==-9999)
{printf("Error, you must use \"--model\" followed by \"lasso\", \"ridge\", \"bolt\", \"bayesr\", \"bayesr-shrink\", \"lasso-sparse\" or \"--mega\" to specify the model type (we generally recommend \"--model bayesr\")\n\n");exit(1);}

if(window_kb==-9999&&window_cm==-9999)
{printf("Error, you must use \"--window-kb\" or \"--window-cm\" to specify the window size in kb or centiMorgans; we recommend using \"--window-cm 1\" (or if genetic distances are not available, \"--window-kb 1000\")\n\n");exit(1);}
if(window_cm!=-9999){window_kb=window_cm*1000;}

if(strcmp(indhers,"blank")==0&&her==-9999)	//can only be here if have provided weights and power, but not her
{
if(ptype==1){her=0.5;}	//her not used for lasso-sparse
else{printf("Error, you must use \"--her\" to specify the heritability contributed by all predictors\n\n");exit(1);}
}

if(strcmp(corsfile,"blank")==0)
{printf("Error, you must use \"--cors\" to specify correlations between predictors in the dataset\n\n");exit(1);}

if(onesums==-9999){onesums=0;}
if(segments==-9999){segments=8;}
if(checksize==-9999){checksize=1;}

if(multiher==-9999){multiher=1;}
if(basicbolt==-9999){basicbolt=0;}
if(ldpred==-9999){ldpred=0;}

if(shrink==-9999){shrink=1;}
}

if(mode==159)	//pseudo-summaries
{
if(strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--summary\" to provide the results from single-predictor analysis\n\n");exit(1);}

if(subprop==-9999)
{printf("Error, you must use \"--training-proportion\" to specify the fraction of samples assumed to be training (we suggest 0.9)\n\n");exit(1);}
}

///////////////////////////

//pca, decompose and adjust-grm

if(mode==161)	//pca
{
if(axes==-9999)
{printf("Error, you must use \"--axes\" to specify the number of principal compoments (20 is usually sufficient)\n\n");exit(1);}
}

if(mode==162)	//calc-pca-loads
{
if(strcmp(pcastem,"blank")==0)
{printf("Error, you must use \"--pcastem\" to provide the results from \"--pca\"\n\n");exit(1);}
}

if(mode==163)	//decompose
{
if(eigenraw==-9999){eigenraw=0;}
}

if(mode==164)	//adjust-grm
{
if(strcmp(covarfile,"blank")==0&&strcmp(envfile,"blank")==0&&strcmp(topfile,"blank")==0)
{printf("Error, you must use \"--covar\", \"--enviro\" or \"--top-preds\" to provide fixed effects (on to which the kinship matrix will be regressed)\n\n");exit(1);}
}

if(mode==166)	//truncate-grm
{
if(cutoff==-9999)
{printf("Error, you must use \"--cutoff\" to specify the truncation threshold\n\n");exit(1);}
}

if(mode==167)	//pca-grm
{
if(axes==-9999)
{printf("Error, you must use \"--axes\" to specify how many principal components to use\n\n");exit(1);}
}

//mode=168 - nothing to do

if(mode==169)	//gxemm-iid
{
if(strcmp(envfile,"blank")==0)
{printf("Error, you must use \"--enviro\" to provide environmental variables\n\n");exit(1);}
}

if(mode==170)	//gxemm-free
{
if(strcmp(envfile,"blank")==0)
{printf("Error, you must use \"--enviro\" to provide environmental variables\n\n");exit(1);}
}

///////////////////////////

//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

//mode=171 - nothing to do

if(mode==172)	//calculating scores
{
if((strcmp(covarfile,"blank")!=0||strcmp(envfile,"blank")!=0)&&strcmp(cofile,"blank")==0)
{printf("Error, if using \"--covar\" or \"--enviro\" to include covariates, you must also use \"--coeffsfile\"\n\n");exit(1);}

if(strcmp(scorefile,"blank")==0)
{printf("Error, you must use \"--scorefile\" to provide a file with at least five columns: \"Predictor\", \"A1\", \"A2\" and \"Centre\" (the average allele count, set to NA if unknown), followed by a column of effect sizes for each profile\n\n");exit(1);}

if(savecounts==-9999){savecounts=0;}
}

if(mode==173)	//making phenotypes
{
if(num_phenos==-9999&&strcmp(causalsfile,"blank")!=0)
{num_phenos=countrows(causalsfile);}
if(num_phenos==-9999&&strcmp(effectsfile,"blank")!=0)
{num_phenos=countrows(effectsfile);}
if(num_phenos==-9999)
{printf("Error, you must use \"--num-phenos\" to specify how many phenotypes to generate\n\n");exit(1);}

if(num_causals==-9999&&strcmp(causalsfile,"blank")!=0)
{num_causals=countcols(causalsfile);}
if(num_causals==-9999&&strcmp(effectsfile,"blank")!=0)
{num_causals=countcols(effectsfile);}
if(num_causals==-9999&&her!=0)
{printf("Error, you must use \"--num-causals\" to specify how many predictors are causal for each phenotype (use \"--num-causals -1\" for the infinitessimal model where every predictor contributes)\n\n");exit(1);}

if(strcmp(causalsfile,"blank")!=0)	//check size
{
count=countrows(causalsfile);
count2=countcols(causalsfile);
if(count!=num_phenos)
{printf("Error, %s should have %d rows (not %d)\n\n", causalsfile, num_phenos, count);exit(1);}
if(count2!=num_causals)
{printf("Error, %s should have %d columns (not %d)\n\n", causalsfile, num_causals, count2);exit(1);}
}
if(strcmp(effectsfile,"blank")!=0)	//check size
{
count=countrows(effectsfile);
count2=countcols(effectsfile);
if(count!=num_phenos)
{printf("Error, %s should have %d rows (not %d)\n\n", effectsfile, num_phenos, count);exit(1);}
if(count2!=num_causals)
{printf("Error, %s should have %d columns (not %d)\n\n", effectsfile, num_causals, count2);exit(1);}
}

if(her==-9999)
{printf("Error, you must use \"--her\" to specify the heritability of each phenotype\n\n");exit(1);}

if(her==0)	//useful to pretend there is one causal
{
if(num_causals!=-9999||strcmp(causalsfile,"blank")!=0||strcmp(effectsfile,"blank")!=0)
{printf("Warning, \"--num_causals\", \"--causals\" and \"--effects\" are redundant with \"--her 0\"\n\n");}
num_causals=1;strcpy(causalsfile,"blank");strcpy(effectsfile,"blank");
}

if(bivar!=-9999&&num_phenos%2==1)
{printf("Error, you can only use \"--bivar\" when \"--num-phenos\" is even (not %d)\n\n", num_phenos);exit(1);}
}

if(mode==174)	//making snps
{
if(num_inds==-9999||num_snps==-9999)
{printf("Error, you must use \"--num-samples\" and \"--num-snps\" to specify how many samples and SNPs to generate\n");exit(1);}

if(maf1==-9999){maf1=0;}
if(maf2==-9999){maf2=0.5;}
}

if(mode==176)	//jackknifing
{
if(strcmp(jackfile,"blank")==0&&strcmp(proffile,"blank")==0)
{printf("Error, you must use either \"--data-pairs\", to provide a file containing predicted and observed values (optionally, this file can have a third column providing regression weights), or \"--profile\", to provide results from \"--calc-scores\"\n\n");exit(1);}

if(strcmp(jackfile,"blank")!=0)
{
if(countcols(jackfile)!=2&&countcols(jackfile)!=3){printf("Error, %s should have either two or three columns (not %d); the first two should provide predicted and observed values, the third (if used) should provide regression weights\n\n", jackfile, countcols(jackfile));exit(1);}
}
else
{
if(countcols(proffile)<6){printf("Error, %s should have at least six columns (not %d), suggesting the file has been changed since creation with \"--calc-scores\"\n\n\n\n", proffile, countcols(proffile));exit(1);}
if(countcols(proffile)%2==1){printf("Error, %s should have an even number of columns (not %d), suggesting the file has been changed since creation with \"--calc-scores\"\n\n\n\n", proffile, countcols(proffile));exit(1);}
}

if(auc==-9999){auc=0;}

if(num_blocks==-9999){printf("Error, you must use \"--num-blocks\" to specify the number of resamplings (we typically use 200)\n\n");exit(1);}
}

if(mode==177)	//cutting folds
{
if(dtype==-9999&&num_kins==0&&strcmp(respfile,"blank")==0)
{printf("Error, there are no IDs; you should provide a set of datafiles, one or more kinships, or phenotypes (if more than one of these is provided, will use the samples common to all)\n\n");exit(1);} 

if(num_folds==-9999)
{printf("Error, you must use \"--num-folds\" to specify the number of cross-validation folds\n\n\n");exit(1);}
}

if(mode==178)	//find gaussian
{
if(strcmp(likefile,"blank")==0)
{printf("Error, you must use \"--likelihoods\" to specify a list of values and corresponding likelihoods\n\n");exit(1);}

count=countcols(likefile);
if(count!=2&&count!=3){printf("Error, %s should have either 2 or 3 columns (not %d)\n\n", likefile, count);exit(1);}
}

///////////////////////////

//making and condensing data

if(mode==181||mode==182||mode==183||mode==184||mode==185)	//do all makings together
{
if(mode==181&&dtype!=1&&genprobs!=0&&threshold==-9999&&minprob==-9999)
{
if(genprobs==1)
{printf("Error, when using \"--make-bed\", you must also use \"--threshold\", followed by a float in [0.5,1], to specify how dosages are converted to hard genotypes\n\n");exit(1);}
printf("Error, when using \"--make-bed\", you must also use \"--threshold\", followed by a float in [0.5,1], or \"--minprob\", followed by 0 or a float in [0.5,1), to specify how genotype probabilities are converted to hard calls\n\n");exit(1);
}

if(usenames==-9999){usenames=0;}
if(comsamps==-9999){comsamps=0;}

if(compreds==-9999)
{
if(usenames==0){compreds=0;}
else{compreds=1;}
}

if(exsame==-9999){exsame=0;}
if(exdups==-9999){exdups=0;}

if(passall==-9999&&(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999)){passall=0;}

if(mode==184&&speedlong==-9999){speedlong=0;}
}

////////

if(mode==186||mode==187||mode==188||mode==189)	//do all condensings together
{
if(strcmp(genefile,"blank")==0&&chunks==-9999&&chunksbp==-9999)
{printf("Error, you should either use \"--genefile\" to divide predictors based on annotations, else use \"--chunks\" or \"--chunks-bp\" to divide into chunks containing a specified number of predictors or basepairs; if dividing based on annotations, the first four columns of the genefile should specify the name, chromosome, start and end basepairs of each gene (using 0-start, half-open positions); if the genefile has a fifth column, this should specify the orientation (+ or -)\n\n");exit(1);}

if(gene_buffer==-9999){gene_buffer=0;}
if(up_buffer==-9999){up_buffer=gene_buffer;down_buffer=gene_buffer;}
if(minweight==-9999){minweight=1e-10;}
if(overlap==-9999){overlap=1;}

if(nonsnp==0&&maxmaf==-9999)
{printf("Error, you must use \"--max-maf\" to specify the maximum MAF allowed\n\n");exit(1);}

if(mode==189&&speedlong==-9999){speedlong=0;}

if(useminor==-9999)
{
if(nonsnp==0){useminor=1;}
else{useminor=0;}
}
}

if(mode==190)	//calculating correlations between data
{
comsamps=1;
compreds=1;

if(exsame==-9999){exsame=0;}
if(exdups==-9999){exdups=0;}
}

///////////////////////////

//bitsize used for 101, 104, 106, 107, 112, 114, 117d, 122, 126, 127, 131, 132, 141, 151, 152, 153, 156, 159, 162, 171, 172, 173, 174, 181-190
if(bitsize==-9999)
{
//for modes 101, 104, 106, 107, 141, 156 will set intelligently (later on)
if(mode==112||mode==114||(mode==117&&extract==1)||mode==122||mode==126||mode==126||mode==127||mode==129||mode==131||mode==132||mode==159||mode==162||mode==171||mode==172||mode==173||mode==174||mode==186||mode==187||mode==188||mode==189||mode==190){bitsize=5000;}
if(mode==151||mode==152||mode==153){bitsize=64;}
if(mode==181||mode==182||mode==183||mode==184||mode==185){bitsize=2000;}
}

//tol used by modes 121, 124, 129, 131, 132, 133, 138, 146, 147, 151, 152, 153, 158, 161, 167
if(tol==-9999)	
{
if(mode==121||mode==124||mode==129||mode==131||mode==132||mode==133||mode==138||mode==146)	//used for likelihood
{tol=0.001;}
if(mode==147)	//sum-cors - used for heritabilities
{tol=0.0001;}
if(mode==151||mode==152||mode==153)	//ridge, bolt and bayesr
{tol=0.01;}
if(mode==158)	//mega-prs
{tol=0.00001;}
if(mode==161||mode==167)	//pca - dsyevx will set value
{tol=-1;}
}

//maxiter used by modes 102, 104, 121, 124, 129, 131, 132, 133, 138, 146, 147, 151, 152, 153, 158
if(maxiter==-9999)
{
if(mode==102||mode==104){maxiter=200000;}
if(mode==121||mode==138)
{
if(num_kins==0){maxiter=100;}
else{maxiter=25;}
}
if(mode==124||mode==129){maxiter=100;}
if(mode==131||mode==132||mode==133){maxiter=25;}
if(mode==146||mode==147){maxiter=100;}
if(mode==151||mode==152||mode==153){maxiter=150;}
if(mode==158){maxiter=50;}
}

//memsave used directly by modes 120, 138, 121, 123, 124, 131 and 133
if(memsave==-9999)	
{
if(mode==120){memsave=0;}	//not allowed memsave when computing correlations
if(mode==121)
{
if(num_kins<=1){memsave=1;}	//if using a kinship, just need the decomposition
else{memsave=0;}
}
if(mode==123||mode==124){memsave=0;}	//will offer memsave
if(mode==131){memsave=1;}	//if using a kinship, just need the decomposition
if(mode==133){memsave=0;}	//will offer memsave
if(mode==138){memsave=1;}	//if using a kinship, just need the decomposition
}

///////////////////////////

