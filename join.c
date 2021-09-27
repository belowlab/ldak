/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//join tagging files (and maybe matrices)

///////////////////////////

//get num_parts, parttype (0=annotations, 1=partitions), addpart, flag (0=MAF, 1=Variance), plus many checks
addpart=0;
for(k=0;k<num_tags;k++)
{
count=countcols(tagstems[k])-9;
if(count<1){printf("Error, %s should have at least 10 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], 9+count);exit(1);}

if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagstems[k]);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring);exit(1);}
if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring3);exit(1);}

if(k==0){flag=(strcmp(readstring2,"Variance")==0);}
if(strcmp(readstring2,"Variance")==0&&flag==0)
{printf("Error, %s was created using \"--hwe-stand NO\", but %s was not\n\n", tagstems[k], tagstems[0]);exit(1);}
if(strcmp(readstring2,"MAF")==0&&flag==1)
{printf("Error, %s was created using \"--hwe-stand NO\", but %s was not\n\n", tagstems[0], tagstems[k]);exit(1);}

for(q=0;q<count;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading element %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", 10+q, tagstems[k]);exit(1);}
}

if(k==0){parttype=(strcmp(readstring,"Base")!=0);}
if(strcmp(readstring,"Base")==0&&parttype==1)
{printf("Error, %s was created using annotations, but %s was not\n\n", tagstems[k], tagstems[0]);exit(1);}
if(strcmp(readstring,"Base")!=0&&parttype==0)
{printf("Error, %s was created using annotations, but %s was not\n\n", tagstems[0], tagstems[k]);exit(1);}

if(parttype==0){addpart2=0;}
else
{
if(strcmp(readstring,"Background")==0){addpart2=1;addpart=1;}
else{addpart2=0;}
}
if(k==0){num_parts=count-addpart2;}
if(count-addpart2!=num_parts)
{
if(parttype==0){printf("Error, the number of annotations used to construct %s (%d) does not match the number used for %s (%d)\n\n", tagstems[k], count-1, tagstems[0], num_parts-1);exit(1);}
else{printf("Error, the number of partitions used to construct %s (%d) does not match the number used for %s (%d)\n\n", tagstems[k], count-addpart2, tagstems[0], num_parts);exit(1);}
}

if(parttype==0)
{
if(count==1){printf("%s provides taggings for the base\n",tagstems[k]);}
if(count==2){printf("%s provides taggings for 1 annotation plus the base\n",tagstems[k]);}
if(count>2){printf("%s provides taggings for %d annotation plus the base\n",tagstems[k], count-1);}
}
else
{
if(count==1){printf("%s provides taggings for 1 partition\n", tagstems[k]);}
else
{
printf("%s provides taggings for %d partitions", tagstems[k], count);
if(addpart2==1){printf(" (including a background partition)");}
printf("\n");
}
}

fclose(input);
}
printf("\n");

if(strcmp(matlist,"blank")!=0)	//check sizes of matrices
{
for(k=0;k<num_tags;k++)
{
if(countcols(matstems[k])!=countcols(tagstems[k])-8)
{printf("Error, %s should have %d columns (not %d), indicating that it does not provide the heritability matrix corresponding to %s\n\n", matstems[k], countcols(tagstems[k])-8, countcols(matstems[k]), tagstems[k]);exit(1);}
}
}

////////

//allocate and set ssums to 0
ssums=malloc(sizeof(double*)*(num_parts+addpart));
for(q=0;q<num_parts+addpart;q++)
{
ssums[q]=malloc(sizeof(double)*(num_parts+addpart+1));
for(q2=0;q2<num_parts+addpart+1;q2++){ssums[q][q2]=0;}
}

if(checkdups==1)	//check no overlap in predictors
{
count=0;for(k=0;k<num_tags;k++){count+=countrows(tagstems[k])-countcols(tagstems[k])+7;}
printf("Checking none of the %d predictors are duplicates; if you are sure this is the case, you can skip this using \"--check-dups NO\"\n", count);

wantpreds=malloc(sizeof(char*)*count);
count=0;
for(k=0;k<num_tags;k++)
{
count2=countrows(tagstems[k])-countcols(tagstems[k])+7;
read_strings(tagstems[k], wantpreds+count, count2, NULL, 1, 1);
count+=count2;
}
check_dups(wantpreds, count, "the tagging files", NULL, 1);
printf("\n");

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
}

//sort out top line of tagging file
sprintf(filename2,"%s.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Neighbours Tagging ");
if(flag==0){fprintf(output2,"Weight MAF Categories Exp_Heritability");}
else{fprintf(output2,"Weight Variance Categories Exp_Heritability");}
for(q=0;q<num_parts;q++)
{
if(parttype==0&&q<num_parts-1){fprintf(output2," Annotation_%d",q+1);}
if(parttype==0&&q==num_parts-1){fprintf(output2," Base");}
if(parttype==1){fprintf(output2," Partition_%d",q+1);}
}
if(addpart==1){fprintf(output2," Background");}
fprintf(output2,"\n");

if(strcmp(matlist,"blank")!=0)	//and also of heritability matrix
{
sprintf(filename4,"%s.matrix",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}

fprintf(output4,"Predictor");
for(q=0;q<num_parts;q++)
{
if(parttype==0&&q<num_parts-1){fprintf(output4," Annotation_%d",q+1);}
if(parttype==0&&q==num_parts-1){fprintf(output4," Base");}
if(parttype==1){fprintf(output4," Partition_%d",q+1);}
}
if(addpart==1){fprintf(output4," Background");}
fprintf(output4,"\n");
}

//now ready to read and print

for(k=0;k<num_tags;k++)
{
//start with tagging file
count=countcols(tagstems[k])-9;
if(count==num_parts+addpart){addpart2=0;}
else{addpart2=1;}
count2=countrows(tagstems[k])-count-2;
if(count2<1)
{printf("Error, %s should have at least %d rows (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], count+3, count2+count+2);exit(1);}

//open and skip first line
if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
printf("Reading tagging details for %d predictors from %s\n", count2, tagstems[k]);

//now read and print next count2 lines, adding on zero if addpart2=1
for(j=0;j<count2;j++)
{
for(q=0;q<9+count;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, tagstems[k]);exit(1);}
fprintf(output2, "%s ", readstring);
}
if(addpart2==1){fprintf(output2, "0.0000");}
fprintf(output2, "\n");
}	//end of j loop

for(q2=0;q2<count;q2++)	//get the sums from the partition lines
{
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\"\n\n", count2+q2+2, tagstems[k]);exit(1);}
for(q=0;q<count;q++)
{
if(fscanf(input, "%lf ", &value)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\"\n\n", q+10, count2+q2+2, tagfile);exit(1);}
ssums[q][q2]+=value;
}
}

//and from the final line (do separately, as must go into num_parts+addpart, and record num preds)
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\"\n\n", count2+count+2, tagstems[k]);exit(1);}
for(q=0;q<count;q++)
{
if(fscanf(input, "%lf ", &value)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\"\n\n", q+10, count2+count+2, tagfile);exit(1);}
ssums[q][num_parts+addpart]+=value;
}

fclose(input);

if(strcmp(matlist,"blank")!=0)	//now read and print heritability matrix
{
count2=countrows(matstems[k])-1;
printf("Reading heritability details for %d predictors from %s\n", count2, tagstems[k]);

//open and skip first line
if((input3=fopen(matstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",matstems[k]);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input3, "%c", &readchar);}

//now read and print next count2 lines
for(j=0;j<count2;j++)
{
for(q=0;q<1+count;q++)
{
if(fscanf(input3, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, matstems[k]);exit(1);}
fprintf(output4, "%s ", readstring);
}
fprintf(output4, "\n");
}

fclose(input3);
}
}	//end of k loop
printf("\n");

//print the sums
for(q2=0;q2<num_parts;q2++)
{
if(parttype==0)
{
if(q2==num_parts-1){fprintf(output2, "The relative contribution of the Base to each category");}
else{fprintf(output2, "The relative contribution of Annotation %d to each category", q2+1);}
}
else{fprintf(output2, "The relative contribution of Partition %d to each category", q2+1);}
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %.4f", ssums[q][q2]);}
fprintf(output2, "\n");
}

if(addpart==1)
{
fprintf(output2, "The relative contribution of the Background to each category");
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %.4f", ssums[q][num_parts]);}
fprintf(output2, "\n");
}

fprintf(output2, "The number of predictors in each of the categories");
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %d", (int)ssums[q][num_parts+addpart]);}
fprintf(output2, "\n");

fclose(output2);
if(strcmp(matlist,"blank")!=0){fclose(output4);}

printf("New tagging file saved in %s", filename2);
if(strcmp(matlist,"blank")!=0){printf(", with new heritability matrix saved in %s", filename4);}
printf("\n\n");

for(q=0;q<num_parts;q++){free(ssums[q]);}free(ssums);

///////////////////////////

