//For each new argument, add checks to consistent.c and possibly required.c

/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Read command line arguments

///////////////////////////

if(argc%2!=1)
{printf("Error, there must be an even number of arguments (not %d); arguments should be provided in pairs\n\n", argc-1);exit(1);}

count=1;
while(count<argc)
{
if(argv[count][0]!='-'||argv[count][1]!='-')
{printf("Error, Argument %d is incorrect (%s); all odd arguments should begin with \"--\"\n\n", count, argv[count]);exit(1);}
if(strcmp(argv[count+1],"blank")==0)
{printf("Error, none of the even arguments can be the word \"blank\"; the universe will probably implode or something\n\n");exit(1);}
count+=2;
}

for(count=1;count<argc;count+=2)
{
for(count2=count+2;count2<argc;count2+=2)
{
if(strcmp(argv[count],argv[count2])==0)
{printf("Error, Arguments %d and %d are the same (%s)\n\n", count, count2, argv[count]);exit(1);}
}
}

if(argc==3){printf("There is one pair of arguments:\n");}
else{printf("There are %d pairs of arguments:\n", (argc-1)/2);}
count=1;
while(count<argc)
{printf("%s %s\n", argv[count], argv[count+1]);count+=2;}
printf("\n");

count=1;
while(count<argc)
{
//see which argument specified, and get the corresponding number/file
found=0;

///////////////////////////
//main arguments

if(strcmp(argv[count],"--cut-weights")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=101;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-weights")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=102;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--join-weights")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=103;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-weights-all")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=104;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--adjust-weights")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=105;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--thin")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=106;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--thin-tops")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=107;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--find-tags")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=108;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--remove-tags")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=109;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--cut-kins")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=111;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-kins")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=112;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--join-kins")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=113;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-kins-direct")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=114;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--filter")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=115;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--add-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=116;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--sub-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=117;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--convert-gz")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=118;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--convert-raw")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=119;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-sim-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=120;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--reml")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=121;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-blups")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=122;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--he")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=123;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--pcgc")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=124;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--fast-he")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=126;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--fast-pcgc")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=127;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--linear")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=131;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--logistic")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=132;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--solve-null")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=133;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--cut-genes")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=136;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-genes-kins")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=137;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-genes-reml")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=138;strcpy(folder2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--join-genes-reml")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=139;strcpy(folder2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--calc-tagging")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=141;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--join-tagging")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=142;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--merge-tagging")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=143;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--reduce-tagging")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=144;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--sum-hers")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=146;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--sum-cors")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=147;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-ind-hers")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=148;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-exps")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=149;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-posts")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=150;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--ridge")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=151;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--bolt")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=152;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--bayesr")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=153;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--calc-cors")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=156;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--join-cors")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=157;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--mega-prs")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=158;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--pseudo-summaries")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=159;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--pca")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=161;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-pca-loads")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=162;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--decompose")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=163;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--adjust-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=164;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--truncate-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=166;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--pca-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=167;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--square-grm")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=168;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--gxemm-iid")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=169;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--gxemm-free")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=170;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--calc-stats")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=171;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-scores")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=172;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--make-phenos")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=173;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--make-snps")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=174;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--jackknife")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=176;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--cut-folds")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=177;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--find-gaussian")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=178;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--make-bed")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=181;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--make-sp")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=182;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--make-sped")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=183;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--make-speed")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=184;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--make-gen")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=185;strcpy(outfile2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--condense-bed")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=186;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--condense-sp")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=187;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--condense-sped")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=188;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--condense-speed")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=189;strcpy(outfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--calc-sim-data")==0)
{
if(mode!=-9999){print_empty();exit(1);}
mode=190;strcpy(outfile2,argv[count+1]);found=1;
}

///////////////////////////
//workdir gets an entry of its own

if(strcmp(argv[count],"--workdir")==0)
{strcpy(workdir2,argv[count+1]);found=1;
}

///////////////////////////
//data input

if(strcmp(argv[count],"--bfile")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(udatafile2,argv[count+1]);
dtype=1;binary=1;found=1;
}

if(strcmp(argv[count],"--sped")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(udatafile2,argv[count+1]);
dtype=3;binary=1;found=1;
}

if(strcmp(argv[count],"--speed")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(udatafile2,argv[count+1]);
dtype=4;binary=1;found=1;
}

////////

if(strcmp(argv[count],"--gen")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(udatafile2,argv[count+1]);
dtype=5;binary=0;found=1;
}

if(strcmp(argv[count],"--chiamo")==0)
{
printf("Sorry, \"--chiamo\" is no longer in use; however, you can achieve the same using \"--gen\", \"--bim\" and \"--fam\"\n\n");exit(1);
}

if(strcmp(argv[count],"--sp")==0)
{
printf("Note that \"--sp %s\" is equivalent to using \"--gen %s.sp\", \"--bim %s.bim\" and \"--fam %s.fam\", with \"--gen-skip 0\", \"--gen-headers 0\" and \"--gen-probs 1\"\nIf the predictor file is gzipped, you should instead use \"--sp-gz\"\n\n", argv[count+1], argv[count+1], argv[count+1], argv[count+1]);
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(unamefile2,argv[count+1]);dtype=11;binary=0;found=1;
}

if(strcmp(argv[count],"--sp-gz")==0)
{
printf("Note that \"--sp-gz %s\" is equivalent to using \"--gen %s.sp.gz\", \"--bim %s.bim\" and \"--fam %s.fam\", with \"--gen-headers 0\" and \"--gen-probs 1\"\n\n", argv[count+1], argv[count+1], argv[count+1], argv[count+1]);
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(unamefile2,argv[count+1]);dtype=12;binary=0;found=1;
}

if(strcmp(argv[count],"--beagle-dose")==0)
{
printf("Note that \"--beagle-dose\" is shorthand for using \"--gen\" with \"--gen-skip 1\", \"--gen-headers 3\" and \"--gen-probs 1\"\n\n");
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(unamefile2,argv[count+1]);dtype=13;binary=0;found=1;
}

if(strcmp(argv[count],"--beagle-probs")==0)
{
printf("Note that \"--beagle-probs\" is shorthand for using \"--gen\" with \"--gen-skip 1\", \"--gen-headers 3\" and \"--gen-probs 3\"\n\n");
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(unamefile2,argv[count+1]);dtype=14;binary=0;found=1;
}

if(strcmp(argv[count],"--haps")==0)
{
printf("Note that \"--haps\" is shorthand for using \"--gen\" with \"--gen-skip 0\", \"--gen-headers 5\" and \"--gen-probs 0\"\nIf %s does not contain predictor details, you should provide these using \"--bim\"\n\n",argv[count]+1);
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
strcpy(unamefile2,argv[count+1]);dtype=15;binary=0;found=1;
}

if(strcmp(argv[count],"--bim")==0)
{strcpy(ubimfile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--fam")==0)
{
if(strcmp(ufamfile2,"blank")!=0)
{printf("Error, you can only use one from \"--fam\" and \"--sample\"\n\n");exit(1);}
strcpy(ufamfile2,argv[count+1]);
famhead=0;found=1;
}

if(strcmp(argv[count],"--sample")==0)
{
if(strcmp(ufamfile2,"blank")!=0)
{printf("Error, you can only use one from \"--fam\" and \"--sample\"\n\n");exit(1);}
if(famhead==0)
{printf("Error, you can not use both \"--sample\" and \"--fam-files YES\"\n\n");exit(1);}
strcpy(ufamfile2,argv[count+1]);
famhead=2;found=1;
}

////////

if(strcmp(argv[count],"--msp")==0)
{
printf("Sorry, \"--msp\" is no longer valid; you should instead use \"--mgen\" with \"--gen-skip 0\", \"--gen-headers 0\" and \"--gen-probs 1\"\n\n");exit(1);
}

if(strcmp(argv[count],"--mbfile")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
dtype=1;binary=1;strcpy(datalist2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--msped")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
dtype=3;binary=1;strcpy(datalist2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--mspeed")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
dtype=4;binary=1;strcpy(datalist2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--mgen")==0)
{
if(dtype!=-9999){printf("Error, you can only use one from \"--bfile\", \"--sped\", \"--speed\", \"--gen\", \"--mbed\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}
dtype=5;binary=0;strcpy(datalist2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--fam-files")==0)
{
if(strcmp(argv[count+1],"YES")==0&&famhead==2)
{printf("Error, you can not use both \"--sample\" and \"--fam-files YES\"\n\n");exit(1);}
famhead=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){famhead=0;}
if(strcmp(argv[count+1],"NO")==0){famhead=2;}
if(famhead==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--oxford-single-chr")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
oxchr=atoi(argv[count+1]);found=1;
if(oxchr<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gen-skip")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
genskip=atoi(argv[count+1]);found=1;
if(genskip<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
if(genskip>1){printf("Warning, %s is typically followed by 0 or 1 (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gen-headers")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
genheaders=atoi(argv[count+1]);found=1;
if(genheaders<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gen-probs")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by 2, 3 or 4 (indicating the number of probabilities provided for each SNP), by 0 (if providing haplotypes) or by 1 (if providing dosages) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
genprobs=atoi(argv[count+1]);found=1;
if(genprobs<0||genprobs>4){printf("Error, %s should be followed by 2, 3 or 4 (indicating the number of probabilities provided for each SNP), by 0 (if providing haplotypes) or by 1 (if providing dosages) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--SNP-data")==0)
{
nonsnp=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){nonsnp=0;}
if(strcmp(argv[count+1],"NO")==0){nonsnp=1;}
if(nonsnp==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--missing-value")==0)
{missingvalue=atof(argv[count+1]);found=1;}

///////////////////////////
//data filtering

if(strcmp(argv[count],"--keep")==0)
{strcpy(bsampfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--remove")==0)
{strcpy(csampfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--subset-number")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_subs=atoi(argv[count+1]);found=1;
if(num_subs<1){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--subset-prefix")==0)
{strcpy(subpref2,argv[count+1]);found=1;
}

////////

if(strcmp(argv[count],"--extract")==0)
{strcpy(bpredfile2,argv[count+1]);extract=1;found=1;}

if(strcmp(argv[count],"--exclude")==0)
{strcpy(cpredfile2,argv[count+1]);extract=1;found=1;}

if(strcmp(argv[count],"--chr")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
onechr=atoi(argv[count+1]);extract=1;found=1;
if(strcmp(argv[count+1],"AUTO")==0){onechr=-1;}
if(onechr<-1){printf("Error, %s should be followed by a non-negative integer or by AUTO (to use only Chr 1-22) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--snp")==0)
{strcpy(onesnp,argv[count+1]);extract=1;found=1;}

////////

if(strcmp(argv[count],"--minmaf")==0||strcmp(argv[count],"--maxmaf")==0||strcmp(argv[count],"--minvar")==0||strcmp(argv[count],"--minobs")==0||strcmp(argv[count],"--mininfo")==0)
{printf("Sorry, \"--minmaf\", \"--maxmaf\", \"--minvar\", \"--minobs\" and \"--mininfo\" have been replaced by \"--min-maf\", \"--max-maf\", \"--min-var\", \"--min-obs\" and \"--min-info\"\n\n");exit(1);}

if(strcmp(argv[count],"--min-maf")==0)
{
minmaf=atof(argv[count+1]);found=1;
if(minmaf<0||minmaf>=0.5){printf("Error, %s should be followed by a float in [0,0.5) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--max-maf")==0)
{
maxmaf=atof(argv[count+1]);found=1;
if(maxmaf<=0||maxmaf>0.5){printf("Error, %s should be followed by a float in (0,.5] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--min-var")==0)
{
minvar=atof(argv[count+1]);found=1;
if(minvar<0){printf("Error, %s should be followed by a non-negative float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--min-obs")==0)
{
minobs=atof(argv[count+1]);found=1;
if(minobs<0||minobs>1){printf("Error, %s should be followed by a non-negative float in [0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--min-info")==0)
{
mininfo=atof(argv[count+1]);found=1;
if(mininfo<0||mininfo>1){printf("Error, %s should be followed by a float in [0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//data scaling (and pvalues) and coding

if(strcmp(argv[count],"--predictor-means")==0)
{strcpy(centresfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--weights")==0)
{strcpy(weightsfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--ignore-weights")==0)
{
ignoreweights=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){ignoreweights=1;}
if(strcmp(argv[count+1],"NO")==0){ignoreweights=0;}
if(ignoreweights==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--power")==0)
{
power=atof(argv[count+1]);found=1;
if(power<-2||power>1){printf("Warning, %s is typically followed by a float within [-2,1] (not %s); the old default was -1, but we now suggest -0.25\n\n", argv[count], argv[count+1]);}
}

if(strcmp(argv[count],"--hwe-stand")==0)
{
hwestand=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){hwestand=1;}
if(strcmp(argv[count+1],"NO")==0){hwestand=0;}
if(hwestand==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--pvalues")==0)
{strcpy(pvafile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--encoding")==0)
{
encoding=-1;found=1;
if(strcmp(argv[count+1],"ADD")==0){encoding=1;}
if(strcmp(argv[count+1],"DOM")==0){encoding=2;}
if(strcmp(argv[count+1],"REC")==0){encoding=3;}
if(strcmp(argv[count+1],"HET")==0){encoding=4;}
if(strcmp(argv[count+1],"MINOR")==0){encoding=5;}
if(strcmp(argv[count+1],"MISS")==0){encoding=6;}
if(encoding==-1){printf("Error, %s should be followed by ADD, DOM, REC, HET, MINOR or MISS (not %s)\nThe additive model codes genotypes 0/1/2, indicating the number of A1 alleles; the dominant, recessive and heterozygous models instead code genotypes 0/2/2, 0/0/2 and 0/2/0, respectively: the missing model codes 1 if the genotype is missing, else 0; the minor allele model matches the additive model, except the counts always correspond to the least-observed allele\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--threshold")==0)
{
if(minprob!=-9999){printf("Error, you can't use both \"--threshold\" and \"--min-prob\"\n\n");exit(1);}
threshold=atof(argv[count+1]);found=1;
if(threshold<0.5||threshold>1){printf("Error, %s should be followed by a float in [0.5,1] (not %s); typical values are 0.9 or 0.95\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--min-prob")==0)
{
if(threshold!=-9999){printf("Error, you cant use both \"--threshold\" and \"--min-prob\"\n\n");exit(1);} 
minprob=atof(argv[count+1]);found=1;
if(minprob!=0&&(minprob<0.5||minprob>=1)){printf("Error, %s should be followed by a float in [0.5,1), or by 0 (sample genotypes based on probabilities) (not %s); typical values are 0.9 or 0.95\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//kinships, regions, responses, summaries and fixed

if(strcmp(argv[count],"--grm")==0)
{
if(strcmp(kinlist2,"blank")!=0){printf("Error, you can't use both \"--grm\" and \"--mgrm\"\n\n");exit(1);}
strcpy(kinname2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--mgrm")==0)
{
if(strcmp(kinname2,"blank")!=0){printf("Error, you can't use both \"--grm\" and \"--mgrm\"\n\n");exit(1);}
strcpy(kinlist2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--kinship-details")==0)
{
kindetails=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kindetails=1;}
if(strcmp(argv[count+1],"NO")==0){kindetails=0;}
if(kindetails==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--region-number")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_regs=atoi(argv[count+1]);found=1;
if(num_regs<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--region-prefix")==0)
{strcpy(regpref2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--region-prune")==0)
{
rprune=atof(argv[count+1]);found=1;
if(rprune<0){printf("Error, %s should be followed by a float within (0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--pheno")==0)
{strcpy(respfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--mpheno")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed a non-negative integer, or by -1 if you want to test each phenotype in turn (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
mpheno=atoi(argv[count+1]);found=1;
if(mpheno<=0&&mpheno!=-1){printf("Error, %s should be followed a non-negative integer, or by -1 if you want to test each phenotype in turn (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--pheno-name")==0)
{
if(strcmp(argv[count+1],"blank")==0)
{printf("Error, you can not use \"--pheno-name blank\"; if a phenotype of this name exists, please instead select it using \"--mpheno\"\n\n");exit(1);}
strcpy(npheno,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--dentist")==0)
{
pad=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){pad=1;}
if(strcmp(argv[count+1],"NO")==0){pad=0;}
if(pad==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--summary")==0)
{strcpy(sumsfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--summary2")==0)
{strcpy(sums2file2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--fixed-n")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
fixn=atoi(argv[count+1]);found=1;
if(fixn<1){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--allow-ambiguous")==0)
{
amb=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){amb=1;}
if(strcmp(argv[count+1],"NO")==0){amb=0;}
if(amb==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--scaling")==0)
{
scaling=atof(argv[count+1]);found=1;
if(scaling<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--scaling2")==0)
{
scaling2=atof(argv[count+1]);found=1;
if(scaling2<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--prevalence")==0)
{
prev=atof(argv[count+1]);found=1;
if(prev<=0||prev>=1){printf("Error, %s should be followed by a float in (0,1) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--ascertainment")==0)
{
ascer=atof(argv[count+1]);found=1;
if(ascer<=0||ascer>=1){printf("Error, %s should be followed by a float in (0,1) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--covar")==0)
{strcpy(covarfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--top-snps")==0)
{printf("Sorry, there was a typo in the supplement; instead of \"--top-snps\" you should use \"--top-preds\"\n\n");exit(1);}

if(strcmp(argv[count],"--enviro")==0)
{strcpy(envfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--top-preds")==0)
{strcpy(topfile2,argv[count+1]);found=1;}

///////////////////////////
//calculating weights, thinning and finding/removing tags

if(strcmp(argv[count],"--no-thin")==0)
{
nothin=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){nothin=1;}
if(strcmp(argv[count+1],"NO")==0){nothin=0;}
if(strcmp(argv[count+1],"DONE")==0){nothin=2;}
if(nothin==-1){printf("Error, %s should be followed by YES, NO or DONE (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--window-prune")==0)
{
wprune=atof(argv[count+1]);found=1;
if(wprune<0||wprune>1){printf("Error, %s should be followed by a float within (0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--window-kb")==0)
{
window_kb=atof(argv[count+1]);found=1;
if(window_kb<0){printf("Error, %s should be followed by a non-negative float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--window-length")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
window_length=atoi(argv[count+1]);found=1;
if(window_length<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--window-cm")==0)
{
window_cm=atof(argv[count+1]);found=1;
if(window_cm<0){printf("Error, %s should be followed by a non-negative float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--section-kb")==0)
{
section_kb=atof(argv[count+1]);found=1;
if(section_kb<=0){printf("Error, %s should be followed by a positive float (not %s)\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--section-length")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count]);exit(1);}
section_length=atoi(argv[count+1]);found=1;
if(section_length<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--section-cm")==0)
{
section_cm=atof(argv[count+1]);found=1;
if(section_cm<=0){printf("Error, %s should be followed by a positive float (not %s)\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--section-buffer")==0)
{
printf("Error, \"--section-buffer\" is no longer valid; you should instead use \"--buffer-kb\" or \"--buffer-length\" to specify the buffer in kb or as a fixed number of predictors\n\n");exit(1);
}

if(strcmp(argv[count],"--buffer-kb")==0)
{
buffer_kb=atof(argv[count+1]);found=1;
if(buffer_kb<0){printf("Error, %s should be followed by a non-negative float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--buffer-length")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
buffer_length=atoi(argv[count+1]);found=1;
if(buffer_length<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--buffer-cm")==0)
{
buffer_cm=atof(argv[count+1]);found=1;
if(buffer_cm<0){printf("Error, %s should be followed by a non-negative float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--section")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
section=atoi(argv[count+1]);found=1;
if(section<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--start-section")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
section_start=atoi(argv[count+1]);found=1;
if(section_start<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--infos")==0)
{strcpy(infosfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--decay")==0)
{
lddecay=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){lddecay=1;}
if(strcmp(argv[count+1],"NO")==0){lddecay=0;}
if(lddecay==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--halflife")==0)
{printf("Sorry, \"--halflife\" has been replaced by \"--half-life\"\n\n");exit(1);}

if(strcmp(argv[count],"--half-life")==0)
{
halflife=atof(argv[count+1]);found=1;
if(halflife<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--quick-weights")==0)
{
fudge=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){fudge=1;}
if(strcmp(argv[count+1],"NO")==0){fudge=0;}
if(fudge==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--simplex")==0)
{
simplex=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){simplex=1;}
if(strcmp(argv[count+1],"NO")==0){simplex=0;}
if(simplex==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--maxtime")==0)
{printf("Sorry, \"--maxtime\" has been replaced by \"--max-time\"\n\n");exit(1);}

if(strcmp(argv[count],"--max-time")==0)
{
maxtime=atof(argv[count+1]);found=1;
if(maxtime<=0){printf("Error, %s should be followed by a positive float (not %s) providing the maximum time in minutes\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--spread")==0)
{
spread=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){spread=1;}
if(strcmp(argv[count+1],"NO")==0){spread=0;}
if(spread==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--targets")==0)
{strcpy(targetfile2,argv[count+1]);found=1;}

///////////////////////////
//calculating and manipulating kinships

if(strcmp(argv[count],"--partition-length")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
part_length=atoi(argv[count+1]);found=1;
if(part_length<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--by-chr")==0)
{
bychr=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){bychr=1;}
if(strcmp(argv[count+1],"NO")==0){bychr=0;}
if(bychr==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--partition-number")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_parts=atoi(argv[count+1]);found=1;
if(num_parts<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--partition-prefix")==0)
{strcpy(partpref2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--check-partitions")==0)
{
checkpart=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){checkpart=1;}
if(strcmp(argv[count+1],"NO")==0){checkpart=0;}
if(checkpart==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--partition")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
partition=atoi(argv[count+1]);found=1;
if(partition<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--kinship-gz")==0)
{
kingz=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kingz=1;}
if(strcmp(argv[count+1],"NO")==0){kingz=0;}
if(kingz==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--kinship-raw")==0)
{
kinraw=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kinraw=1;}
if(strcmp(argv[count+1],"NO")==0){kinraw=0;}
if(kinraw==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--single")==0)
{
single=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){single=1;}
if(strcmp(argv[count+1],"NO")==0){single=0;}
if(single==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--max-rel")==0)
{
maxrel=atof(argv[count+1]);found=1;
if(maxrel<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--min-rel")==0)
{
minrel=atof(argv[count+1]);found=1;
if(minrel<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--kin-stand")==0)
{
kinstand=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kinstand=1;}
if(strcmp(argv[count+1],"NO")==0){kinstand=0;}
if(kinstand==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//reml, blup and he/pcgc

if(strcmp(argv[count],"--shortcut")==0)
{
shortcut=-1;found=1;
if(strcmp(argv[count+1],"NO")==0){shortcut=0;}
if(shortcut==-1){printf("Error, %s must be followed by NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--compact")==0)
{
compact=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){compact=1;}
if(strcmp(argv[count+1],"NO")==0){compact=0;}
if(compact==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--subgroups")==0)
{
discenv=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){discenv=1;}
if(strcmp(argv[count+1],"NO")==0){discenv=0;}
if(discenv==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--liability")==0)
{
printf("Sorry, you can't (currently) use \"--liability\"\n\n");exit(1);
liab=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){liab=1;}
if(strcmp(argv[count+1],"NO")==0){liab=0;}
if(liab==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--remlfile")==0)
{strcpy(remlfile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--adjusted")==0)
{
adjusted=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){adjusted=1;}
if(strcmp(argv[count+1],"NO")==0){adjusted=0;}
if(adjusted==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--repetitions")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_vects=atoi(argv[count+1]);found=1;
if(num_vects<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//association analysis

if(strcmp(argv[count],"--exact")==0)
{
exact=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){exact=1;}
if(strcmp(argv[count+1],"NO")==0){exact=0;}
if(exact==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--sample-weights")==0)
{strcpy(sampwfile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--genefile")==0)
{
if(strcmp(genefile2,"blank")!=0||chunks!=-9999||chunksbp!=-9999)
{printf("Error, you can only use one from \"--genefile\", \"--chunks\" or \"--chunks-bp\"\n\n");exit(1);}
strcpy(genefile2,argv[count+1]);found=1;
}

if(strcmp(argv[count],"--chunks")==0)
{
if(strcmp(genefile2,"blank")!=0||chunks!=-9999||chunksbp!=-9999)
{printf("Error, you can only use one from \"--genefile\", \"--chunks\" or \"--chunks-bp\"\n\n");exit(1);}
chunks=atof(argv[count+1]);found=1;
if(chunks<1){printf("Error, %s should be followed by a float at least 1 (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--chunks-bp")==0)
{
if(strcmp(genefile2,"blank")!=0||chunks!=-9999||chunksbp!=-9999)
{printf("Error, you can only use one from \"--genefile\", \"--chunks\" or \"--chunks-bp\"\n\n");exit(1);}
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
chunksbp=atoi(argv[count+1]);found=1;
if(chunksbp<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gene-buffer")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
gene_buffer=atoi(argv[count+1]);found=1;
if(gene_buffer<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--up-buffer")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
up_buffer=atoi(argv[count+1]);found=1;
if(up_buffer<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--down-buffer")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
down_buffer=atoi(argv[count+1]);found=1;
if(down_buffer<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--min-weight")==0)
{
minweight=atof(argv[count+1]);found=1;
if(minweight<0){printf("Error, %s should be followed by a positive float or zero (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--overlap")==0)
{
overlap=-1;found=1;
if(strcmp(argv[count+1],"NO")==0){overlap=0;}
if(strcmp(argv[count+1],"YES")==0){overlap=1;}
if(overlap==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--gene-permutations")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
gene_perms=atoi(argv[count+1]);found=1;
if(gene_perms<0){printf("Error, %s should be followed by a non-negative integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gene-prune")==0)
{
gprune=atof(argv[count+1]);found=1;
if(gprune<0){printf("Error, %s should be followed by a float within (0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--sig1")==0)
{
cut1=atof(argv[count+1]);found=1;
if(cut1<=0||cut1>=1){printf("Error, %s should be followed by a float within (0,1) (not %s)\n\n", argv[count],  argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--sig2")==0)
{
cut2=atof(argv[count+1]);found=1;
if(cut2<=0||cut2>1){printf("Error, %s should be followed by a float within (0,1] (not %s)\n\n", argv[count],  argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gamma-fraction")==0)
{
gamp=atof(argv[count+1]);found=1;
if(gamp<=0||gamp>1){printf("Error, %s should be followed by a float within (0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gamma-alpha")==0)
{
gam1=atof(argv[count+1]);found=1;
if(gam1<=0){printf("Error, %s should be followed by positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--gamma-beta")==0)
{
gam2=atof(argv[count+1]);found=1;
if(gam2<=0){printf("Error, %s should be followed by positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//sumher

if(strcmp(argv[count],"--annotation-number")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_anns=atoi(argv[count+1]);found=1;
if(num_anns<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--annotation-prefix")==0)
{strcpy(annpref2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--background")==0)
{
backpart=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){backpart=1;}
if(strcmp(argv[count+1],"NO")==0){backpart=0;}
if(backpart==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--all-one")==0)
{
allone=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){allone=1;}
if(strcmp(argv[count+1],"NO")==0){allone=0;}
if(allone==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--reduce")==0)
{
reduce=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){reduce=1;}
if(strcmp(argv[count+1],"NO")==0){reduce=0;}
if(reduce==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--regression-predictors")==0)
{strcpy(printfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--heritability-predictors")==0)
{strcpy(herfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--unbiased")==0)
{
unbias=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){unbias=1;}
if(strcmp(argv[count+1],"NO")==0){unbias=0;}
if(unbias==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--save-matrix")==0)
{
savemat=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){savemat=1;}
if(strcmp(argv[count+1],"NO")==0){savemat=0;}
if(savemat==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--coverage")==0)
{
cover=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){cover=1;}
if(strcmp(argv[count+1],"NO")==0){cover=0;}
if(cover==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--taglist")==0)
{strcpy(taglist2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--matlist")==0)
{strcpy(matlist2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--check-dups")==0)
{
checkdups=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){checkdups=1;}
if(strcmp(argv[count+1],"NO")==0){checkdups=0;}
if(checkdups==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--tagfile")==0)
{strcpy(tagfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--alternative-tags")==0)
{strcpy(altfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--cv-predictors")==0)
{strcpy(cvsfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--categories")==0)
{strcpy(catfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--taus")==0)
{strcpy(taufile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--check-sums")==0)
{
checksums=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){checksums=1;}
if(strcmp(argv[count+1],"NO")==0){checksums=0;}
if(checksums==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--genomic-control")==0)
{
gcon=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){gcon=1;}
if(strcmp(argv[count+1],"NO")==0){gcon=0;}
if(gcon==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--intercept")==0)
{
cept=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){cept=1;}
if(strcmp(argv[count+1],"NO")==0){cept=0;}
if(cept==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--LDSC")==0)
{
ldsc=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){ldsc=1;}
if(strcmp(argv[count+1],"NO")==0){ldsc=0;}
if(ldsc==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--chisq-solver")==0)
{
chisol=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){chisol=1;}
if(strcmp(argv[count+1],"NO")==0){chisol=0;}
if(chisol==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--tag-one")==0)
{
tagone=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){tagone=1;}
if(strcmp(argv[count+1],"NO")==0){tagone=0;}
if(tagone==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--divisions")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
divide=atoi(argv[count+1]);found=1;
if(divide<1){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--update-taus")==0)
{
uptaus=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){uptaus=1;}
if(strcmp(argv[count+1],"NO")==0){uptaus=0;}
if(uptaus==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--powerfile")==0)
{strcpy(powfile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--pleiotrophy")==0)
{
plet=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){plet=1;}
if(strcmp(argv[count+1],"NO")==0){plet=0;}
if(plet==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--matrix")==0)
{strcpy(matfile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--props")==0)
{strcpy(propfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--expectations")==0)
{strcpy(expfile2,argv[count+1]);found=1;}

///////////////////////////
//individual-level data prediction, then megaprs

if(strcmp(argv[count],"--ind-hers")==0)
{strcpy(indhers2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--her-scale")==0)
{
herscale=atof(argv[count+1]);found=1;
if(herscale<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--restart")==0)
{
restart=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){restart=1;}
if(strcmp(argv[count+1],"NO")==0){restart=0;}
if(restart==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--cv-proportion")==0)
{
cvprop=atof(argv[count+1]);found=1;
if(cvprop<=0||cvprop>=1){printf("Error, %s should be followed by a float within (0,1) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--cv-samples")==0)
{strcpy(bvsfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--skip-cv")==0)
{
skipcv=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){skipcv=1;}
if(strcmp(argv[count+1],"NO")==0){skipcv=0;}
if(skipcv==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--point-mass")==0)
{
pointmass=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){pointmass=1;}
if(strcmp(argv[count+1],"NO")==0){pointmass=0;}
if(pointmass==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--full-space")==0)
{
fullspace=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){fullspace=1;}
if(strcmp(argv[count+1],"NO")==0){fullspace=0;}
if(fullspace==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--corslist")==0)
{strcpy(corslist2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--cors")==0)
{strcpy(corsfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--model")==0)
{
ptype=-1;found=1;
if(strcmp(argv[count+1],"mega")==0){ptype=0;}
if(strcmp(argv[count+1],"lasso-sparse")==0){ptype=1;}
if(strcmp(argv[count+1],"lasso")==0){ptype=2;}
if(strcmp(argv[count+1],"ridge")==0){ptype=3;}
if(strcmp(argv[count+1],"bolt")==0){ptype=4;}
if(strcmp(argv[count+1],"bayesr")==0){ptype=5;}
if(strcmp(argv[count+1],"bayesr-shrink")==0){ptype=6;}
if(ptype==-1){printf("Error, %s should be followed by \"lasso\", \"ridge\", \"bolt\", \"bayesr\", \"bayesr-shrink\", \"lasso-sparse\"  or \"--mega\" (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--one-sums")==0)
{
onesums=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){onesums=1;}
if(strcmp(argv[count+1],"NO")==0){onesums=0;}
if(onesums==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--segments")==0)
{
if(log2(atof(argv[count+1]))!=(int)log2(atof(argv[count+1])))
{printf("Error, %s should be followed by a power of two (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
segments=atoi(argv[count+1]);found=1;
if(segments<=0){printf("Error, %s should be followed by a power of two (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--check-sizes")==0)
{
checksize=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){checksize=1;}
if(strcmp(argv[count+1],"NO")==0){checksize=0;}
if(checksize==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--multi-hers")==0)
{
multiher=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){multiher=1;}
if(strcmp(argv[count+1],"NO")==0){multiher=0;}
if(multiher==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--basic-bolt")==0)
{
basicbolt=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){basicbolt=1;}
if(strcmp(argv[count+1],"NO")==0){basicbolt=0;}
if(basicbolt==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--LDpred")==0)
{
ldpred=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){ldpred=1;}
if(strcmp(argv[count+1],"NO")==0){ldpred=0;}
if(ldpred==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--training-proportion")==0)
{
subprop=atof(argv[count+1]);found=1;
if(subprop<=0||subprop>=1){printf("Error, %s should be followed by a float within (0,1) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//pca, decompose and adjust-grm

if(strcmp(argv[count],"--axes")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
axes=atoi(argv[count+1]);found=1;
if(axes<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--pcastem")==0)
{strcpy(pcastem2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--eigen-raw")==0)
{
eigenraw=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){eigenraw=1;}
if(strcmp(argv[count+1],"NO")==0){eigenraw=0;}
if(eigenraw==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--eigen")==0)
{strcpy(eigenfile2,argv[count+1]);found=1;}

///////////////////////////
//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

if(strcmp(argv[count],"--scorefile")==0)
{strcpy(scorefile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--coeffsfile")==0)
{strcpy(cofile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--save-counts")==0)
{
savecounts=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){savecounts=1;}
if(strcmp(argv[count+1],"NO")==0){savecounts=0;}
if(savecounts==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--final-effects")==0)
{strcpy(finalfile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--num-phenos")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_phenos=atoi(argv[count+1]);found=1;
if(num_phenos<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--num-causals")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer or by -1 (every predictor causal) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_causals=atoi(argv[count+1]);found=1;
if(num_causals<=0&&num_causals!=-1){printf("Error, %s should be followed by a positive integer or by -1 (every predictor causal) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--her")==0)
{
her=atof(argv[count+1]);found=1;
if(her>1){printf("Error, %s should be followed by a float within [0,1] (not %s); if you are sure you want a value greater than 1, use \"--her-big\"\n\n", argv[count], argv[count+1]);exit(1);}
if(her<0){printf("Error, %s should be followed by a float within [0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--her-big")==0)
{
her=atof(argv[count+1]);found=1;
if(her<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--bivar")==0)
{
bivar=atof(argv[count+1]);found=1;
if(bivar<-1||bivar>1){printf("Error, %s should be followed by a float within [-1,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--causals")==0)
{strcpy(causalsfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--effects")==0)
{strcpy(effectsfile2,argv[count+1]);found=1;}

////////

if(strcmp(argv[count],"--num-samples")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_inds=atoi(argv[count+1]);found=1;
if(num_inds<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--num-snps")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_snps=atoi(argv[count+1]);found=1;
if(num_snps<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--maf-low")==0)
{
maf1=atof(argv[count+1]);found=1;
if(maf1<0||maf1>=0.5){printf("Error, %s should be followed by a float in [0,0.5) (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--maf-high")==0)
{
maf2=atof(argv[count+1]);found=1;
if(maf2<=0||maf2>0.5){printf("Error, %s should be followed by a float in (0,0.5] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}
////////

if(strcmp(argv[count],"--data-pairs")==0)
{strcpy(jackfile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--profile")==0)
{strcpy(proffile2,argv[count+1]);found=1;}

if(strcmp(argv[count],"--AUC")==0)
{
auc=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){auc=1;}
if(strcmp(argv[count+1],"NO")==0){auc=0;}
if(auc==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--num-folds")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by an integer greater than 1 (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_folds=atoi(argv[count+1]);found=1;
if(num_folds<2){printf("Error, %s should be followed by an integer greater than 1 (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--likelihoods")==0)
{strcpy(likefile2,argv[count+1]);found=1;}

///////////////////////////
//making and condensing data

if(strcmp(argv[count],"--use-names")==0)
{
usenames=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){usenames=1;}
if(strcmp(argv[count+1],"NO")==0){usenames=0;}
if(usenames==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--common-samples")==0)
{
comsamps=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){comsamps=1;}
if(strcmp(argv[count+1],"NO")==0){comsamps=0;}
if(comsamps==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--common-preds")==0)
{
compreds=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){compreds=1;}
if(strcmp(argv[count+1],"NO")==0){compreds=0;}
if(compreds==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--exclude-same")==0)
{
exsame=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){exsame=1;}
if(strcmp(argv[count+1],"NO")==0){exsame=0;}
if(exsame==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--exclude-dups")==0)
{
exdups=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){exdups=1;}
if(strcmp(argv[count+1],"NO")==0){exdups=0;}
if(exdups==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--pass-all")==0)
{
passall=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){passall=1;}
if(strcmp(argv[count+1],"NO")==0){passall=0;}
if(passall==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--speed-long")==0)
{
speedlong=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){speedlong=1;}
if(strcmp(argv[count+1],"NO")==0){speedlong=0;}
if(speedlong==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--count-minor")==0)
{
useminor=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){useminor=1;}
if(strcmp(argv[count+1],"NO")==0){useminor=0;}
if(useminor==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////
//common options

if(strcmp(argv[count],"--check-root")==0)
{
checkroot=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){checkroot=1;}
if(strcmp(argv[count+1],"NO")==0){checkroot=0;}
if(checkroot==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--mincor")==0)
{printf("Sorry, \"--mincor\" has been replaced by \"--min-cor\"\n\n");exit(1);}

if(strcmp(argv[count],"--min-cor")==0)
{
mincor=atof(argv[count+1]);found=1;
if(mincor<0||mincor>1){printf("Error, %s should be followed by a float in [0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--cutoff")==0)
{cutoff=atof(argv[count+1]);found=1;}

if(strcmp(argv[count],"--constrain")==0)
{
constrain=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){constrain=1;}
if(strcmp(argv[count+1],"NO")==0){constrain=0;}
if(constrain==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--num-blocks")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer, or by -1 if you want to do leave-one-sample-out jackknifing (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
num_blocks=atoi(argv[count+1]);found=1;
if(num_blocks<=2&&num_blocks!=-1){printf("Error, %s should be followed by an integer greater than 2, or by -1 if you want to do leave-one-sample-out jackknifing (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--permute")==0)
{
permute=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){permute=1;}
if(strcmp(argv[count+1],"NO")==0){permute=0;}
if(permute==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--shrink")==0)
{
shrink=atof(argv[count+1]);found=1;
if(shrink<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--strip")==0)
{
strip=atof(argv[count+1]);found=1;
if(strip<=0||strip>1){printf("Error, %s should be followed by a float within (0,1] (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

////////

if(strcmp(argv[count],"--bit-size")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
bitsize=atoi(argv[count+1]);found=1;
if(bitsize<1){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--tolerance")==0)
{
tol=atof(argv[count+1]);found=1;
if(tol<=0){printf("Error, %s should be followed by a positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--maxiter")==0)
{printf("Sorry, \"--maxiter\" has been replaced by \"--max-iter\"\n\n");exit(1);}

if(strcmp(argv[count],"--max-iter")==0)
{
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
maxiter=atoi(argv[count+1]);found=1;
if(maxiter<=0){printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--memory-save")==0)
{
memsave=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){memsave=1;}
if(strcmp(argv[count+1],"NO")==0){memsave=0;}
if(memsave==-1){printf("Error, %s should be followed by YES or NO (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////

//threading option

if(strcmp(argv[count],"--max-threads")==0)
{
//if(MKL==0){printf("Sorry, you can not use \"--max-threads\" in this version of LDAK; please omit or ask Doug for a revised version\n\n");exit(1);}
if(atof(argv[count+1])!=atoi(argv[count+1]))
{printf("Error, %s should be followed by a positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
maxthreads=atoi(argv[count+1]);found=1;
if(maxthreads<=0){printf("Error, %s should be followed positive integer (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

///////////////////////////

if(found==0){printf("Error, %s is not a recognised argument\n",argv[count]);exit(1);}
count+=2;
}

///////////////////////////

if(mode==-9999)
{
printf("Error, each command requires a main argument; here is a list of main arguments\n\n");

//currently do not advertise adjust-weights, truncate-grm, pca-grm or square-grm

printf("Heritability analysis using individual-level data:\n\n");

printf("--cut-weights <folder> - cut predictors into sections ready for calculating weightings\n\
--calc-weights <folder> - calculate weightings for one section\n\
--join-weights <folder> - join weightings across sections\n\
--calc-weights-all <folder> - calculate weightings for all sections in one step\n\n");

printf("--cut-kins <folder> - cut predictors into partitions ready for calculating kinships\n\
--calc-kins <folder> - calculate kinship matrix for one partition\n\
--join-kins <folder> - join kinship matrices across partitions\n\
--calc-kins-direct <output> - calculate kinship matrix in one step\n\n");

printf("--reml <output> - regress phenotype on kinships and/or regions\n\
--calc-blups <output> - calculate BLUP effect size estimates (using results from --reml)\n\
--he <output> - perform Haseman-Elston Regression (an alternative to --reml)\n\
--pcgc <output> - perform PCGC Regression (recommended instead of --reml for binary traits)\n\
--fast-he <output> - perform fast (but approximate) Haseman-Elston Regression\n\
--fast-pcgc <output> - perform fast (but approximate) PCGC Regression\n\n");

printf("Association testing:\n\n");

printf("--linear <output> - perform standard or mixed-model linear regression\n\
--logistic <output> - perform logistic regression\n\
--solve-null <output> - perform two-part mixed-model linear regression (for complex models)\n\n");

printf("--cut-genes <folder> - cut predictors into genes (or fixed-length chunks)\n\
--calc-genes-kins <folder> - calculate kinships for each gene/chunk in one partition\n\
--calc-genes-reml <folder> - calculate association test for each gene in one partition\n\
--join-genes-reml <folder> - join association test results across partitions (also computes more accurate p-values)\n\n");

printf("Heritability analysis using summary statistics:\n\n");

printf("--calc-tagging <output> - calculate tagging file for use with --sum-hers or --sum-cors\n\
--join-tagging <output> - join tagging files across predictors\n\
--merge-tagging <output> - combine tagging files across categories\n\
--reduce-tagging <output> - extract categories from a tagging file\n\n");

printf("--sum-hers <output> - estimate heritabilities from summary statistics\n\
--sum-cors <output> - estimate genetic correlations from summary statistics\n\
--calc-ind-hers <output> - calculate expected heritability contributed by predictors (using results from --sum-hers)\n\
--calc-exps <output> - calculate expected heritability tagged by predictors (using results from --sum-hers)\n\
--calc-posts <output> - calculate posterior estimate of predictor effect sizes (using results from --calc-exps)\n\n");

printf("Prediction:\n\n");

printf("--ridge <output> - construct a ridge regression prediction model\n\
--bolt <output> - construct a Bolt prediction model\n\
--bayesr <output> -  construct a BayesR prediction model\n\n");

printf("--calc-cors <output> - calculate predictor-predictor correlations for use with --mega-prs\n\
--join-cors <output> - join correlations across predictors\n\
--mega-prs <output> - construct lasso, ridge regression, Bolt and BayesR prediction models from summary statistics\n\
--pseudo-summaries <output> - generate partial summary statistics (used to train MegaPRS prediction models)\n\n");

printf("Other features:\n\n");

printf("--thin <output> - prune predictors based on pairwise correlations\n\
--thin-tops <output> - prune highly-associated predictors based on pairwise correlations\n\
--find-tags <output> - search for the predictors most correlated with those in a scorefile\n\
--remove-tags <output> - search for all predictors correlated with those in a target file\n\n");

printf("--filter <output> - fiter samples based on relatedness\n\
--add-grm <output> - combine kinship matrices\n\
--sub-grm <output> - subtract from one kinship matrix the contributions of other kinship matrices\n\
--convert-gz <output> - convert kinship matrix saved in gz (old GCTA) format\n\
--convert-raw <output> - convert kinship matrix saved as a text file\n\
--calc-sim-grms <output> - calculate correlations between pairs of kinship matrices\n\n");

printf("--pca <output> - compute the principal component axes of a kinship matrix\n\
--calc-pca-loads <output> - calculate the principal component predictor loadings (using results from --pca)\n\
--decompose <output> - eigen-decompose a kinship matrix\n\
--adjust-grm <output> - adjust a kinship matrix for covariates\n\
--gxemm-iid / --gxemm-free <output> - calculate environmental kinship matrices\n\n");

printf("--calc-stats <output> - calculate predictor allele frequencies, call-rates and possibly information scores\n\
--calc-scores <output> - calculate one or more polygenic risk scores\n\
--make-phenos <output> - simulate phenotypes\n\
--make-snps <output> - simulate SNP data\n\n");

printf("--jackknife <output> - measure prediction accuracy, obtaining standard deviations via block jackknifing\n\
--cut-folds <output> - divide samples in preparation for performing K-fold cross-validation\n\
--find-gaussian <output> - estimate the selection-related parameter alpha (using results from --sum-hers)\n\n");

printf("--make-bed / --make-sp / --make-sped / --make-speed / --make-gen <outfile> - convert to bed / sp / sped / speed / gen format\n\
--condense-bed / condense-sp / condense-sped / condense-speed <outfile> - condense to bed / sp / sped / speed format\n\
--calc-sim-data <outfile> - calculate concordance between two datasets\n\n");

exit(1);
}

///////////////////////////

/*
//to come (maybe)

if(strcmp(argv[count],"--calc-genes-bivar")==0)
{}

if(strcmp(argv[count],"--join-genes-bivar")==0)
{}

if(strcmp(argv[count],"--decompose-bivar")==0)
{}


if(strcmp(argv[count],"--calc-blups-bivar")==0)
{}

if(strcmp(argv[count],"--reml-bivar-mcmc")==0)
{}

if(strcmp(argv[count],"--bivar-null")==0)
{}

*/

