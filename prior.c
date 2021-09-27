//I previously allowed for priors when performing gene-based analysis and REML
//Here are the related scripts

//some declarations in declare.c
double priora=-9999, priorb=-9999, priorc=-9999;

//command line stuff from readargs.c
if(strcmp(argv[count],"--priora")==0)
{
priora=atof(argv[count+1]);found=1;
if(priora<=0){printf("Error, %s should be followed by positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--priorb")==0)
{
priorb=atof(argv[count+1]);found=1;
if(priorb<=0){printf("Error, %s should be followed by positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

if(strcmp(argv[count],"--priorc")==0)
{
priorc=atof(argv[count+1]);found=1;
if(priorc<=0){printf("Error, %s should be followed by positive float (not %s)\n\n", argv[count], argv[count+1]);exit(1);}
}

//requirements from required.c
if(num_kins+num_regs==0)
{
if(priorc!=-9999){printf("Error, you can only use \"--priorc\" if providing kinship matrices or regions\n\n");exit(1);}
if(priora!=-9999||priorb!=-9999)
{
if(priora==-9999||priorb==-9999)
{printf("Error, to perform Bayesian analysis, you must use both \"--priora\" and \"--priorb\"\n\n");exit(1);}
}
}
else
{
if(priora!=-9999||priorb!=-9999||priorc!=-9999)
{
if(priora==-9999||priorb==-9999||priorc==-9999)
{printf("Error, to perform Bayesian analysis, you must use \"--priora\", \"--priorb\" and \"--priorc\"\n\n");exit(1);}
}
}

//told in param.c
if(priora!=-9999)
{printf("Will assume the prior Beta(%.4f,%.4f) (note that for binary traits, this corresponds to observed heritability, not liability heritability)\n\n", priora, priorb);}

//appended in append.c
(void)append_check(projfile,projfile2,workdir);

////////

//some checks in consistent.c

if((priora!=-9999||priorb!=-9999||priorc!=-9999)&&mode!=128)
{printf("Error, you can only use \"--priora\", \"--priorb\" or \"--priorc\" with \"--calc-genes-reml\"\n\n");exit(1);}

///////////////////////////

