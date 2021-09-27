/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Declare variables

///////////////////////////

//parameters set directly by user

//mode indicates which main argument is called (see mode.c for list)
int mode=-9999;

//data input

char folder2[500]="blank", folder[500]="";
char outfile2[500]="blank", outfile[500]="";

char workdir2[500]="blank", workdir[500]="";

char udatafile2[500]="blank", udatafile[500]="";
char ubimfile2[500]="blank", ubimfile[500]="";
char ufamfile2[500]="blank", ufamfile[500]="";
char unamefile2[500]="blank", unamefile[500]="";
char datalist2[500]="blank", datalist[500]="";

int oxchr=-9999;
int genskip=-9999;
int genheaders=-9999;
int genprobs=-9999;
int nonsnp=-9999;

double missingvalue=-9999;

//data filtering

char bsampfile2[500]="blank", bsampfile[500]="";
char csampfile2[500]="blank", csampfile[500]="";
int num_subs=-9999;
char subpref2[500]="blank", subpref[500];

char bpredfile2[500]="blank", bpredfile[500]="";
char cpredfile2[500]="blank", cpredfile[500]="";
int onechr=-9999;
char onesnp[500]="blank";

double minmaf=-9999;
double maxmaf=-9999;
double minvar=-9999;
double minobs=-9999;
double mininfo=-9999;

//data scaling (and pvalues) and coding

char centresfile2[500]="blank", centresfile[500]="";
char weightsfile2[500]="blank", weightsfile[500]="";
int ignoreweights=-9999;
double power=-9999;
int hwestand=-9999;
char pvafile2[500]="blank", pvafile[500]="";

int encoding=-9999;
double threshold=-9999;
double minprob=-9999;

//kinships, regions, responses, summaries and fixed

char kinname2[500]="blank", kinname[500]="";
char kinlist2[500]="blank", kinlist[500]="";
int kindetails=-9999;

int num_regs=-9999;
char regpref2[500]="blank", regpref[500];
double rprune=-9999;

char respfile2[500]="blank", respfile[500]="";
int mpheno=-9999;
char npheno[500]="blank";
int pad=-9999;

char sumsfile2[500]="blank", sumsfile[500]="";
char sums2file2[500]="blank", sums2file[500]="";
int fixn=-9999;
int amb=-9999;
double scaling=-9999, scaling2=-9999;

double prev=-9999;
double ascer=-9999;

char covarfile2[500]="blank", covarfile[500]="";
char topfile2[500]="blank", topfile[500]="";
char envfile2[500]="blank", envfile[500]="";

//calculating weights, thinning and finding/removing tags

int nothin=-9999;
double wprune=-9999;
double window_kb=-9999, window_kb2;
int window_length=-9999, window_length2;
double window_cm=-9999;

double section_kb=-9999;
int section_length=-9999;
double section_cm=-9999;
double buffer_kb=-9999;
int buffer_length=-9999;
double buffer_cm=-9999;

int section=-9999;
int section_start=-9999;

char infosfile2[500]="blank", infosfile[500]="";
int lddecay=-9999;
double halflife=-9999;

int fudge=-9999;
int simplex=-9999;
double maxtime=-9999;
int spread=-9999;

char targetfile2[500]="blank", targetfile[500];

//calculating and manipulating kinships (partitions also used for gene/chunk-based reml and sumher)

int part_length=-9999;
int bychr=-9999;
int num_parts=-9999;
char partpref2[500]="blank", partpref[500];
int checkpart=-9999;

int partition=-9999;
int kingz=-9999;
int kinraw=-9999;
int single=-9999;

double maxrel=-9999, minrel=-9999;
int kinstand=-9999;

//reml, blup and he/pcgc (shortcut used indirectly by other functions)

int shortcut=-9999;
int compact=-9999;
int discenv=-9999;
int liab=-9999;

char remlfile2[500]="blank", remlfile[500]="";

int adjusted=-9999;

int num_vects=-9999;

//association analysis (genes also used for condensing data)

int exact=-9999;
char sampwfile2[500]="blank", sampwfile[500]="";

char genefile2[500]="blank", genefile[500]="";
double chunks=-9999;
int chunksbp=-9999;
int gene_buffer=-9999;
int up_buffer=-9999, down_buffer=-9999;
double minweight=-9999;
int overlap=-9999;

int gene_perms=-9999;
double gprune=-9999;

double cut1=-9999, cut2=-9999;
double gamp=-9999, gam1=-9999, gam2=-9999;

//sumher (annotations also used by fast he and pcgc

int num_anns=-9999;
char annpref2[500]="blank", annpref[500];
int backpart=-9999;
int allone=-9999;
int reduce=-9999;

char printfile2[500]="blank", printfile[500]="";
char herfile2[500]="blank", herfile[500]="";
int unbias=-9999;
int savemat=-9999;
int cover=-9999;

char taglist2[500]="blank", taglist[500]="";
char matlist2[500]="blank", matlist[500]="";
int checkdups=-9999;

char tagfile2[500]="blank", tagfile[500]="";
char altfile2[500]="blank", altfile[500]="";
char cvsfile2[500]="blank", cvsfile[500]="";
char catfile2[500]="blank", catfile[500]="";
char taufile2[500]="blank", taufile[500]="";

int checksums=-9999;
int gcon=-9999, cept=-9999;

int ldsc=-9999;
int chisol=-9999;
int tagone=-9999;

int divide=-9999;
int uptaus=-9999;
char powfile2[500]="blank", powfile[500]="";

int plet=-9999;

char matfile2[500]="blank", matfile[500]="";

char propfile2[500]="blank", propfile[500]="";
char expfile2[500]="blank", expfile[500]="";

//individual-level data prediction, then megaprs

char indhers2[500]="blank", indhers[500]="";
double herscale=-9999;

int restart=-9999;
double cvprop=-9999;
char bvsfile2[500]="blank", bvsfile[500]="";
int skipcv=-9999;

int pointmass=-9999;
int fullspace=-9999;

char corslist2[500]="blank", corslist[500]="";

char corsfile2[500]="blank", corsfile[500]="";
int ptype=-9999;
int onesums=-9999;
int segments=-9999;
int checksize=-9999;

int multiher=-9999;
int basicbolt=-9999;
int ldpred=-9999;

double subprop=-9999;

//pca, decompose and adjust-grm

int axes=-9999;
char pcastem2[500]="blank", pcastem[500]="";

int eigenraw=-9999;
char eigenfile2[500]="blank", eigenfile[500]="";

//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

char scorefile2[500]="blank", scorefile[500]="";
char cofile2[500]="blank", cofile[500]="";
int savecounts=-9999;
char finalfile2[500]="blank", finalfile[500]="";

int num_phenos=-9999;
int num_causals=-9999;
double her=-9999;
double bivar=-9999;
char causalsfile2[500]="blank", causalsfile[500];
char effectsfile2[500]="blank", effectsfile[500];

int num_inds=-9999;
int num_snps=-9999;
double maf1=-9999, maf2=-9999;

char jackfile[500], jackfile2[500]="blank";
char proffile[500], proffile2[500]="blank";
int auc=-9999;

int num_folds=-9999;

char likefile[500], likefile2[500]="blank";

//making and condensing data

int usenames=-9999;
int comsamps=-9999, compreds=-9999;
int exsame=-9999, exdups=-9999;
int passall=-9999;
int speedlong=-9999;

int useminor=-9999;

//common options

int checkroot=-9999;
double mincor=-9999;
double cutoff=-9999;
int constrain=-9999;
int num_blocks=-9999;
int permute=-9999;
double shrink=-9999, strip=-9999;

int bitsize=-9999, bitsize2;
double tol=-9999;
int maxiter=-9999;
int memsave=-9999;

//threading option

int maxthreads=8;

///////////////////////////

//variables which are used for a specific purpose (and set fairly early on)

//storing data

int dtype=-9999, binary=-9999, famhead=-9999, extract=0, use_data, num_files=0;
char **datastems, **bimstems, **famstems;
char datafile[500], bimfile[500], famfile[500];

int num_samples=-9999, num_samples_use=-9999, *keepsamps, **subindex;
int num_preds=-9999, num_preds_use=-9999, *keeppreds, data_length=-9999, *keeppreds_use;

int maxpreds=500000;	//the apriori max number of predictors in gen files

char **ids1, **ids2, **ids3, **allids1, **allids2, **allids3;
int *chr, *allchr;
double *cm, *bp, *cmbp, *allcm, *allbp;
char **preds, **allpreds;
char *al1, *al2, *allal1, *allal2;

unsigned char **data_char;
float *data_single;
double *data, *data2;
double **bytelookup;

float *speedstarts, *speedscales, **ps, **ps2;
double *centres, *mults, *sqdevs, *datasqs;
double *weights, *pvalues;

int num_pows;
double *powers, *powers2;

//kinships, regions, responses, summaries and fixed

int num_kins, *kinnums;
float *kins_single, *kins_single2, **mkins_single;
double *kins, *kins2, **mkins, *kintraces, *kinsums;
char **kinstems;

int *kindex;
double *kcentres, *kmults, *kweights;
char **kpreds, *kal1, *kal2;

int rnum_preds_use, *rkeeppreds, **regindex, rdata_length;
char **rpreds, *ral1, *ral2;
double *rdata;
double *rcentres, *rmults, *rsqdevs;
double *rweights;

int num_resps, num_resps_use=0, *keepresps, *respcounts;
double *resp;

double *nss, *nss2, *rnss, *tnss, *snss, *snss2;
double *chis, *chis2, *rchis, *tchis, *schis, *schis2;
double *rhos, *rhos2, *rrhos, *trhos, *srhos, *srhos2;

int num_covars, num_envs, num_fixed;
double *covar, *thetas, *thetasds, *thetapvas, topher, covher;

int num_tops, *tkeeppreds;
int *tchr;
double *tbp, *tcentres, *tvars;
char **tpreds, *tal1, *tal2;

//calculating weights, thinnings and finding/removing tags

int num_sections;
double decay;
int *sstarts, *sends, *sstarts2, *sends2;
double *infos, *cors, *tally1, *tally2, *tally3, *tally4, *cors2;
int *replace;

int **pindexes, addpart, addpart2, *windex, *vindex;
double **pweights;

int num_tags, num_anals, num_reds, parttype=-9999;
char **tagstems, **matstems;
int *pcounts, *keepparts, *keepparts2;

int num_seek, *sindex, *stypes;

//calculating and manipulating kinships

int *pstarts, *pends;
int *highs, *losts;
double *maxes;

//reml, blup and he/pcgc

char blupfile2[500]="blank", blupfile[500]="";
char regfile2[500]="blank", regfile[500]="";

int **nums, *nummiss;
double **blupcentres, **blupfactors, *bluprands, *blupvalues, *blupmids, *blupprojs;
double **guesses;

int *bstarts, *bends;
double *R, *RTdata, *RTdata2, *mkinsD, *mkinsD2, *mkinsR, *mkinsR2;
double *KKtraces, *KKtraces2, *KYtraces, *KYtraces2;

//association analysis

int *tindex;
double cuts[6]={.1,.01,.001,.0001,.00001,5e-8};
double *sweights;

int num_genes, genemax;
int *gchr, *gstrand, *gstarts, *gends, *gparts;
double *gbp1, *gbp2, *gpvas;
char **gnames;
double *permlike;

//tagging, sumher and bayes factors

double *exps;

int num_parts_use;
int ncv, *cvindex;
double *stags, **svars, **svars2, **ssums, **ssums2, **ssums_use, *cvexps;
char **spreds,  *sal1, *sal2;

int num_means, num_sds, topm, tops;
double *influs, *powlikes, *mvals, *svals, **perfs;

double *props, *pmeans, *pvars, *bfs;

//individual-level data prediction, then megaprs

int num_train, num_test, num_try, restage, recount, *keepboth, *keepboth2;
double *lambdas, *lambdas2, *lambdas3, *lambdas4, *effs, *residuals, *changes, *effs2;

int num_cors;
char **corstems;

int *trytypes; 
double *tryhers, *trylams, *tryscales, *tryps, *tryp2s, *tryp3s, *tryp4s, *tryf2s;
double *loads, *loads2;

size_t *scumsums;
int *maxnums, *actnums, *convs, **bigs;
float **rjks, *rjksums;

double *randnorms, *subrhos, *restrhos, *datarands;

//pca, decompose and adjust-grm

double *kinZ, *ZTkinZ, *W, *WZTkinZ;

//stats, scores, making phenotypes, snps, folds, jackknifing and cors

int num_scores=0;
double **unders, **phens, **liabs;
double *predcors, *predsds, *sP, *sPTP, *sPTP2, *sPTY;
int *sZ1, *sZ2;
double *sX, *sXTX, *sXTX2, *sW;

//making and condensing data

int *Xcurrent;
int *Xnall, *Xnuse, **Xks, **Xks2;
int *XNall, *XNuse, **Xkp, **Xkp2;
char ***Xids1, ***Xids2, ***Xids3;
int **Xchr;
double **Xcm, **Xbp;
char ***Xpreds, **Xal1, **Xal2;
gzFile *Xdatainputgz;

int erra, errb, errc, errd;
char **pid, **mid, **schar, **pchar;
unsigned char onechar;
unsigned short oneshort;
float minfloat, maxfloat;

//common options

int bit, bittotal, bitstart, bitstart2, bitend, bitend2, bitlength, bitlength2, bitmax, step;

double *Y, *Y2, *YTdata;
double *Z, *ZTZ, *ZTZ2, *ZTZ3, *ZZT, *ZTdata;
double *U, *U2, *E, *UTY, *UTZ, *UTdata;

int Xtotal, *Xstarts, *Xends;
double *X, *Xsums, *Xnss, *Xrhos;

///////////////////////////

//generic working variables

size_t scount, smax;
int i, i2, j, j2, j3, k, k2, g, m, p, q, q2, r, s, count, count2, count3, count4;
int current, head, found, total, total2, total3, token, indcount, ecount, wcount, xcount;
int shuffle, start, end, best, mark, gen, gen2, flag, cflag, *order;
double sum, sum2, sumsq, sumsq2, sumsq3, mean, mean2, var, value, value2, value3, value4;
double last, min, max, maf, varphen, gif, postmean, minmse, unifrand;
double *hers, *hersds, *shares, *sharesds, *cohers, *jacks;
double likenull, like, lrtstat, lrtpva;
double neff, neff2, diff, *pens, *likes, *likesold, factor;
double weightsum, minpvalue, minpvalue2;
double *vstarts, *stats, *stats2, **effects;

int one=1, three=3, info, lwork, *iwork, *ifail;
float alpha_single, beta_single;
double alpha, beta, wkopt, *work, vl, vu;

int idshead, *idsorder, *predorder, *usedids, *usedpreds, *indexer, *indexer2, *retain;
char idsfile[500]="blank", **kinids, **kinids2, **kinids3, **wantids;
char **wantpreds, **wantpreds2, **kinpreds;

int readint, readint2, *writeints;
float writefloat;
double readdouble, readdouble2, *writedoubles;
char readchar, readchar2, *rs, readstring[500], readstring2[500], readstring3[500], readstring4[500];
char writestring[500], cmd[500];

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500];
FILE *input, *input2, *input3, **inputs;
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7, *output8, *output9;

DIR *dir;
gzFile datainputgz;
struct stat statstruct;

struct sorting_double *dptrs;

///////////////////////////

