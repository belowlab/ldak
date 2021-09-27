/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//////////////////////////

//set a few parameters - in particular, those used by multiple modes and/or required for parsefiles.c

///////////////////////////

#if MKL==1
mkl_set_num_threads(maxthreads);
#endif

///////////////////////////

if(num_subs==-9999){num_subs=0;}

if(hwestand==-9999)
{
if(nonsnp==0){hwestand=1;}
else{hwestand=0;}
}

if(encoding==-9999){encoding=1;}

if(kindetails==-9999){kindetails=1;}

if(num_regs==-9999){num_regs=0;}

if(pad==-9999){pad=0;}

if(kingz==-9999){kingz=0;}
if(kinraw==-9999){kinraw=0;}

if(discenv==-9999){discenv=0;}

if(allone==-9999){allone=1;}
 
if(strcmp(indhers,"blank")!=0){ignoreweights=1;power=-1;}

if(checkroot==-9999){checkroot=1;}

///////////////////////////

