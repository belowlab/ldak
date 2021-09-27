/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//List of possible modes (D/d indicates data required/possible - K/k, R/r, T/t, S/s same for kinships, regions, tops, summaries; P indicates phenotypes used)

///////////////////////////

//101 - cut-weights D
//102 - calc-weights D
//103 - join-weights D
//104 - calc-weights-all D
//105 - adjust-weights D

//106 - thin D
//107 - thin-tops D
//108 - find tags D
//109 - remove tags D

//111 - cut-kins D
//112 - calc-kins D
//113 - join-kins K
//114 - calc-kins-direct D

//115 - filter K P
//116 - add-grm K
//117 - sub-grm K d
//118 - convert-gz K
//119 - convert-raw K
//120 - calc-sim-grms K

//121 - reml krts P
//122 - calc-blups D k P
//123 - he krt P
//124 - pcgc krt P

//126 - fast-he D P
//127 - fast-pcgc D P

//131 - linear D krts P
//132 - logistic D ts P
//133 - solve-null K rt P

//136 - cut-genes D
//137 - calc-genes-kins D
//138 - calc-genes-reml D krts P
//139 - join-genes-reml d

//141 - calc-tagging D
//142 - join-tagging
//143 - merge-tagging
//144 - reduce-tagging

//146 - sum-hers S
//147 - sum-cors S
//148 - calc-ind-hers
//149 - calc-exps
//150 - calc-posts S

//151 - ridge D P
//152 - bolt D P
//153 - bayesr D P

//156 - calc-cors D
//157 - join-cors
//158 - mega-prs DS P
//159 - pseudo-summaries DS

//161 - pca K
//162 - calc-pca-loads DK
//163 - decompose K
//164 - adjust-grm K

//166 - truncate-grm K
//167 - pca-grm K
//168 - square-grm K
//169 - gxemm-iid K
//170 - gxemm-free K

//171 - calc-stats D
//172 - calc-scores D s P
//173 - make-phenos D
//174 - make-snps

//176 - jackknife
//177 - cut-folds dk
//178 - find-gaussian

X//181 - make-bed D+
X//182 - make-sp D+
X//183 - make-sped D+
X//184 - make-speed D+
X//185 - make-gen D+

X//186 - condense-bed D
X//187 - condense-sp D
X//188 - condense-sped D
X//189 - condense-speed D
X//190 - calc-sim-data DD

///////////////////////////

