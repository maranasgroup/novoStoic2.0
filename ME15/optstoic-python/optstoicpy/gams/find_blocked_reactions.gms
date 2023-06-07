***Identify blocked reactions in database to be excluded for glycolytic pathways
$INLINECOM /*  */
$ONEMPTY

$set dbpath data/optstoic_v3_
***----------------------------------------------------------
Options
        limrow = 5000
        optCR = 1E-6
        optCA = 1E-6
        iterlim = 100000
        decimals = 8
        reslim = 1200
        work = 50000000
        mip = cplex;
***----------------------------------------------------------
Sets
*Metabolites
i
/
$include "%dbpath%metabolites.txt"
/

*Reactions
j
/
$include "%dbpath%reactions.txt"
'EX_glc'
'EX_nad'
'EX_adp'
'EX_phosphate'
'EX_pyruvate'
'EX_nadh'
'EX_atp'
'EX_h2o'
'EX_h+'
'EX_nadp'
'EX_nadph'
/


methglx(j)
/
'R00203'
'R00205'
'R01016'
'R02260'
'R02527'
'R02528'
'R02529'
'R02530'
'R02531'
'R07183'
'R09796'
'R10049'
'R10050'
/


otherUndesirableRxn(j)
/
*atpcycle
'R00088'
'R00122'
'R00123'
'R00125'
'R00127'
*Bicarbonate and pyrrole cycle
'R09794'
'R09795'
*Undesirable glucose uptake loop
'R00305'
'R00874'
'R07359'
'R00837'
'R09749'
'R03075'
'R02985'
'R02558'
'R01555'
'R02727'
'R00010'
'R02778'
'R08946'
'R00306'
'R01139'
'R01138'
/


*GTP/CTP/UTP/ITP/AMP involving reactions
NTP_AMP_rxn(j)
/
$include "data/NTP_and_AMP_reactions.txt"
/

*Reactions involving cofactors only
cofactorOnlyRxn(j)
/
$include "data/cofactor_only_reactions.txt"
**added after analyzing previous result (2016/01/20 6.56pm)
'R10092'
/


atp_irreversible_forward(j)
/
$include "%dbpath%ATP_irreversible_forward_rxns.txt"
/

atp_irreversible_backward(j)
/
$include "%dbpath%ATP_irreversible_backward_rxns.txt"
/
;

Parameters
S(i,j)
/
$include "%dbpath%Sij.txt"
'C00031'.'EX_glc'       -1
'C00003'.'EX_nad'       -1
'C00008'.'EX_adp'       -1
'C00009'.'EX_phosphate' -1
'C00022'.'EX_pyruvate'  -1
'C00004'.'EX_nadh'      -1
'C00002'.'EX_atp'       -1
'C00001'.'EX_h2o'       -1
'C00080'.'EX_h+'        -1
'C00006'.'EX_nadp'      -1
'C00005'.'EX_nadph'     -1
/
*'C00267'.'EX_glc'       -1

rxntype(j)
/
$include "%dbpath%reactiontype.txt"
'EX_glc'       4
'EX_nad'       4
'EX_adp'       4
'EX_phosphate' 4
'EX_pyruvate'  4
'EX_nadh'      4
'EX_atp'       4
'EX_h2o'       4
'EX_h+'        4
'EX_nadp'      4
'EX_nadph'     4
/


LB(j)

UB(j)

dummy(j)
v_up(j)
v_lo(j)

all_excluded_rxns(j)
;

all_excluded_rxns(j) = NTP_AMP_rxn(j) + methglx(j) + cofactorOnlyRxn(j) + otherUndesirableRxn(j);


Scalar
M
/1000/
;

Variables
*** objective function
z
;

Variables
*** flux
v(j)
;


***----------------------------------------------------------
***Setting manually curated rxntype constraints
rxntype(j)$atp_irreversible_forward(j) = 0;
rxntype(j)$atp_irreversible_backward(j) = 2;

*Irreversible forward
LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = M;
*Irreversible backward
LB(j)$(rxntype(j) = 2) = -M;
UB(j)$(rxntype(j) = 2) = 0;
*Reversible
LB(j)$(rxntype(j) = 1) = -M;
UB(j)$(rxntype(j) = 1) = M;

LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = 0;


*** Fix stoichiometry of source/sink metabolites
LB('EX_glc') = -1;
UB('EX_glc') = -1;

LB('EX_pyruvate') = 2;
UB('EX_pyruvate') = 2;

LB('EX_nad') = -2;
UB('EX_nad') = 0;

LB('EX_nadh') = 0;
UB('EX_nadh') = 2;

LB('EX_nadp') = -2;
UB('EX_nadp') = 0;

LB('EX_nadph') = 0;
UB('EX_nadph') = 2;

*generating 0.5 to 5 ATP
LB('EX_adp') = -5;
UB('EX_adp') = -0.5;

LB('EX_phosphate') = -5;
UB('EX_phosphate') = -0.5;

LB('EX_atp') = 0.5;
UB('EX_atp') = 5;

LB('EX_h2o') = 0.5;
UB('EX_h2o') = 5;

LB('EX_h+') = -10;
UB('EX_h+') = 10;

v.lo(j) = LB(j);
v.up(j) = UB(j);

v.fx(j)$all_excluded_rxns(j) = 0;

***----------------------------------------------------------

Equations
obj
stoic
nadphcons1
nadphcons2
nadphcons3
nadphcons4
;

obj.. z =e= sum(j, dummy(j)*v(j));
*obj.. z =e= v('EX_pyruvate');
stoic(i).. sum(j, S(i,j)*v(j)) =e= 0;
nadphcons1.. v('EX_nadph') + v('EX_nadh') =e= 2;
nadphcons2.. v('EX_nadp') + v('EX_nad') =e= -2;
nadphcons3.. v('EX_nadh') + v('EX_nad') =e= 0;
nadphcons4.. v('EX_nadph') + v('EX_nadp') =e= 0;


Model
fva
/
obj
stoic
nadphcons1
nadphcons2
nadphcons3
nadphcons4
/
;
fva.optfile = 1;
fva.holdfixed = 1;

PARAMETER
epsilon /1e-8/

;

alias(j, j1);

file outputfile /optstoic_v3_blocked_reactions_0to5ATP.txt/;
put outputfile;
loop(j1,
         dummy(j) = 0;
         dummy(j1) = 1;
         Solve fva using MIP maximizing z;
         v_up(j1) = z.l;
         Solve fva using MIP minimizing z;
         v_lo(j1) = z.l;
         if ((v_up(j1) < epsilon and v_lo(j1) > -epsilon),
             put  j1.tl:0:30/;
         );
);

putclose outputfile;

