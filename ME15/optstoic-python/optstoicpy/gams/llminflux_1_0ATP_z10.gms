***This is use to identify glycolytic pathway (C00031 to C00022)
$INLINECOM /*  */
$ONEMPTY
$set dbpath data/optstoic_v3_
***------------------CHANGE_THIS_SECTION---------------------
Scalar
nATP /1.0/
zlb /10/
;

$set outputfname llminflux_1_0ATP_z10
FILE outputfile /%outputfname%.json/;

*///---------------------------------------------------------

***----------------------------------------------------------
Options
        limrow = 5000
        optCR = 1E-6
        optCA = 1E-6
        iterlim = 100000
        decimals = 8
        reslim = 1800
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

blocked(j)
/
$include "%dbpath%blocked_reactions_0to5ATP_excluderxns_20161024.txt"
/

atp_irreversible_forward(j)
/
$include "%dbpath%ATP_irreversible_forward_rxns.txt"
/

atp_irreversible_backward(j)
/
$include "%dbpath%ATP_irreversible_backward_rxns.txt"
/

jint(j)
/
$include "%dbpath%reactions.txt"
/

l
/
$include "%dbpath%loops_nocofactor_20161025.txt"
/

jloop(j)

*GTP/CTP/UTP/ITP/AMP involving reactions
;

*reaction involves in internal cycle/loops
jloop(j) = jint(j) - blocked(j);

***------------------INTEGER_CUT_SECTION---------------------
***Mode 1: No integer cut

Set
k
/1*1000/
;

Parameters
**For adding integer cut from previous run
maxiter /1000/
**Iter is the starting iteration
iter /1/
**For adding integer cut from previous run
store(k,j)
onConstraint(k)

*initial objective value  (after first iteration)
*which is also the lower bound of objective value

*changes in lower bound of z
*dz /0/
*dz should change every time changeFlag = xx
*xx /200/

*flag change lb
*changeFlag /1/
;

store(k,j) = 0;
onConstraint(k)=0;

$ontext

***Mode 2: With integer cut

Set
k
/1*1000/
**For adding integer cut from previous run
previousk(k)
/1*570/
;

Parameters
**For adding integer cut from previous run
maxiter /1000/
*Iter is the starting iteration
iter /101/
**For adding integer cut from previous run
store(k,j)
/
$include "integer_cut/1ATP_integer_cut.txt"
/
onConstraint(k)
;
**For adding integer cut from previous run
onConstraint(k)$previousk(k) = 1;
$offtext
*///---------------------------------------------------------

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

Nint(l, j)
/
$include "%dbpath%null_sij_nocofactor_20161025.txt"
/

y(j)
;

y(j) = 0;

Scalar
M /1000/
flag /0/
count /0/
epsilon /0.5/

;

Variables
*** objective function
z
G(j)
;

Integer Variables
*** flux
v(j)
vf(j)
vb(j)
;

Binary variables
yf(j)
yb(j)
a(j)
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

LB(j)$blocked(j) = 0;
UB(j)$blocked(j) = 0;

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

LB('EX_adp') = -nATP;
UB('EX_adp') = -nATP;

LB('EX_phosphate') = -nATP;
UB('EX_phosphate') = -nATP;

LB('EX_atp') = nATP;
UB('EX_atp') = nATP;

LB('EX_h2o') = nATP;
UB('EX_h2o') = nATP;

LB('EX_h+') = -10;
UB('EX_h+') = 10;

v.lo(j) = LB(j);
v.up(j) = UB(j);

*Fixing the bound for vf and vb
vf.lo(j) = 0;
vf.up(j) = M;

vb.lo(j) = 0;
vb.up(j) = M;

G.lo(j) = -M;
G.up(j) = M;

*Turn off all undesirable reactions
v.fx(j)$blocked(j) = 0;
yf.fx(j)$blocked(j) = 0;
yb.fx(j)$blocked(j) = 0;
vf.fx(j)$blocked(j) = 0;
vb.fx(j)$blocked(j) = 0;

*Irreversible forward
vb.fx(j)$(rxntype(j) = 0) = 0;
yb.fx(j)$(rxntype(j) = 0) = 0;

*Irreversible backward
vf.fx(j)$(rxntype(j) = 2) = 0;
yf.fx(j)$(rxntype(j) = 2) = 0;

***----------------------------------------------------------

Equations
obj
objmin
flux
stoic
cons1
cons2
cons3
cons4
cons5
integercut
loopless
consl1
consl2
consl3
consl4
nadphcons1
nadphcons2
nadphcons3
nadphcons4
;


obj.. z =e= sum(j$(rxntype(j) ne 4), vf(j) + vb(j));
objmin.. z =e= zlb;
flux(j).. v(j) =e= vf(j) - vb(j);
stoic(i).. sum(j, S(i,j)*v(j)) =e= 0;
cons1(j).. vf(j) =g= yf(j) * epsilon;
cons2(j).. vf(j) =l= yf(j) * M;
cons3(j).. vb(j) =g= yb(j) * epsilon;
cons4(j).. vb(j) =l= yb(j) * M;
cons5(j).. yf(j) + yb(j) =l= 1;

***Turn on integer cut constraint only when the solution is found for that iteration
integercut(k)$(onConstraint(k)=1).. sum(j$(store(k,j)=1), 1-yf(j)-yb(j)) =g= 1;

***Alternate integer cut
*integercut(k)$(onConstraint(k)=1).. sum(j$(store(k,j)=1 and (rxntype(j) ne 4)), 1-yf(j)-yb(j)) + sum(j$(store(k,j)=0 and (rxntype(j) ne 4)), yf(j) + yb(j)) =g= 1;

***Constraint to force z to be at least a certain number
*cons3.. z =g= 10;

***Loopless constraints
loopless(l).. sum(j$jloop(j), Nint(l,j) * G(j)) =e= 0;
consl1(j)$jloop(j).. G(j) =g= -M*a(j) + (1-a(j));
consl2(j)$jloop(j).. G(j) =l= -a(j) + M*(1-a(j));
consl3(j)$jloop(j).. v(j) =g= -M *(1-a(j));
consl4(j)$jloop(j).. v(j) =l= M*a(j);

***Fix nad(p)h production and consumption
nadphcons1.. v('EX_nadph') + v('EX_nadh') =e= 2;
nadphcons2.. v('EX_nadp') + v('EX_nad') =e= -2;
nadphcons3.. v('EX_nadh') + v('EX_nad') =e= 0;
nadphcons4.. v('EX_nadph') + v('EX_nadp') =e= 0;


Model
findpath
/
obj
objmin
flux
stoic
cons1
cons2
cons3
cons4
cons5
integercut
loopless
consl1
consl2
consl3
consl4
nadphcons1
nadphcons2
nadphcons3
nadphcons4
/
;
findpath.optfile = 1;
findpath.holdfixed = 1;

scalar xcount /0/;
***---------------------------------------------------------
put outputfile;
put "{"/;

*
*

while(iter <= maxiter and flag = 0,
     count = 0;
     xcount = 0;

     SOLVE findpath USING MIP MINIMIZING z;

     if ((findpath.modelstat ne 1 and findpath.modelstat ne 8),
          put @4, '"', @5, iter:0:0, '" : {'/;
          put @8, '"comment" : "No feasible solution found. Current iteration terminated with modelstat: ', findpath.modelstat:0:0'"'/;
          put @4, '}' /;
          flag = 1;
     else
         put @4, '"', @5, iter:0:0, '" : {'/;
         put @8, '"modelstat" : ', findpath.modelstat:0:0, ','/;
         put @8, '"solvestat" : ', findpath.solvestat:0:0, ','/;
         put @8, '"pathway" : {'/;

         loop(j$(v.l(j) ne 0),
              if((count ne 0),
                put ','/;
              );
              store(k,j)$(ord(k) = iter) = 1;
              put @12, '"', j.tl:0:30,'" : ', v.l(j):0:8;
              count = count + 1;
         );
         put /@8, '}, '/;

         put @8, '"pathway_y" : {'/;

         y(j) = yf.l(j) + yb.l(j);

         loop(j$(y(j) ne 0),
              if((xcount ne 0),
                put ','/;
              );
              put @12, '"', j.tl:0:30,'" : [', v.l(j):0:8, ',', vf.l(j):0:8, ',', vb.l(j):0:8, ']';
              xcount = xcount + 1;
         );
         put /@8, '}, '/;
         put @8, '"num_y" : ', xcount:0:0, ','/;

         put @8, '"num_reaction" : ', count:0:0, ','/;
         put @8, '"total_flux_no_exchange" : ', z.l:0:2, ','/;
         put @8, '"time" : ', findpath.resUsd:0:3 /;
         put @4, '},' /;
         flag = 0;
         onConstraint(k)$(ord(k) = iter) = 1;
         iter = iter + 1;
*store integer cut
         Execute_Unload '%outputfname%.gdx', store, iter, onConstraint;
     );
);
put "}"/;
putclose outputfile;