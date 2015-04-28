#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define p printf
#define q return
#define t struct
#define o void
#define M malloc
#define D free
#define xi extern u
#define xt extern
#define ls else
#define f if
#define w while
#define y for
#define K break
#define z sizeof
#define db double
#define tp typedef
tp const char cc;
cc* demo="duck.println('hello world')"
"\nfor i=1 to 20 do \n if i/2*2 != i then \n stars=\"\" \n for j=1 to i do \n"
" stars=stars+\"*\" \n loop \n duck.println(stars) \n end \nloop \n";
cc* es="Error\n";
tp int u;
tp unsigned char u8;tp unsigned long u32;
tp unsigned int uint;
#define SBS 4096
#define RRHS_SIZE 24
#define KV (256*16-1)
#define KS (256*16)
#define KT (256*32)
#define AE 0
#define SH 1
#define RD 2
#define AA 3
#define TL 1
#define TST 2
#define DOUBLE_EPSILON (1.0E-12)
#define PAGE_SIZE 4096
#define HTMC 2048
#define HTRF 8
#define GC_COLLECT_LIMIT 8192
enum {NIL=0, VP, VS, VR, VF, FP, DT};tp t R{u lhs;u rhsL;u* rhs;}R;tp t GT{u nS;char* sy;u nt;char* tokens;u nR;R* rules;}GT;tp t LI{u pr;u dot;u lookahead;}LI;tp t SET
{LI item;t SET* x;}SET;tp t IC{SET* set;t IC*x;}IC;tp t L{u tk;char* sr;u length;u line;t L* x;}L;tp t N{u ty;u v;}N;tp t LR{u nt;u nS;u numStates;N* actionTable;
u* gotoTable;}LR;tp t T{o* tk;u ty;u stt;t T* x;t T* prev;}T;t AST_TAG;tp t S{u tk;u pr;cc* sr;u length;u line;t AST_TAG *tag;u numC;t S** c;}S;
o pg(GT tb);o lab(L** str,L** h,L* x);L* eof(u line);L* newl(u line);L* IT(char* src,u str,u h,u line);L* ST(char* src,u str,u h,u line);
L* TO(char* src,u str,u h,u line,GT g);L* ID(char* src,u str,u h,u line,GT g);u isAlpha(char c);u iN(char c);u iS(char c);u iG(char c);L* LexSB(cc* sb,char** sp,GT tb);o FL(L* ex,char* B);
N At(LR tb,u stt,u tk);S* PSc(L* in,LR pr,GT g);u PS(L* in,LR pr,GT g);o FPT(S* syntaxTree);o PrintParseTree(S* ast,
GT g);o PrintParseTreeFormat(S* ast,GT g);u IN(S* no);o PrintNode(S* no);xt GT CFG;u* TAG;N* ACT;t C;t P;tp t F{
t P* pm;S* b;uint bi;u (*functor)(u);t C* cl;u rc;cc* fm;}F;t H;tp t V{u ty;union{u m;db fp;cc* sr;F* fn;t C* ref;t H* dt;}d;u const_string;}V;tp t P
{cc* id;V v;t P* x;}P;tp t C{P* list;u rc;t C* par;}C;tp t CS{cc* fm;F* fn;t CS* x;}CS;db TF(V v);u TI(V v);u IsD(V ty);V Cop(V sr);
o FE();V LN(cc* id);V CrF(u (*fn)(u));o AP(V functor,cc* argn);o LinkFn(V ref_nm,cc* id,V fn);o LCP(V ref_nm,cc* id,u v);o LCF(V ref_nm,cc* id,db v);o LCS(V ref_nm,
cc* id,cc* sr);xt C* gb;xt C* gCC;xt V A;xt P* gpl;xt t H* gDI;xi arr;xt P* ar;xt CS gStackTrace;tp t SK{C* cx;t SK* prev;t SK* x;}SK;
xt SK* Ex;o Push(C* cx);u Pop();xt cc* gv;xt C* gC;xt V ix;xt t H* gl;xi ai;xi ret;xi br;xi ct;xi hl;xi skd;xi gsd;xi le;xt S* pfp;xi gccc;xi gcic;o PCx(C* cx);
V get(cc* id,C* cx);o store(cc* id,V v,C* cx);db TF(V v);u TI(V v);u Inter(S*);o RdPAST(S**);o CCS(CS*);o PCS(CS*,F*);o PST();xt cc* ErrorMessage(u e);tp t MP
{o* ptr[PAGE_SIZE];uint page_n;t MP* np;}MP;xt MP* gBM;xt MP* gWM;o FE();V get(cc* id,C* cx);o store(cc* id,V v,C* cx);o RR(cc* id,C* cx);
tp t KVP{V key;V v;}KVP;tp t H{u cap;u sz;KVP* tb;u rc;}H;V HG(V kd,H* tb);o hs(V kd,V store,H* tb);uint HashFn(V v);H* cht();H* chtN(u sz);o fh(H* tb);
o ResizeHashTable(H* tb);u ET(V v);V TV();V Fal();V Con(V a,V b);V Add(V a,V b);V Sub(V a,V b);V Mu(V a,V b);V Div(V a,V b);V Mod(V a,V b);V CE(V a,V b);V CI(V a,V b);V LcXOR(V a,V b);V CLT(V a,V b);
V CGT(V a,V b);V LTE(V a,V b);V GTE(V a,V b);tp t LS{L* tokens;char* B;t LS* x;}LS;tp t AST_{S* pT;t AST_* x;}AST_;tp t GCD{LS* lexings;AST_* pT;C** ctx;
uint ctx_cap;uint ctx_size;F** fns;uint func_cap;uint func_size;char** strings;uint str_cap;uint str_size;H** tabs;uint tab_cap;uint tab_size;}GCD;GCD InitGC(GCD);
o FGCMO(GCD);u CGCRC(C*,GCD*);u CGCRD(H*,GCD*);u CGCRF(F*,GCD*);o CallGCTraceRoot(C*,V);o ClearAllGC();u GCAd(V,GCD*);
u GCAdS(char*,GCD*);u GCAdC(C*,GCD*);u GCAddFn(F*,GCD*);u GCAdD(H*,GCD*);o GCAdL(L* ex,char* B);o GCAdP(S* syntax);
xt GCD gGC;o ClearFn(F*);o CSt(char*);o CCx(C*);o CD(H*);o PrintValue(V);o PrintFn(F*);o PrintObject(C*);o PD(t H*);o BSL();
o BindStringLibrary();u Eval(u ac);u Parses(u ac);u DP(u ac);u DPt(u ac);u Type(u ac);u DInt(u ac);u DFl(u ac);u Len(u ac);u DQ(u ac);
u ET(V v){q ((v.ty==VP&&v.d.m)||(v.ty==FP&&fabs(v.d.fp)>DOUBLE_EPSILON)||(v.ty==VS||v.ty==VR||v.ty==VF||v.ty==DT));}
V TV(){V truth;truth.ty=VP;truth.d.m=1;q truth;}V Fal(){V truth;truth.ty=VP;truth.d.m=0;q truth;}
V Con(V a,V b){V v;char* string_a=NULL;char* string_b=NULL;char buffer1[32];char buffer2[32];u length;char* new_string;f(a.ty==VS){
string_a=(char*)a.d.sr;}ls{f(a.ty==VP){sprintf(buffer1,"%i",a.d.m);string_a=buffer1;}ls f(a.ty==FP){sprintf(buffer1,"%.16g",a.d.fp);string_a=buffer1;}
ls{sprintf(buffer1,"null");string_a=buffer1;}}f(b.ty==VS){string_b=(char*)b.d.sr;}ls{f(b.ty==VP){sprintf(buffer2,"%i",b.d.m);string_b=buffer2;}ls f(b.ty==FP){
sprintf(buffer2,"%.16g",b.d.fp);string_b=buffer2;}ls{sprintf(buffer2,"null");string_b=buffer2;}}
length=strlen(string_a)+strlen(string_b)+1;new_string=(char*)M(z(char)*length);GCAdS(new_string,&gGC);sprintf(new_string,"%s%s",string_a,string_b);v.ty=VS;v.d.sr=new_string;v.const_string=0;q v;}
V Add(V a,V b){V r;r.ty=NIL;r.d.m=0;f(a.ty==FP||b.ty==FP){r.ty=FP;r.d.fp=TF(a)+TF(b);}ls f(a.ty==VP||b.ty==VP){r.ty=VP;r.d.m=TI(a)+TI(b);}q r;}V Sub(V a,V b){V r;r.ty=NIL;r.d.m=0;f(a.ty==FP||b.ty==FP)
{r.ty=FP;r.d.fp=TF(a)-TF(b);}ls f(a.ty==VP||b.ty==VP){r.ty=VP;r.d.m=TI(a)-TI(b);}q r;}V Mu(V a,V b){V r;r.ty=NIL;r.d.m=0;f(a.ty==FP||b.ty==FP){r.ty=FP;r.d.fp=TF(a)*TF(b);}ls f(a.ty==VP||
b.ty==VP){r.ty=VP;r.d.m=TI(a)*TI(b);}q r;}V Div(V a,V b){V r;r.ty=NIL;r.d.m=0;f(a.ty==FP||b.ty==FP){db divisor=TF(b);r.ty=FP;f(divisor != 0.0){r.d.fp=TF(a) / divisor;}ls{r.ty=NIL;r.d.m=0;}}ls f(a.ty==VP||b.ty==VP){u divisor=TI(b);
r.ty=VP;f(divisor != 0){r.d.m=TI(a) / divisor;}ls{r.ty=NIL;r.d.m=0;}}q r;}V Mod(V a,V b){V r;u base;r.ty=NIL;r.d.m=0;base=TI(b);f(base != 0){r.ty=VP;r.d.m=TI(a) % base;}q r;}V CE(V a,V b){
V v;f((a.ty==FP&&b.ty==VP)||(a.ty==VP&&b.ty==FP)){float a_val,b_val;a_val=TF(a);b_val=TF(b);a.ty=b.ty=FP;a.d.fp=a_val;b.d.fp=b_val;}f(a.ty==b.ty){v.ty=VP;switch (a.ty){case VP:v.d.m=(a.d.m==b.d.m);K;
case FP:v.d.m=(fabs(a.d.fp-b.d.fp) < DOUBLE_EPSILON);K;case VS:v.d.m=strcmp(a.d.sr,b.d.sr)?0:1;K;case VR:v.d.m=(a.d.ref==b.d.ref);K;case VF:v.d.m=(a.d.fn==b.d.fn);
K;case DT:v.d.m=(a.d.dt==b.d.dt);K;}}ls{v.ty=VP;v.d.m=0;}q v;}V CI(V a,V b){V r=CE(a,b);r.ty=VP;r.d.m=(r.d.m?0:1);q r;}V LcAND(V a,V b){f(ET(a)&&ET(b)) q b;ls{V value_false;value_false.ty=VP;value_false.d.m=0;q value_false;}
}V LcOR(V a,V b){f(ET(a)) q a;ls f(ET(b)) q b;ls{V value_false;value_false.ty=VP;value_false.d.m=0;q value_false;}}V LcNOR(V a,V b){f(ET(a)){V value_false;value_false.ty=VP;value_false.d.m=0;q value_false;}
ls f(ET(b)){V value_false;value_false.ty=VP;value_false.d.m=0;q value_false;}ls{V value_true;value_true.ty=VP;value_true.d.m=1;q value_true;}}
V LcXOR(V a,V b){f(ET(a)&&!ET(b)){q a;}ls f(!ET(a)&&ET(b)){q b;}ls{V value_false;value_false.ty=VP;value_false.d.m=0;q value_false;}}V CLT(V a,V b)
{V r;f(a.ty==VP&&b.ty==VP){r.ty=VP;r.d.m=(a.d.m < b.d.m);}ls f(a.ty==FP||b.ty==FP){db afloat=TF(a);db bfloat=TF(b);
r.ty=VP;r.d.m=(afloat < bfloat);}ls{f(a.ty==VS&&b.ty==VS){r.ty=VP;r.d.m=(strcmp(a.d.sr,b.d.sr) < 0);}ls{r.ty=NIL;r.d.m=0;}}
r.d.m=r.d.m?1:0;q r;}V CGT(V a,V b){V r;f(a.ty==VP&&b.ty==VP){r.ty=VP;r.d.m=(a.d.m > b.d.m);}ls f(a.ty==FP||b.ty==FP){db afloat=TF(a);db bfloat=TF(b);r.ty=VP;r.d.m=(afloat > bfloat);}
ls{f(a.ty==VS&&b.ty==VS){r.ty=VP;r.d.m=(strcmp(a.d.sr,b.d.sr) > 0);}ls{r.ty=NIL;r.d.m=0;}}r.d.m=r.d.m?1:0;q r;}
V LTE(V a,V b){V r;f(a.ty==VP&&b.ty==VP){r.ty=VP;r.d.m=(a.d.m <= b.d.m);}ls f(a.ty==FP||b.ty==FP){db afloat=TF(a);db bfloat=TF(b);r.ty=VP;r.d.m=(afloat <= bfloat);}
ls{f(a.ty==VS&&b.ty==VS){r.ty=VP;r.d.m=(strcmp(a.d.sr,b.d.sr) <= 0);}ls{r.ty=NIL;r.d.m=0;}}r.d.m=r.d.m?1:0;q r;}V GTE(V a,V b){V r;f(a.ty==VP&&b.ty==VP){
r.ty=VP;r.d.m=(a.d.m >= b.d.m);}ls f(a.ty==FP||b.ty==FP){db afloat=TF(a);db bfloat=TF(b);r.ty=VP;r.d.m=(afloat >= bfloat);}ls{f(a.ty==VS&&b.ty==VS){r.ty=VP;r.d.m=(strcmp(a.d.sr,b.d.sr) >= 0);}ls{
r.ty=NIL;r.d.m=0;}}r.d.m=r.d.m?1:0;q r;}u RdProgram(S* no){S* sl1=no->c[0];u e=0;e=IN(sl1);q e;}u RdStListA(S* no){S* stmt1=no->c[0];S* sl1=no->c[1];u e=0;f(hl==0){
w (no->numC){f(ret==0&&br==0&&ct==0&&e==0&&hl==0){e=IN(no->c[0]);}ls{f(hl){ret=1;br=1;}q e;}f(e==0&&ret==0){f((++gcic)==GC_COLLECT_LIMIT){gcic=0;gccc++;CallGCTraceRoot(gb,A);}
}no=no->c[1];}}ls{ret=1;br=1;}q e;}u RdStA(S* no){u e=0;V library=get(no->c[1]->sr,gCC);f(library.ty==VR){C* nm=library.d.ref;P* itr=nm->list;w (itr){store(itr->id,itr->v,gb);itr=itr->x;}
}q e;}u RdStJ(S* no){S* r=no->c[1];u e=0;e=IN(r);ret=1;q e;}u RdStK(S* no){u e=0;br=1;q e;}u RdStL(S* no){u e=0;ct=1;q e;}u RdFnDef(S* no){V rr;S* id1=no->c[1];S* parameters1=no->c[2];S* sl1=no->c[4];
u e=0;gpl=NULL;e=IN(parameters1);f(e==0){rr.ty=VF;rr.d.fn=(F*)M(z(F));rr.d.fn->pm=gpl;rr.d.fn->b=sl1;rr.d.fn->cl=gCC;rr.d.fn->bi=0;rr.d.fn->functor=NULL;rr.d.fn->fm=id1->sr;store(id1->sr,rr,gCC);
GCAddFn(rr.d.fn,&gGC);}q e;}u RdPmDeclA(S* no){S* id1=no->c[0];u e=0;gpl=(P*)M(z(P));gpl->id=id1->sr;gpl->v.ty=NIL;gpl->x=NULL;q e;}u RdPmDeclB(S* no){S* param_decl1=no->c[0];S* id1=no->c[2];P* last;u e=0;
e=IN(param_decl1);last=gpl;w (last->x){last=last->x;}last->x=(P*)M(z(P));last->x->id=id1->sr;last->x->v.ty=NIL;last->x->x=NULL;q e;}u RdIf(S* no){S* cd1=no->c[1];S* sl1=no->c[4];S* else_if1=no->c[5];
u e=0;e=IN(cd1);f(e==0){f(ET(A)){e=IN(sl1);}ls{e=IN(else_if1);}}q e;}u RdForLoop(S* no){S* id1=no->c[1];S* ar1=no->c[3];S* ar2=no->c[5];S* sl1=no->c[8];
V str,h;cc* id=id1->sr;u e=0;str.ty=NIL;h.ty=NIL;e=IN(ar1);str=A;f(e==0){e=IN(ar2);h=A;}f(str.ty==FP&&h.ty==FP){w (str.ty==FP&&h.ty==FP&&str.d.fp <= h.d.fp&&e==0&&br==0)
{store(id,str,gCC);e=IN(sl1);f(br==0){ct=0;str.d.fp=str.d.fp+1.0;}}}ls f(str.ty==VP&&h.ty==VP){w (str.ty==VP&&h.ty==VP&&str.d.m <= h.d.m&&e==0&&br==0){store(id,str,gCC);e=IN(sl1);f(br==0)
{ct=0;str.d.m++;}}}ls{e=25;br=0;q e;}br=0;q e;}u RdWhileLoop(S* no){S* cd1=no->c[1];S* sl1=no->c[4];u e=0;e=IN(cd1);w (A.ty != NIL&&ET(A)&&e==0&&br==0){
e=IN(sl1);f(e==0){e=IN(cd1);}ct=0;}br=0;q e;}u RdAsA(S* no){S* l_value1=no->c[0];S* assignment1=no->c[2];V xp;u e=0;e=IN(assignment1);xp=A;f(e) q e;e=IN(l_value1);f(e) q e;f(ai)
{V previous=HG(ix,gl);hs(ix,xp,gl);ai=0;}ls{V previous=get(gv,gC);f(ix.ty==VS&&ix.const_string==0){ix=Cop(ix);}store(gv,xp,gC);}A=xp;q e;}u RdAsB(S* no)
{S* l_value1=no->c[0];S* cd1=no->c[2];V xp;u e=0;e=IN(cd1);xp=A;f(e) q e;e=IN(l_value1);f(e) q e;f(ai){V previous=HG(ix,gl);f(ix.ty==VS&&ix.const_string==0){ix=Cop(ix);}hs(ix,xp,gl);ai=0;}ls{
V previous=get(gv,gC);store(gv,xp,gC);}A=xp;q e;}u RdLValueA(S* no){S* id1=no->c[0];u e=0;ai=0;gv=id1->sr;gC=gCC;ix.ty=NIL;ix.d.m=0;gl=NULL;q e;}u RdLValueC(S* no){
S* ref1=no->c[0];S* id1=no->c[2];u e=0;e=IN(ref1);f(e) q e;f(A.ty==VR){ai=0;gv=id1->sr;gC=A.d.ref;ix.ty=NIL;ix.d.m=0;gl=NULL;}ls f(A.ty==DT){
ai=1;ix.ty=VS;ix.d.sr=id1->sr;ix.const_string=1;gl=A.d.dt;gv=NULL;gC=NULL;}ls{e=31;ai=0;gv=NULL;gC=NULL;ix.ty=NIL;ix.d.m=0;gl=NULL;pfp=no;}q e;}u RdLValueD(S* no){
S* ref1=no->c[0];S* r=no->c[2];V ref;V n;u e=0;e=IN(ref1);ref=A;f(e) q e;e=IN(r);n=A;f(e) q e;f(ref.ty==VR&&(n.ty==VS||n.ty==VP)){f(n.ty==VP||n.ty==FP)
{V ns;ns.ty=VS;ns.d.sr="";n=Con(n,ns);}gv=n.d.sr;gC=ref.d.ref;ix.ty=NIL;ix.d.m=0;gl=NULL;ai=0;}ls f(ref.ty==DT){ix=n;gl=ref.d.dt;gv=NULL;gC=NULL;ai=1;}
ls{e=32;gv=NULL;gC=NULL;ix.ty=NIL;ix.d.m=0;gl=NULL;ai=0;pfp=no;}q e;}u RdCnA(S* no){S* cd1=no->c[0];S* logic1=no->c[2];V cd;u e=0;e=IN(cd1);cd=A;f(e) q e;f(ET(cd)){e=IN(logic1);f(e) q e;}ls{A.ty=VP;A.d.m=0;}
q e;}u RdCnB(S* no){S* cd1=no->c[0];S* logic1=no->c[2];V cd;V logic;u e=0;e=IN(cd1);cd=A;f(e) q e;f(ET(cd)){A=cd;q e;}e=IN(logic1);logic=A;f(e) q e;A=logic;q e;}u RdCnC(S* no)
{S* cd1=no->c[0];S* logic1=no->c[2];V cd;V logic;u e=0;e=IN(cd1);cd=A;f(e) q e;f(ET(cd)){A.ty=VP;A.d.m=0;q e;}e=IN(logic1);logic=A;f(e) q e;f(ET(logic)){A.ty=VP;A.d.m=0;q e;}ls{A.ty=VP;A.d.m=1;q e;}q e;}
u RdCnD(S* no){S* cd1=no->c[0];S* logic1=no->c[2];V cd;V logic;u e=0;e=IN(cd1);cd=A;f(e) q e;e=IN(logic1);logic=A;f(e) q e;A=LcXOR(cd,logic);q e;}u RdLcA(S* no){S* cm1=no->c[1];u e=0;e=IN(cm1);
f(A.ty==VP){A.d.m=A.d.m?0:1;}ls f(A.ty==NIL){A.ty=VP;A.d.m=1;}ls{A.ty=NIL;A.d.m=0;}q e;}u RdCmA(S* no){S* cm1=no->c[0];S* ar1=no->c[2];V cm;V ar;u e=0;e=IN(cm1);cm=A;f(e) q e;e=IN(ar1);ar=A;f(e) q e;
A=CE(cm,ar);q e;}u RdCmB(S* no){S* cm1=no->c[0];S* ar1=no->c[2];V cm;V ar;u e=0;e=IN(cm1);cm=A;f(e) q e;e=IN(ar1);ar=A;f(e) q e;A=CI(cm,ar);q e;}u RdCmC(S* no){S* cm1=no->c[0];S* ar1=no->c[2];V cm;V ar;
u e=0;e=IN(cm1);cm=A;f(e) q e;e=IN(ar1);ar=A;f(e) q e;A=CLT(cm,ar);q e;}u RdCmD(S* no){S* cm1=no->c[0];S* ar1=no->c[2];V cm;V ar;u e=0;e=IN(cm1);cm=A;f(e) q e;e=IN(ar1);ar=A;f(e) q e;
A=CGT(cm,ar);q e;}u RdCmE(S* no){S* cm1=no->c[0];S* ar1=no->c[2];V cm;V ar;u e=0;e=IN(cm1);cm=A;f(e) q e;e=IN(ar1);ar=A;f(e) q e;A=LTE(cm,ar);q e;}u RdCmF(S* no){S* cm1=no->c[0];S* ar1=no->c[2];V ar;V cm;
u e=0;e=IN(cm1);cm=A;f(e) q e;e=IN(ar1);ar=A;f(e) q e;A=GTE(cm,ar);q e;}u RdArithmeticA(S* no){S* ar1=no->c[0];S* term1=no->c[2];V ar;V term;u e=0;e=IN(ar1);ar=A;A.ty=NIL;f(e) q e;e=IN(term1);term=A;f(e) q e;
f(ar.ty==VS||term.ty==VS){A=Con(ar,term);}ls{A=Add(ar,term);}q e;}u RdAB(S* no){S* ar1=no->c[0];S* term1=no->c[2];V ar;V term;u e=0;e=IN(ar1);ar=A;f(e) q e;e=IN(term1);term=A;f(e) q e;
A=Sub(ar,term);q e;}u RdArithmeticC(S* no){S* ar1=no->c[0];S* term1=no->c[2];V ar;V term;u e=0;e=IN(ar1);ar=A;f(e) q e;e=IN(term1);term=A;f(e) q e;f(ar.ty==VP&&term.ty==VP){A.ty=VP;A.d.m =
(ar.d.m & term.d.m);}ls{A.ty=NIL;}q e;}u RdArithmeticD(S* no){S* ar1=no->c[0];S* term1=no->c[2];V ar;V term;u e=0;e=IN(ar1);ar=A;f(e) q e;e=IN(term1);term=A;f(e) q e;f(ar.ty==VP&&term.ty==VP){A.ty=VP;A.d.m =
(ar.d.m | term.d.m);}ls{A.ty=NIL;}q e;}u RdTermA(S* no){S* term1=no->c[0];S* factor1=no->c[2];V term;V factor;u e=0;e=IN(term1);term=A;f(e) q e;e=IN(factor1);factor=A;f(e) q e;
A=Mu(term,factor);q e;}u RdTermB(S* no){S* term1=no->c[0];S* factor1=no->c[2];V term;V factor;u e=0;e=IN(term1);term=A;f(e) q e;e=IN(factor1);factor=A;f(e) q e;
A=Div(term,factor);q e;}u RdFactorA(S* no){S* factor1=no->c[1];u e=0;e=IN(factor1);f(A.ty==VP){A.d.m=-A.d.m;}ls f(A.ty==FP){A.d.fp=-A.d.fp;}ls{A.ty=NIL;}q e;}u RdFB(S* no){
S* factor1=no->c[1];u e=0;e=IN(factor1);f(A.ty==VP){A.d.m=!A.d.m;}ls f(A.ty==NIL){A.ty=VP;A.d.m=1;}ls{A.ty=NIL;}q e;}u RdFlC(S* no){S* integer1=no->c[0];u e=0;A.ty=VP;A.d.m=atoi(integer1->sr);q e;}
u RdFlD(S* no){S* float1=no->c[0];u e=0;A.ty=FP;A.d.fp=atof(float1->sr);q e;}u RdFlE(S* no){S* string1=no->c[0];u e=0;A.ty=VS;A.d.sr=string1->sr;A.const_string=1;q e;}u RdRefA(S* no){S* l_value1=no->c[0];u e=0;e=IN(l_value1);
f(ai){A=HG(ix,gl);ai=0;}ls{A=get(gv,gC);}q e;}u RdRefB(S* no){S* ref1=no->c[0];V fn;u e=0;e=IN(ref1);fn=A;f(e) q e;f(fn.ty==VF){C* fc=(C*)M(z(C));C* cr=gCC;
fc->list=NULL;fc->par=fn.d.fn->cl;gCC=fc;Push(cr);GCAdC(fc,&gGC);f(fn.d.fn->bi){e=fn.d.fn->functor(0);}ls{e=IN(fn.d.fn->b);}f(e){PCS(&gStackTrace,fn.d.fn);}ls{}Pop();gCC=cr;ret=0;}
ls{e=12345;pfp=no;}q e;}u RdRefC(S* no){S* ref1=no->c[0];S* a1=no->c[2];V fn;u ac;u e=0;e=IN(ref1);fn=A;f(e) q e;f(fn.ty==VF){P* param;P* arg;C* fc;C* cr;P* last;P* an;
e=IN(a1);param=fn.d.fn->pm;arg=ar;f(e) q e;fc=(C*)M(z(C));fc->list=NULL;fc->par=fn.d.fn->cl;cr=gCC;gCC=fc;Push(cr);GCAdC(fc,&gGC);last=NULL;ac=0;w (param&&arg){A.ty=NIL;f(last){
last->x=(P*)M(z(P));last=last->x;}ls{gCC->list=last=(P*)M(z(P));}last->id=param->id;last->v=arg->v;last->x=NULL;ac++;an=arg->x;D(arg);param=param->x;arg=an;}f(fn.d.fn->bi){
e=fn.d.fn->functor(0);}ls{e=IN(fn.d.fn->b);}f(e){PCS(&gStackTrace,fn.d.fn);}ls{}Pop();gCC=cr;ret=0;}ls{e=12345;pfp=no;}q e;}u RdAA(S* no){S* a1=no->c[0];S* r=no->c[2];P* as;P* list;u e=0;e=IN(a1);
as=ar;e=e?e:IN(r);list=as;w (list->x){list=list->x;}list->x=(P*)M(z(P));list->x->id=NULL;list->x->v=A;list->x->x=NULL;ar=as;q e;}u RdArgumentsB(S* no){S* r=no->c[0];u e=0;e=IN(r);
ar=(P*)M(z(P));ar->id=NULL;ar->v=A;ar->x=NULL;q e;}u RdObjectA(S* no){u e=0;A.ty=DT;A.d.dt=cht();GCAdD(A.d.dt,&gGC);A.d.dt->rc=1;q e;}u RdObjectB(S* no){S* array_init1=no->c[1];u e=0;e=IN(array_init1);
A.ty=DT;A.d.dt=gDI;GCAdD(gDI,&gGC);A.d.dt->rc=1;q e;}u RdObjectC(S* no){S* dic=no->c[1];u e=0;e=IN(dic);A.ty=DT;A.d.dt=gDI;GCAdD(gDI,&gGC);A.d.dt->rc=1;q e;}u RdAIA(S* no){S* array_init1=no->c[0];S* r=no->c[2];H* dt;V xp;V key;u e=0;u n;
e=IN(array_init1);dt=gDI;n=arr;e=(e?e:IN(r));xp=A;key.ty=VP;key.d.m=n;hs(key,xp,dt);arr=++n;gDI=dt;q e;}u RdArrayInitB(S* no){S* r=no->c[0];V xp;u e=0;V n;e=IN(r);xp=A;gDI=cht();arr=0;n.ty=VP;n.d.m=arr;hs(n,xp,gDI);arr++;q e;}
u RdDIA(S* no){S* dic=no->c[0];S* id1=no->c[2];S* r=no->c[4];H* dt;u e=0;V xp;V key;V dupe;e=IN(dic);dt=gDI;e=e?e:IN(r);xp=A;key.ty=VS;key.d.sr=id1->sr;key.const_string=1;dupe=HG(key,dt);f(dupe.ty != NIL)
{e=76;q e;}hs(key,xp,dt);gDI=dt;q e;}u RdDI(S* no){S* id1=no->c[0];S* r=no->c[2];V xp;u e=0;V key;e=IN(r);xp=A;gDI=cht();key.ty=VS;key.d.sr=id1->sr;key.const_string=1;hs(key,xp,gDI);q e;}u RdBooleanA(S* no){u e=0;A.ty=VP;A.d.m=1;q e;}
u RdB(S* no){u e=0;A.ty=VP;A.d.m=0;q e;}u IN(S* no){f(no==NULL||no->pr==0){p(es);q 1;}le=no->line;pfp=no;switch (no->pr){
case 0x01:q RdProgram(no);case 0x02:q RdStListA(no);case 0x04:q RdStA(no);
case 0x0D:q RdStJ(no);case 0x0E:q RdStK(no);case 0x0F:q RdStL(no);
case 0x10:q RdFnDef(no);case 0x14:q RdPmDeclA(no);case 0x15:q RdPmDeclB(no);
case 0x16:q RdIf(no);case 0x19:q RdForLoop(no);case 0x1A:q RdWhileLoop(no);
case 0x1B:q RdAsA(no);case 0x1C:q RdAsB(no);case 0x1D:q RdLValueA(no);
case 0x1F:q RdLValueC(no);case 0x20:q RdLValueD(no);case 0x22:q RdCnA(no);
case 0x23:q RdCnB(no);case 0x24:q RdCnC(no);case 0x25:q RdCnD(no);
case 0x27:q RdLcA(no);case 0x29:q RdCmA(no);case 0x2A:q RdCmB(no);
case 0x2B:q RdCmC(no);case 0x2C:q RdCmD(no);case 0x2D:q RdCmE(no);
case 0x2E:q RdCmF(no);case 0x30:q RdArithmeticA(no);case 0x31:q RdAB(no);
case 0x32:q RdArithmeticC(no);case 0x33:q RdArithmeticD(no);case 0x35:q RdTermA(no);
case 0x36:q RdTermB(no);case 0x38:q RdFactorA(no);case 0x39:q RdFB(no);
case 0x3D:q RdFlC(no);case 0x3E:q RdFlD(no);case 0x3F:q RdFlE(no);
case 0x42:q RdRefA(no);case 0x43:q RdRefB(no);case 0x44:q RdRefC(no);
case 0x45:q RdAA(no);case 0x46:q RdArgumentsB(no);case 0x47:q RdObjectA(no);
case 0x48:q RdObjectB(no);case 0x49:q RdObjectC(no);case 0x4A:q RdAIA(no);
case 0x4B:q RdArrayInitB(no);case 0x4C:q RdDIA(no);case 0x4D:q RdDI(no);
case 0x4F:q RdBooleanA(no);case 0x50:q RdB(no);
case 0xFF:q 0;default:p(es);q 1;}}
GCD gGC;u APA(o*** array,uint* cap,uint* sz,o* entry)
{u r=0;f(entry){u i;y(i=0;i < *sz;i++)f((*array)[i]==entry)K;f(i==*sz){r=1;(*array)[*sz]=entry;(*sz)++;f(*sz==*cap){*cap=*cap*2;*array=realloc((o*)(*array),z(o*)**cap);
f(array==0) p("Error:Out of memory.\n");}}}q r;}u CPA(o** array,uint sz,o* entry){uint i;y(i=0;i < sz;i++){f(array[i]==entry)q 1;}q 0;}GCD
InitGC(GCD gcStore){gcStore.lexings=NULL;gcStore.pT=NULL;gcStore.ctx_size=0;gcStore.ctx_cap=128;gcStore.ctx=(C**)M(z(C*)*gcStore.ctx_cap);gcStore.func_size=0;gcStore.func_cap=128;gcStore.fns=(F**)M(z(F*)*gcStore.func_cap);gcStore.str_size=0;gcStore.str_cap=128;gcStore.strings=(char**)M(z(char*)*gcStore.str_cap);gcStore.tab_size=0;gcStore.tab_cap=128;gcStore.tabs=(H**)M(z(H*)*gcStore.tab_cap);q gcStore;}
o FGCMO(GCD gcStore){D(gcStore.ctx);D(gcStore.fns);D(gcStore.strings);D(gcStore.tabs);}u CGCRC(C* cx,GCD* mgmt){
P* itr;u r=0;w (cx){itr=cx->list;f(GCAdC(cx,mgmt)){r=1;w (itr){f(itr->v.ty==VR){CGCRC(itr->v.d.ref,
mgmt);}ls f(itr->v.ty==VF){CGCRF(itr->v.d.fn,mgmt);}ls f(itr->v.ty==DT){CGCRD(itr->v.d.dt,mgmt);}ls f(itr->v.ty==VS&&itr->v.const_string==0){GCAdS((char*)itr->v.d.sr,mgmt);}itr=itr->x;}}
cx=cx->par;}q r;}u CGCRD(H* tb,GCD* mgmt){u r=0;f(tb&&GCAdD(tb,mgmt)){uint n;r=1;y(n=0;n < tb->cap;n++)
{f(tb->tb[n].key.ty != NIL){V v=tb->tb[n].v;f(v.ty==VR){CGCRC(v.d.ref,mgmt);}ls f(v.ty==VF){CGCRF(v.d.fn,mgmt);}ls f(v.ty==DT){CGCRD(v.d.dt,
mgmt);}ls f(v.ty==VS&&v.const_string==0){GCAdS((char*)v.d.sr,mgmt);}}}}q r;}u CGCRF(F* fn,GCD* mgmt){u r=0;f(fn&&GCAddFn(fn,mgmt)){r=1;CGCRC(fn->cl,mgmt);}
q r;}o CallGCTraceRoot(C* root,V cE){GCD mgmt;SK* sk;uint n;mgmt=InitGC(mgmt);GCAd(cE,&mgmt);CGCRC(gCC,&mgmt);sk=Ex;w (sk){CGCRC(sk->cx,&mgmt);sk=sk->x;}
y(n=0;n < gGC.ctx_size;n++){C* cx=gGC.ctx[n];f(!CPA((o**)mgmt.ctx,mgmt.ctx_size,(o*)cx)){CCx(cx);gGC.ctx[n]=NULL;gGC.ctx[n]=gGC.ctx[gGC.ctx_size-1];gGC.ctx_size--;n--;}}y(n=0;n < gGC.func_size;n++){F* fn=gGC.fns[n];
f(!CPA((o**)mgmt.fns,mgmt.func_size,(o*)fn)){ClearFn(fn);gGC.fns[n]=NULL;gGC.fns[n]=gGC.fns[gGC.func_size-1];gGC.func_size--;n--;}}y(n=0;n < gGC.tab_size;n++){H* tb=gGC.tabs[n];f(!CPA((o**)mgmt.tabs,mgmt.tab_size,(o*)tb)){CD(tb);gGC.tabs[n]=NULL;
gGC.tabs[n]=gGC.tabs[gGC.tab_size-1];gGC.tab_size--;n--;}}y(n=0;n < gGC.str_size;n++){char* sr=gGC.strings[n];f(!CPA((o**)mgmt.strings,mgmt.str_size,(o*)sr)){CSt(sr);gGC.strings[n]=NULL;
gGC.strings[n]=gGC.strings[gGC.str_size-1];gGC.str_size--;n--;}}FGCMO(mgmt);}o ClearAllGC(){uint n;AST_ *cur_ast,*nast;LS *cur_lex,*next_lex;y(n=0;n<gGC.ctx_size;n++){C* cx=gGC.ctx[n];CCx(cx);}
y(n=0;n<gGC.func_size;n++){F* fn=gGC.fns[n];ClearFn(fn);}y(n=0;n<gGC.tab_size;n++){H* tb=gGC.tabs[n];CD(tb);}
y(n=0;n<gGC.str_size;n++){char* sr=gGC.strings[n];CSt(sr);}cur_ast=gGC.pT;w(cur_ast){
nast=cur_ast->x;FPT(cur_ast->pT);D(cur_ast);cur_ast=nast;}cur_lex=gGC.lexings;w(cur_lex){next_lex=cur_lex->x;FL(cur_lex->tokens,cur_lex->B);D(cur_lex);cur_lex=next_lex;}}
u GCAd(V v,GCD* gc){f(v.ty==VS&&v.const_string==0&&v.d.sr){q GCAdS((char*)v.d.sr,gc);}ls f(v.ty==VR&&v.d.ref){q GCAdC(v.d.ref,gc);}ls f(v.ty==VF&&
v.d.fn){q GCAddFn(v.d.fn,gc);}ls f(v.ty==DT&&v.d.dt){q GCAdD(v.d.dt,gc);}q 0;}u GCAdS(char* sr,GCD* gc){q APA((o***)&gc->strings,&(gc->str_cap),&(gc->str_size),(o*)sr);}
u GCAdC(C* cx,GCD* gc){q APA((o***)&gc->ctx,&(gc->ctx_cap),&(gc->ctx_size),(o*)cx);}u GCAddFn(F* fn,GCD* gc){q APA((o***)&gc->fns,&(gc->func_cap),&(gc->func_size),(o*)fn);}
u GCAdD(H* tb,GCD* gc){q APA((o***)&gc->tabs,&(gc->tab_cap),&(gc->tab_size),(o*)tb);}o GCAdL(L* ex,char* B){f(gGC.lexings==NULL){gGC.lexings=(LS*)M(z(LS));gGC.lexings->tokens=ex;gGC.lexings->B=B;
gGC.lexings->x=NULL;}ls{LS* cr=gGC.lexings;w (cr->x) cr=cr->x;cr->x=(LS*)M(z(LS));cr=cr->x;cr->tokens=ex;cr->B=B;cr->x=NULL;}}
o GCAdP(S* syntax){f(gGC.pT==NULL){gGC.pT=(AST_*)M(z(AST_));gGC.pT->pT=syntax;gGC.pT->x=NULL;}ls{AST_* cr=gGC.pT;w (cr->x) cr=cr->x;cr->x=(AST_*)M(z(AST_));cr=cr->x;cr->pT=syntax;
cr->x=NULL;}}o ClearFn(F* fn){f(fn){P *itr,*itr_next;itr=fn->pm;w (itr){itr_next=itr->x;D(itr);itr=itr_next;}D(fn);}}o CSt(char* sr){D(sr);}o CCx(C* cx){
f(cx){P *itr,*itr_next;itr=cx->list;w (itr){itr_next=itr->x;D(itr);itr=itr_next;}D(cx);}}o CD(H* tb){fh(tb);}u gSR=0;u gGoal=KS+1;u ep=KT+1;u gSymbolEOF=KT+2;u gEnd=KT+3;
u gInt=KT+4;u gSF=KT+5;u gSS=KT+6;u gSI=KT+7;char* GetE(u i,GT g){char* sc;char* h;u n;f(i & KS){sc=g.sy;h=g.sy+SBS;i ^= KS;}ls f(i & KT){
sc=g.tokens;h=g.tokens+SBS;i ^= KT;}ls{q 0;}n=1;w (n < i&&sc < h){w (sc[0] != 0&&sc != h)sc++;sc++;n++;}f(sc < h)q sc;q 0;}R* GetRule(u lhs,u n,GT g)
{u r;y(r=0;r < g.nR;r++){f(g.rules[r].lhs==lhs){f(n==0)q g.rules+r;n--;}}q NULL;}u AS(char* sr,u len,GT* g){
u j;char* sc=g->sy;char* h=g->sy+SBS;u n=1;w (n <= g->nS&&sc < h){y(j=0;j < len;j++){f(sr[j] != sc[j])K;}f(j==len&&sc[j]==0){q n | KS;}w (sc[0] != 0&&sc != h)
sc++;sc++;n++;}f(sc < h){y(j=0;j < len&&sc < h;j++){sc[0]=sr[j];sc++;}sc[0]=0;g->nS++;q n | KS;}q 0;}u AT(char* sr,u len,GT* g){
u j;char* sc=g->tokens;char* h=g->tokens+SBS;u n=1;w (n <= g->nt&&sc < h){y(j=0;j < len;j++){f(sr[j]=='\\') j++;f(sr[j] != sc[j])K;}f(j==len&&sc[j]==0)
{q n | KT;}w (sc[0] != 0&&sc != h)sc++;sc++;n++;}f(sc < h){y(j=0;j < len&&sc < h;j++){f(sr[j]=='\\') j++;sc[0]=sr[j];sc++;}sc[0]=0;g->nt++;q n | KT;}
q 0;}u FT(char* sr,u len,GT g){u j;char* sc=g.tokens;char* h=g.tokens+SBS;u n=1;w (n <= g.nt&&
sc < h){y(j=0;j < len;j++){f(sr[j]=='\\') j++;f(sr[j] != sc[j])K;}f(j==len&&sc[j]==0){q n | KT;}
w (sc[0] != 0&&sc != h)sc++;sc++;n++;}q 0;}o fg(GT* tb){u i;D(tb->sy);D(tb->tokens);tb->nS=0;tb->nt=0;y(i=0;i < tb->nR;i++)D(tb->rules[i].rhs);D(tb->rules);tb->nR=0;}
C* gb;C* gCC;V A;P* gpl;H* gDI;u arr;P* ar;cc* gv;C* gC;V ix;H* gl;u ai;u ret;u br;u ct;u hl;CS gStackTrace;u le;S* pfp;u gsd;u skd;u gccc;u gcic;o PCx(C* cx)
{w (cx){P* list=cx->list;p("{");w (list){V v=list->v;p("%s:",list->id);f(v.ty==VP)p("%i",v.d.m);ls f(v.ty==VS)p("%s",v.d.sr);ls f(v.ty==VR)
{p("[");PCx(v.d.ref);p("]");}ls f(v.ty==VF)p("f()");ls p("nil");p(",");list=list->x;}cx=cx->par;p("}");f(cx) p("\n");}}
SK* Ex=NULL;o Push(C* cx){f(Ex==NULL){Ex=(SK*)M(z(SK));Ex->cx=cx;Ex->x=NULL;Ex->prev=NULL;}ls{
SK* tail=Ex;w (tail->x){tail=tail->x;}tail->x=(SK*)M(z(SK));tail->x->cx=cx;tail->x->x=NULL;tail->x->prev=tail;}
}u Pop(){SK* tail=Ex;w (tail&&tail->x){tail=tail->x;}f(tail){gCC=tail->cx;tail=tail->prev;f(tail){D(tail->x);tail->x=NULL;}ls{D(Ex);Ex=NULL;}
q 1;}q 0;}V LN(cc* id){V ref_nm;ref_nm.ty=VR;ref_nm.d.ref=(C*)M(z(C));ref_nm.d.ref->par=NULL;ref_nm.d.ref->list=NULL;ref_nm.d.ref->rc=-1;store(id,ref_nm,gb);GCAdC(ref_nm.d.ref,&gGC);q ref_nm;}
o LinkFn(V ref_nm,cc* id,V fn){f(fn.ty==VF){fn.d.fn->fm=id;store(id,fn,ref_nm.d.ref);}}o LCP(V ref_nm,cc* id,u v){V cn;cn.ty=VP;cn.d.m=v;store(id,cn,ref_nm.d.ref);}
o LCF(V ref_nm,cc* id,db v){V cn;cn.ty=FP;cn.d.fp=v;store(id,cn,ref_nm.d.ref);}o LCS(V nm,cc* id,cc* sr){V cn;cn.ty=VS;cn.d.sr=sr;cn.const_string=1;store(id,cn,nm.d.ref);}
V CrF(u (*fn)(u)){V rr;rr.ty=VF;rr.d.fn=(F*)M(z(F));rr.d.fn->pm=NULL;rr.d.fn->b=NULL;rr.d.fn->cl=gCC;rr.d.fn->bi=1;rr.d.fn->functor=fn;rr.d.fn->rc=-1;rr.d.fn->fm="fn";GCAddFn(rr.d.fn,&gGC);q rr;}
o AP(V functor,cc* argn){P* itr;P* pm=(P*)M(z(P));pm->id=argn;pm->v.ty=NIL;pm->v.d.m=0;pm->x=NULL;f(functor.d.fn->pm){y(itr=functor.d.fn->pm;itr->x;itr=itr->x);itr->x=pm;q;}
functor.d.fn->pm=pm;}o CCS(CS* sk){sk->x=0l;sk->fm=NULL;}o PCS(CS* sk,F* func){f(sk&&func){cc* id=(func?func->fm:NULL);f(sk->fm==NULL){(*sk).fm=id;(*sk).fn=func;(*sk).x=0l;}ls{CS* x;x=(CS*)M(z(CS));
x->fm=sk->fm;x->fn=sk->fn;x->x=sk->x;(*sk).fm=id;(*sk).fn=func;sk->x=x;}}}o PST(){CS* sk=&gStackTrace;p("Program halted on line %i.\n",(le+1));f(pfp){
p("Failed production:");PrintParseTreeFormat(pfp,CFG);p("\n");}f(sk->fm){u depth=0;u i;p("Printing stack trace:\n\n");w(sk){y(i=0;i<depth;i++)p(" ");
PrintFn(sk->fn);p("\n");sk=sk->x;depth++;}}}u Inter(S* tree){u e; gGC=InitGC(gGC);gcic=0;gccc=0;gCC=gb=(C*)M(z(C));gCC->par=NULL;gCC->list=NULL;gCC->rc=-1;GCAdC(gb,&gGC);
A.ty=NIL;A.d.m=0;gpl=NULL;gDI=NULL;arr=0;ar=NULL;gv=NULL;gC=NULL;ix.ty=NIL;ix.d.m=0;gl=NULL;ai=0;BSL();gpl=NULL;ret=0;br=0;ct=0;hl=0;gsd;skd;CCS(&gStackTrace);le=0;pfp=NULL;e=IN(tree);w (Pop());
q e;}o RdPAST(S** g){const uint empty_production=0xFF;uint c0p[]={7,8,9,10,11,12,33,38,40,47,52,55,58,60,64,65};uint c1p[] ={19,24,30,59};uint c2p[]={23};uint mt[]={3,5,6,17,18,78};uint i;
uint rd=1;w (rd){rd=0;y(i=0;i < z(mt)/z(u);i++){f((*g)->pr==mt[i]){(*g)->pr=empty_production;rd=1;}}y(i=0;i < z(c0p)/z(u);i++){f((*g)->pr==c0p[i]){
S* cr=*g;uint ch;y(ch=1;ch < cr->numC;ch++){FPT(cr->c[ch]);}*g=(*g)->c[0];D(cr->c);D(cr);rd=1;}}y(i=0;i < z(c1p)/z(u);i++){f((*g)->pr==c1p[i]){
S* cr=*g;uint ch;y(ch=0;ch < cr->numC;ch++){f(ch != 1){FPT(cr->c[ch]);}}*g=(*g)->c[1];D(cr->c);D(cr);rd=1;}}y(i=0;i < z(c2p)/z(u);i++){f((*g)->pr==c2p[i]){S* cr=*g;uint ch;y(ch=0;ch < cr->numC;ch++){
f(ch != 2){FPT(cr->c[ch]);}}*g=(*g)->c[2];D(cr->c);D(cr);rd=1;}}}y(i=0;i < (*g)->numC;i++){RdPAST(&(*g)->c[i]);}}
o lab(L** str,L** h,L* x){f(x){f(x->sr){char* B=(char*)M(z(char)*(x->length+1));u i;y(i=0;i < x->length;i++)B[i]=x->sr[i];B[i]='\0';x->sr=B;}f(*str)(*h)->x=x;ls*str=x;*h=x;}}
L* eof(u line){L* bl=(L*)M(z(L));bl->tk=gSymbolEOF;bl->sr=0;bl->length=0;bl->line=line;bl->x=0;q bl;}L* newl(u line)
{L* bl=(L*)M(z(L));bl->tk=gEnd;bl->sr=0;bl->length=0;bl->line=line;bl->x=0;q bl;}L* IT(char* src,u str,u h,u line){
L* bl=(L*)M(z(L));bl->tk=gInt;bl->sr=src+str;bl->length=h-str;bl->line=line;bl->x=0;q bl;}L* FLOAT(char* src,
u str,u h,u line){L* bl=(L*)M(z(L));bl->tk=gSF;bl->sr=src+str;bl->length=h-str;bl->line=line;bl->x=0;q bl;}
L* ST(char* src,u str,u h,u line){L* bl=(L*)M(z(L));bl->tk=gSS;bl->sr=src+str+1;bl->length=h-str-2;bl->line=line;bl->x=0;q bl;}
L* TO(char* src,u str,u h,u line,GT g){u tk;L* bl=(L*)M(z(L));tk=FT(&src[str],h-str,g);f(tk){
bl->tk=tk;bl->sr=src+str;bl->length=h-str;bl->line=line;bl->x=0;}ls{D(bl);bl=NULL;p("Syntax error on line %i,",line+1);p("Error, illegal identifier:");w (str < h)
p("0x%X",src[str++]);p("\n");}q bl;}L* ID(char* src,u str,u h,u line,GT g){L* bl;u tk;bl =(L*)M(z(L));tk=FT(&src[str],h-str,g);f(tk){
bl->tk=tk;bl->sr=src+str;bl->length=h-str;bl->line=line;bl->x=0;}ls{bl->tk=gSI;bl->sr=src+str;bl->length=h-str;bl->line=line;bl->x=0;}
q bl;}u isAlpha(char c){q (('a' <= c&&c <= 'z')||('A' <= c&&c <= 'Z'));}u iN(char c){q ('0' <= c&&c <= '9');}
u iS(char c){q (c==' '||c=='\t'||c=='\r'||c=='\n');}u iG(char c){q (!isAlpha(c)&&!iN(c)&&!iS(c));}L* LexSB(cc* sb,char** sp,
GT tb){cc* B;char* format;u* lineN;u line;u sz;u read;u i,j;u a,b;L* str;L* h;B=sb;sz=strlen(B);
lineN=(u*)M((sz+1)*z(u));line=0;y(i=0;i <= sz;i++){lineN[i]=line;f(B[i]=='\n')line++;}format=(char*)M(sz+1);j=0;y(i=0;i < sz;i++){f(B[i]=='/'&&B[i+1]=='/'){
w (i < sz&&B[i] != '\n') i++;i--;}ls f(B[i]=='/'&&B[i+1]=='*'){i++;i++;w (i < sz-1&&(B[i] != '*'||
B[i+1] != '/')){i++;}i++;}ls f(B[i]=='#'){w (i < sz&&B[i] != '\n') i++;i--;}ls f(B[i]==';'){
w (i < sz&&B[i] != '\n') i++;i--;}ls f(iS(B[i])){format[j]=' ';lineN[j]=lineN[i];w (i < sz&&iS(B[i]))
{f(B[i]=='\n'){format[j]='\n';lineN[j]=lineN[i];}i++;}i--;j++;}ls f(B[i]=='"'||B[i]=='\''){char quote=B[i];format[j++]=B[i++];w (i < sz-1&&B[i] != quote)
{f(B[i]=='\\'){lineN[j]=lineN[i];format[j++]=B[i++];}lineN[j]=lineN[i];format[j++]=B[i++];}lineN[j]=lineN[i];format[j++]=B[i];}ls{lineN[j]=lineN[i];format[j]=B[i];j++;}
}format[j]=0;sz=j;str=NULL;h=NULL;y(i=0;i < sz;i++){line=lineN[i];f(isAlpha(format[i])){a=i;
w (i < sz&&(isAlpha(format[i])||iN(format[i])||format[i]=='_')){i++;}b=i;lab(&str,&h,ID(format,a,b,line,tb));i--;}
ls f(iN(format[i])){a=i;w (i < sz&&iN(format[i])){i++;}f(format[i]=='.'){i++;w (i < sz&&iN(format[i])){i++;}
b=i;lab(&str,&h,FLOAT(format,a,b,line));i--;}ls{b=i;lab(&str,&h,IT(format,a,b,line));i--;}}ls f(format[i]=='"'||
format[i]=='\''){char delimiter=format[i];a=i++;w (i < sz&&format[i] != delimiter){i++;}b=i++;lab(&str,&h,ST(format,a,b+1,line));i--;}
ls f(iG(format[i])){a=i;b=0;w (i < sz&&iG(format[i])){i++;f(FT(&format[a],i-a,tb))b=i;}b=b?b:i;
lab(&str,&h,TO(format,a,b,line,tb));i=b-1;}ls f(format[i]=='\n'){lab(&str,&h,newl(line));}
ls f(format[i]==' '){}ls{p("Syntax error on line %i:",line+1);p("Illegal token %i.\n",format[i]);FL(str,format);*sp=0;q NULL;}
}lab(&str,&h,newl(line));lab(&str,&h,eof(line));D(lineN);*sp=format;q str;}o FL(L* ex,
char* B){L* bl;w (ex){bl=ex;ex=ex->x;f(bl->sr){D(bl->sr);}D(bl);}f(B) D(B);}MP* gBM;MP* gWM;o FE(){ClearAllGC();FGCMO(gGC);}V get(cc* id,C* cx)
{V nill;w (cx){P* itr=cx->list;w (itr){f(strcmp(id,itr->id)==0){q itr->v;}itr=itr->x;}cx=cx->par;}
nill.ty=NIL;nill.d.m=0;q nill;}o store(cc* id,V v,C* cx){C* top=cx;w (cx){P* itr=cx->list;w (itr){f(strcmp(id,itr->id)==0)
{itr->v=v;q;}itr=itr->x;}cx=cx->par;}f(top){f(top->list){P* itr=top->list;w (itr->x){itr=itr->x;}itr->x=(P*)M(z(P));itr->x->v=v;itr->x->id=id;itr->x->x=NULL;}ls{
top->list=(P*)M(z(P));top->list->v=v;top->list->id=id;top->list->x=NULL;}}}o RR(cc* id,C* cx){P* prev=NULL;w (cx){P* itr=cx->list;w (itr)
{f(strcmp(id,itr->id)==0){f(prev){prev->x=itr->x;}ls{cx->list=itr->x;}D (itr);q;}prev=itr;itr=itr->x;}cx=cx->par;}}uint HashFn(V v){uint hash;cc* x;switch (v.ty){
case NIL:q 0;case VP:q (uint)v.d.m;case FP:q (uint)v.d.fp;case VR:q (unsigned long)v.d.ref;case VF:q (unsigned long)v.d.fn;case DT:q (unsigned long)v.d.dt;case VS:
hash=0;x=v.d.sr;w (*x){hash=(hash << 7) ^ hash;hash += (*x);x++;}q hash;}q 0;}H* cht()
{H* ht=(H*)M(z(H));ht->cap=HTMC;ht->sz=0;ht->tb=(KVP*)M(HTMC*z(KVP));memset((o*)ht->tb,0,HTMC*z(KVP));q ht;}H* chtN(u htsize)
{H* ht;f(htsize <= 2) htsize=HTMC;ht=(H*)M(z(H));ht->cap=htsize;ht->sz=0;ht->tb=(KVP*)M(htsize*z(KVP));memset((o*)ht->tb,0,htsize*z(KVP));q ht;}o fh(H* tb){D(tb->tb);D(tb);}o ResizeHashTable(H* tb)
{u itr;u cap;uint n;KVP pair;u new_cap;KVP* new_tab;new_cap=tb->cap*HTRF;new_tab=(KVP*)M(new_cap*z(KVP));memset((o*)new_tab,0,new_cap*z(KVP));cap=tb->cap;
y(itr=0;itr < cap;itr++){f(tb->tb[itr].key.ty != NIL){pair=tb->tb[itr];n=HashFn(pair.key) % new_cap;w (new_tab[n].key.ty != NIL)
{n=(n+1) % new_cap;}new_tab[n]=pair;}}tb->cap=new_cap;D(tb->tb);tb->tb=new_tab;}V HG(V kd,H* tb)
{uint n;KVP* te;n=HashFn(kd);n=n % tb->cap;te=tb->tb;f(kd.ty==VS){w (te[n].key.ty != NIL&&!(te[n].key.ty==VS&&
strcmp(te[n].key.d.sr,kd.d.sr)==0)){n=(n+1) % tb->cap;}}ls{w (te[n].key.ty != NIL&&
!(kd.ty==te[n].key.ty&&kd.d.m==te[n].key.d.m)){n=(n+1) % tb->cap;}}f(te[n].key.ty==NIL){
V null;null.ty=NIL;null.d.m=0;q null;}q te[n].v;}o hs(V kd,V store,H* tb){uint n;KVP* te;f(kd.ty==NIL) q;f(tb->sz > 5*tb->cap/7)
{ResizeHashTable(tb);}n=HashFn(kd);n=n % tb->cap;te=tb->tb;f(kd.ty==VS){w (te[n].key.ty != NIL&&!(te[n].key.ty==VS&&
strcmp(te[n].key.d.sr,kd.d.sr)==0)){n=(n+1) % tb->cap;}}ls{w (te[n].key.ty != NIL&&
!(kd.ty==te[n].key.ty&&kd.d.m==te[n].key.d.m)){n=(n+1) % tb->cap;}}f(te[n].key.ty==NIL) tb->sz++;te[n].key=kd;te[n].v=store;q;
}N At(LR tb,u stt,u tk){N ac;ac.ty=AE;f(tk & KT)tk=tk ^ KT;f(stt >= 0&&stt < tb.numStates&&tk > 0&&tk <= tb.nt){ac=tb.actionTable[stt*tb.nt+tk-1];}q ac;}
u Go(LR tb,u stt,u sym){u nS=-1;f(sym & KS)sym=sym ^ KS;f(stt >= 0&&stt < tb.numStates&&sym > 0&&sym <= tb.nS){nS=tb.gotoTable[stt*tb.nS+sym-1];}q nS;}T* gps=NULL;T* gpt=NULL;
o StackPushState(u stt){f(gpt){gpt->x=(T*)M(z(T));gpt->x->prev=gpt;gpt=gpt->x;}ls{gps=gpt=(T*)M(z(T));gpt->prev=NULL;}gpt->x=NULL;gpt->tk=NULL;gpt->ty=0;gpt->stt=stt;}
o StackPushToken(L* tk){f(gpt){gpt->x=(T*)M(z(T));gpt->x->prev=gpt;gpt=gpt->x;}ls{gps=gpt=(T*)M(z(T));gpt->prev=NULL;}gpt->x=NULL;gpt->tk=tk;gpt->ty=TL;gpt->stt=-1;}
o StackPushTree(S* tree){f(gpt){gpt->x=(T*)M(z(T));gpt->x->prev=gpt;gpt=gpt->x;}ls{gps=gpt=(T*)M(z(T));gpt->prev=NULL;}gpt->x=NULL;gpt->tk=tree;gpt->ty=TST;gpt->stt=-1;}
T StackPeek(){T r;r.tk=NULL;r.ty=0;r.stt=-1;f(gpt){r=*gpt;r.x=NULL;r.prev=NULL;}q r;}T StackPop(){T r;r.tk=NULL;r.ty=0;r.stt=-1;
f(gpt){r=*gpt;r.x=NULL;r.prev=NULL;f(gpt->prev){gpt=gpt->prev;D(gpt->x);gpt->x=NULL;}ls{D(gpt);gpt=NULL;gps=NULL;}}q r;}o PE(L* in,GT g){u i;f(in){p("Syntax error on line %i:",(in->line+1));f(in->sr)
{y(i=0;i < in->length;i++){p("%c",in->sr[i]);}}ls{p("%s",GetE(in->tk,g));}p("\n");}}u PS(L* in,LR pr,GT g){L* ip=in;S* ast=NULL;u successful=0;
StackPushState(0);y(;;){N ac;T s=StackPeek();f(s.tk||s.stt==-1||s.stt >= pr.numStates){successful=-1;K;}f(ip==NULL){successful=0;K;}ac=At(pr,s.stt,ip->tk);f(ac.ty==SH)
{StackPushToken(ip);StackPushState(ac.v);ip=ip->x;}ls f(ac.ty==RD){S* no;R r;u rhs;f(ac.v < 0||ac.v >= g.nR){
successful=-1;K;}r=g.rules[ac.v];y(rhs=0;rhs < r.rhsL;rhs++){StackPop();StackPop();}s=StackPeek();f(s.tk||s.stt==-1||s.stt >= pr.numStates){successful=-1;K;}StackPushTree(NULL);
StackPushState(Go(pr,s.stt,r.lhs));}ls f(ac.ty==AA){successful=1;K;}ls{f(ip->tk==gSymbolEOF)successful=0;ls
successful=-1;K;}}w (gps)StackPop();q successful;}S* PSc(L* in,LR pr,GT g){S* ast=NULL;L* ip=in;StackPushState(0);y(;;){N ac;T s=StackPeek();f(s.tk||s.stt==-1
||s.stt >= pr.numStates||ip==NULL){f(ip)p(es);ls p(es);q 0;}ac=At(pr,s.stt,ip->tk);f(ac.ty==SH)
{StackPushToken(ip);StackPushState(ac.v);ip=ip->x;}ls f(ac.ty==RD){S* no;R r;u rhs;f(ac.v < 0||ac.v >= g.nR){
p(es);PE(ip,g);q 0;}r=g.rules[ac.v];no=(S*)M(z(S));no->tk=r.lhs;no->pr=ac.v;no->sr=NULL;no->length=0;no->line=0;
no->c=(S**)M(r.rhsL*z(S*));no->numC=r.rhsL;y(rhs=0;rhs < r.rhsL;rhs++){T sym;u ch;StackPop();sym=StackPop();f(sym.tk==NULL) {p(es);PE(ip,g);q 0;} ch=r.rhsL-rhs-1;f(sym.ty==TL)
{L* tk=(L*)sym.tk;no->c[ch]=(S*)M(z(S));no->c[ch]->tk=tk->tk;no->c[ch]->pr=0;no->c[ch]->sr=tk->sr;no->c[ch]->length=tk->length;no->c[ch]->line=tk->line;no->c[ch]->c=NULL;no->c[ch]->numC=0;}
ls f(sym.ty==TST){no->c[ch]=(S*)sym.tk;}ls{p(es);PE(ip,g);}}f(r.rhsL) no->line=no->c[0]->line;s=StackPeek();f(s.tk||s.stt==-1||s.stt >= pr.numStates){p(es);PE(ip,g);q 0;}
f(r.lhs==gGoal){ast=no;}StackPushTree(no);StackPushState(Go(pr,s.stt,r.lhs));}ls f(ac.ty==AA){K;}ls{p(es);PE(ip,g);q 0;}}w (gps)StackPop();q ast;}o FPT(S* ast){f(ast&&ast->c){u i;y(i=0;i < ast->numC;i++){
FPT(ast->c[i]);}D(ast->c);D(ast);}ls f(ast){D(ast);}}o PrintParseTreeNode(S* ast,GT g){f(!ast) q;f(ast->sr&&ast->tk != gEnd){u i;y(i=0;i < ast->length;i++)p("%c",ast->sr[i]);}ls
{p("%s",GetE(ast->tk,g));}}o PrintParseTreeRecurse(S* ast,u indent,GT g){f(ast->c){u i;y(i=0;i < indent*2;i++){p(" ");}PrintParseTreeNode(ast,g);p("::= ");
y(i=0;i < ast->numC;i++){PrintParseTreeNode(ast->c[i],g);p(" ");}p("\n");y(i=0;i < ast->numC;i++){PrintParseTreeRecurse(ast->c[i],indent+1,g);}}}
o PrintParseTree(S* ast,GT g){p("\n\nAbstract Syntax Tree\n");p("====================\n\n");PrintParseTreeRecurse(ast,0,g);}
o PrintParseTreeFormat(S* ast,GT g){f(ast->numC > 0){u i;y(i=0;i < ast->numC;i++)
PrintParseTreeFormat(ast->c[i],g);}ls{f(ast->sr&&ast->tk != gEnd){u i;y(i=0;i < ast->length;i++)p("%c",ast->sr[i]);}ls{p("%s",GetE(ast->tk,g));}}}
char GRAMMAR_SYMBOLS[283] ="<g>\0<stmt list>\0<stmt>\0<ref>\0<xp>\0<assignment>\0<function def>\0<if>\0<for loop>\0<while loop>\0<parameters>\0"
"<param decl>\0<cd>\0<else if>\0<arithmetic>\0<l-value>\0<logic>\0<comparison>\0<term>\0<factor>\0"
"<final>\0<boolean>\0<object>\0<as>\0<array init>\0<dictionary init>\0";
char GRAMMAR_TOKENS[218] =
"<epsilon>\0<$>\0<endl>\0<integer>\0<float>\0<string>\0<id>\0import\0"
"call\0return\0break\0continue\0function\0end\0(\0)\0,\0if\0then\0else\0for\0"
"=\0to\0do\0loop\0while\0.\0[\0]\0and\0or\0nor\0xor\0not\0==\0!=\0<\0>\0<=\0>=\0"
"+\0-\0&\0|\0*\0/\0!\0:\0true\0false\0";
u R0[]={4097,0};u R1[]={4098,0};u R2[]={4099,4098,0};u R3[]={8193,0};u R4[]={8200,8199,8195,0};u R5[]={8201,4100,8195,0};u R6[]={8195,0};
u R7[]={4101,8195,0};u R8[]={4102,8195,0};u R9[]={4103,8195,0};u R10[]={4104,8195,0};u R11[]={4105,8195,0};u R12[]={4106,8195,0};
u R13[]={8202,4101,8195,0};u R14[]={8203,8195,0};u R15[]={8204,8195,0};u R16[]={8205,8199,4107,8195,4098,8206,0};u R17[]={8193,0};u R18[]={8207,8208,0};
u R19[]={8207,4108,8208,0};u R20[]={8199,0};u R21[]={4108,8209,8199,0};u R22[]={8210,4109,8211,8195,4098,4110,0};u R23[]={8212,8195,4098,8206,0};
u R24[]={8212,4104,0};u R25[]={8213,8199,8214,4111,8215,4111,8216,8195,4098,8217,0};u R26[]={8218,4109,8216,8195,4098,8217,0};u R27[]={4112,8214,4102,0};
u R28[]={4112,8214,4109,0};u R29[]={8199,0};u R30[]={8207,4112,8208,0};u R31[]={4100,8219,8199,0};u R32[]={4100,8220,4101,8221,0};u R33[]={4109,0};
u R34[]={4109,8222,4113,0};u R35[]={4109,8223,4113,0};u R36[]={4109,8224,4113,0};u R37[]={4109,8225,4113,0};u R38[]={4113,0};u R39[]={8226,4114,0};
u R40[]={4114,0};u R41[]={4114,8227,4111,0};u R42[]={4114,8228,4111,0};u R43[]={4114,8229,4111,0};u R44[]={4114,8230,4111,0};u R45[]={4114,8231,4111,0};
u R46[]={4114,8232,4111,0};u R47[]={4111,0};u R48[]={4111,8233,4115,0};u R49[]={4111,8234,4115,0};u R50[]={4111,8235,4115,0};u R51[]={4111,8236,4115,0};
u R52[]={4115,0};u R53[]={4115,8237,4116,0};u R54[]={4115,8238,4116,0};u R55[]={4116,0};u R56[]={8234,4116,0};u R57[]={8239,4116,0};
u R58[]={4117,0};u R59[]={8207,4101,8208,0};u R60[]={4118,0};u R61[]={8196,0};u R62[]={8197,0};u R63[]={8198,0};u R64[]={4119,0};
u R65[]={4100,0};u R66[]={4112,0};u R67[]={4100,8207,8208,0};u R68[]={4100,8207,4120,8208,0};u R69[]={4120,8209,4101,0};u R70[]={4101,0};
u R71[]={8220,8221,0};u R72[]={8220,4121,8221,0};u R73[]={8220,4122,8221,0};u R74[]={4121,8209,4101,0};u R75[]={4101,0};u R76[]={4122,8209,8199,8240,4101,0};
u R77[]={8199,8240,4101,0};u R78[]={8206,0};u R79[]={8241,0};u R80[]={8242,0}; R GRAMMAR_RS[81] ={
0,1,R0,4097,1,R1,4098,2,R2,4098,0,R3,4099,3,R4,4099,3,R5,4099,1,R6,4099,2,R7,4099,2,R8,
4099,2,R9,4099,2,R10,4099,2,R11,4099,2,R12,4099,3,R13,4099,2,R14,4099,2,R15,4103,6,R16,4107,0,R17,
4107,2,R18,4107,3,R19,4108,1,R20,4108,3,R21,4104,6,R22,4110,4,R23,4110,2,R24,4105,10,R25,4106,6,R26,
4102,3,R27,4102,3,R28,4112,1,R29,4112,3,R30,4112,3,R31,4112,4,R32,4101,1,R33,4109,3,R34,4109,3,R35,
4109,3,R36,4109,3,R37,4109,1,R38,4113,2,R39,4113,1,R40,4114,3,R41,4114,3,R42,4114,3,R43,4114,3,R44,
4114,3,R45,4114,3,R46,4114,1,R47,4111,3,R48,4111,3,R49,4111,3,R50,4111,3,R51,4111,1,R52,4115,3,R53,
4115,3,R54,4115,1,R55,4116,2,R56,4116,2,R57,4116,1,R58,4117,3,R59,4117,1,R60,4117,1,R61,4117,1,R62,
4117,1,R63,4117,1,R64,4117,1,R65,4100,1,R66,4100,3,R67,4100,4,R68,4120,3,R69,4120,1,R70,4119,2,R71,
4119,3,R72,4119,3,R73,4121,3,R74,4121,1,R75,4122,5,R76,4122,3,R77,4110,1,R78,4118,1,R79,4118,1,R80,
};GT CFG ={26,GRAMMAR_SYMBOLS,50,GRAMMAR_TOKENS,81,GRAMMAR_RS}; u GOTO_TABLE[]= {
1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,-1,2,11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,56,43,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,-1,2,11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,
20,1,-1,630,71,1,-1,11,72,1,-1,13,4,1,74,1,-1,7,11,1,-1,1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,84,4,1,78,1,-1,7,11,1,-1,1,12,1,79,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,8,80,1,-1,
1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,32,4,1,-1,8,82,1,-1,1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,83,1,-1,7,11,1,-1,1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,1,84,1,
85,1,-1,3,4,1,-1,10,12,1,72,1,-1,1,88,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,11,72,1,-1,3,89,1,18,1,19,1,20,1,-1,6,4,1,-1,11,72,1,-1,3,90,1,18,1,19,1,20,1,-1,84,4,1,91,1,-1,7,11,1,-1,1,12,1,72,1,14,1,
15,1,16,1,17,1,18,1,19,1,20,1,92,1,-1,31,4,1,95,1,-1,7,11,1,-1,1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,162,4,1,-1,10,12,1,72,1,96,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,12,1,72,1,97,1,
15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,12,1,72,1,98,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,12,1,72,1,99,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,11,72,1,-1,2,100,1,17,1,18,1,19,1,20,1,
-1,6,4,1,-1,11,72,1,-1,2,101,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,11,72,1,-1,2,102,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,11,72,1,-1,2,103,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,1,104,1,-1,6,105,1,-1,1,12,1,13,1,14,
1,15,1,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,106,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,107,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,108,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,
6,4,1,-1,10,109,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,110,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,10,111,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,6,4,1,-1,11,72,1,-1,3,112,1,18,1,19,
1,20,1,-1,6,4,1,-1,11,72,1,-1,3,113,1,18,1,19,1,20,1,-1,84,116,1,-1,11,79,1,-1,98,118,1,-1,1118,135,1,-1,95,4,1,-1,10,139,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,32,4,1,141,1,-1,7,11,1,-1,1,12,1,72,1,14,
1,15,1,16,1,17,1,18,1,19,1,20,1,-1,84,4,1,143,1,-1,7,11,1,-1,1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,32,4,1,144,1,-1,7,11,1,-1,1,12,1,72,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,30,145,1,3,1,4,1,
5,1,6,1,7,1,8,1,9,1,10,1,-1,2,11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1, 82,148,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,-1,2,11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,30,150,1,
3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,-1,2,11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,198,154,1,-1,15,4,1,-1,10,157,1,72,1,-1,2,16,1,17,1,18,1,19,1,20,1,-1,32,4,1,159,1,-1,7,11,1,-1,1,12,1,72,1,
14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,114,160,1,-1,123,163,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,-1,2,11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,56,166,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,-1,2,
11,1,-1,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19,1,20,1,-1,81}; u ACT_A[]={
0,1,2,1,1,11,2,1,1,1,0,2,1,1,0,1,2,1,1,1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,1,3,1,0,49,2,1,0,49,2,1,1,11,2,1,1,1,0,2,1,1,0,1,2,1,1,1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,
2,2,1,0,11,1,1,2,2,0,1,2,1,0,3,2,2,0,2,1,2,2,5,0,1,2,12,0,6,1,1,0,49,1,1,0,49,1,1,0,49,1,1,0,49,1,1,0,49,1,1,0,49,2,1,0,12,2,2,0,11,2,1,1,4,0,19,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,2,6,1,4,0,8,2,1,0,11,
2,3,0,1,2,1,0,2,1,1,2,2,0,2,2,7,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,19,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,1,6,0,12,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,10,1,2,0,6,2,1,0,12,2,2,0,
1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,55,2,14,0,2,2,1,0,1,2,2,
0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,2,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,
0,11,2,3,0,1,2,1,0,2,2,3,0,2,2,7,0,1,2,12,0,10,1,1,0,49,1,1,0,7,1,1,0,38,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,1,1,0,49,1,1,0,53,1,1,0,46,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,
1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,6,1,1,0,46,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,2,0,4,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,
1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,2,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,
2,2,0,4,2,5,0,1,2,12,0,5,2,1,0,11,2,1,0,5,2,1,0,4,2,1,0,28,1,4,0,7,1,2,0,11,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,6,1,1,0,46,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,1,2,14,0,2,2,1,0,1,2,2,0,
3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,
1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,
3,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,
1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,
1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,
1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,
1,2,0,3,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,2,1,1,0,49,1,1,0,11,1,1,0,11,1,2,0,24,2,1,0,11,2,3,0,1,2,1,0,3,2,2,0,2,2,7,0,1,2,12,0,10,1,1,0,7,1,1,0,37,1,1,0,48,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,
1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,2,2,1,0,11,1,1,0,50,1,1,0,36,2,1,0,11,2,1,1,1,2,1,0,1,2,1,0,3,2,2,0,2,2,7,0,1,2,12,0,22,1,1,0,10,
1,4,0,38,1,1,0,51,1,1,0,5,1,4,0,33,2,1,0,11,2,1,0,37,1,1,0,11,1,1,0,37,1,1,0,11,1,1,0,23,2,1,0,11,2,3,0,1,2,1,0,2,2,3,0,2,2,7,0,1,2,12,0,1,1,1,0,4,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,
2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,1,6,0,12,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,19,2,2,0,48,1,2,0,35,2,1,0,11,2,3,0,1,2,1,0,3,2,2,0,2,2,7,0,1,2,12,
0,6,2,1,0,11,2,3,0,1,2,1,0,2,2,3,0,2,2,7,0,1,2,12,0,32,1,1,0,23,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,19,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,19,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,19,2,1,0,12,2,2,0,1,
2,1,0,4,2,1,0,4,2,5,0,19,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,10,1,2,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,10,1,2,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,10,1,2,0,6,2,1,0,12,2,2,0,1,2,
1,0,3,2,2,0,4,2,5,0,1,2,10,1,2,0,6,2,1,0,49,2,1,0,26,1,4,0,19,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,2,6,1,4,0,8,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,2,6,1,4,0,8,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,
1,2,6,1,4,0,8,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,2,6,1,4,0,8,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,2,6,1,4,0,8,2,1,0,12,2,2,0,1,2,1,0,4,2,1,0,4,2,5,0,1,2,6,1,4,0,8,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,
2,5,0,1,2,12,0,6,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,5,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,1,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,
2,0,14,1,1,0,11,1,2,0,23,2,14,0,2,2,1,0,1,2,2,0,3,2,2,0,1,2,1,0,5,2,1,0,7,2,1,0,4,2,1,0,1,2,2,0,2,1,1,0,103,1,1,0,8,1,1,0,36,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,6,2,1,0,11,2,3,0,1,2,1,0,2,2,3,0,2,
2,7,0,1,2,12,0,6,1,1,0,50,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,2,1,1,0,50,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,2,1,0,12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,10,1,1,0,45,2,1,0,
12,2,2,0,1,2,1,0,3,2,2,0,4,2,5,0,1,2,12,0,7,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,2,1,0,11,2,3,0,1,2,1,0,3,2,2,0,2,2,7,0,1,2,12,0,7,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,
2,1,0,11,2,3,0,1,2,1,0,2,2,3,0,2,2,7,0,1,2,12,0,5,2,1,1,11,2,1,1,1,0,2,1,1,0,1,2,1,1,1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,15,1,2,0,48,2,2,0,35,2,1,0,48,2,1,1,11,2,1,1,1,0,2,1,1,0,1,2,1,1,
1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,22,1,1,0,17,1,4,0,7,2,1,1,11,2,1,1,1,0,2,1,1,0,1,2,1,1,1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,16,2,1,0,11,2,1,0,68,1,1,0,18,2,1,0,11,
2,1,0,36,2,2,0,46,1,1,0,38,2,1,0,53,1,1,0,56,1,1,0,5,1,1,0,33,1,4,0,7,1,1,0,12,1,1,0,13,1,1,0,4,1,1,0,1,1,2,0,24,1,1,0,28,1,4,0,7,1,1,0,12,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,2,1,0,62,2,2,0,35,2,1,0,49,
2,1,0,49,1,1,0,14,1,1,0,55,1,1,0,16,1,4,0,8,2,1,0,63,2,1,0,11,2,1,0,23,2,1,0,48,2,1,1,11,2,1,1,1,0,2,1,1,0,1,2,1,1,1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,1,1,0,60,1,1,0,37,2,1,1,11,2,1,1,
1,0,2,1,1,0,1,2,1,1,1,0,3,2,1,1,1,0,1,1,1,0,5,1,1,0,7,1,1,0,4,1,1,0,1,1,2,0,2,2,1,0,71,1,1,0,27,2,1,0,47
};
u ACT_B[]={
0,1,3,1,22,1,23,1,24,1,25,1,26,1,27,1,28,1,29,1,30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,51,1,1,0,49,3,1,22,1,23,1,24,1,25,1,26,1,27,1,28,
1,29,1,30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,65,1,0,11,44,1,65,2,0,1,65,1,0,3,65,2,0,2,45,1,46,1,65,5,0,1,65,12,0,6,47,1,0,49,48,1,0,
49,49,1,0,49,50,1,0,49,51,1,0,49,52,1,0,49,33,1,0,12,33,2,0,11,33,1,53,1,54,1,55,1,56,1,0,19,47,1,0,12,47,2,0,1,47,1,0,4,47,1,0,4,47,5,0,1,47,6,57,1,58,1,59,1,60,1,0,8,66,1,0,11,66,3,0,1,66,1,0,2,61,1,66,2,0,
2,66,7,0,1,66,12,0,6,38,1,0,12,38,2,0,1,38,1,0,4,38,1,0,4,38,5,0,19,40,1,0,12,40,2,0,1,40,1,0,4,40,1,0,4,40,5,0,1,62,1,63,1,64,1,65,1,66,1,67,1,0,12,52,1,0,12,52,2,0,1,52,1,0,3,52,2,0,4,52,5,0,1,52,10,68,1,69,
1,0,6,55,1,0,12,55,2,0,1,55,1,0,3,55,2,0,4,55,5,0,1,55,12,0,6,58,1,0,12,58,2,0,1,58,1,0,3,58,2,0,4,58,5,0,1,58,12,0,6,60,1,0,12,60,2,0,1,60,1,0,3,60,2,0,4,60,5,0,1,60,12,0,6,64,1,0,12,64,2,0,1,64,1,0,3,64,2,0,
4,64,5,0,1,64,12,0,55,6,14,0,2,6,1,0,1,6,2,0,3,6,2,0,1,6,1,0,5,6,1,0,7,6,1,0,4,6,1,0,1,6,2,0,2,61,1,0,12,61,2,0,1,61,1,0,3,61,2,0,4,61,5,0,1,61,12,0,6,62,1,0,12,62,2,0,1,62,1,0,3,62,2,0,4,62,5,0,1,62,12,0,6,63,
1,0,12,63,2,0,1,63,1,0,3,63,2,0,4,63,5,0,1,63,12,0,6,29,1,0,11,29,3,0,1,29,1,0,2,29,3,0,2,29,7,0,1,29,12,0,10,70,1,0,49,26,1,0,7,73,1,0,38,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,
41,1,42,1,0,2,75,1,0,49,76,1,0,53,77,1,0,46,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,
6,81,1,0,46,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,86,1,0,7,33,1,0,12,37,1,87,1,0,4,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,
1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,2,79,1,0,12,
79,2,0,1,79,1,0,3,79,2,0,4,79,5,0,1,79,12,0,6,80,1,0,12,80,2,0,1,80,1,0,3,80,2,0,4,80,5,0,1,80,12,0,5,2,1,0,11,2,1,0,5,2,1,0,4,2,1,0,28,23,1,24,1,25,1,26,1,0,7,33,1,93,1,0,11,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,
41,1,42,1,0,6,94,1,0,46,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,1,7,14,0,2,7,1,0,1,7,2,0,3,7,2,0,1,7,1,0,5,7,1,0,7,7,1,0,4,7,1,0,1,7,2,0,1,8,14,0,2,8,1,0,1,8,2,0,3,8,2,
0,1,8,1,0,5,8,1,0,7,8,1,0,4,8,1,0,1,8,2,0,1,9,14,0,2,9,1,0,1,9,2,0,3,9,2,0,1,9,1,0,5,9,1,0,7,9,1,0,4,9,1,0,1,9,2,0,1,10,14,0,2,10,1,0,1,10,2,0,3,10,2,0,1,10,1,0,5,10,1,0,7,10,1,0,4,10,1,0,1,10,2,0,1,11,14,0,2,
11,1,0,1,11,2,0,3,11,2,0,1,11,1,0,5,11,1,0,7,11,1,0,4,11,1,0,1,11,2,0,1,12,14,0,2,12,1,0,1,12,2,0,3,12,2,0,1,12,1,0,5,12,1,0,7,12,1,0,4,12,1,0,1,12,2,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,
0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,
1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,
1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,
1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,
12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,
26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,3,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,2,
114,1,0,49,115,1,0,11,44,1,0,11,45,1,46,1,0,24,66,1,0,11,66,3,0,1,66,1,0,3,66,2,0,2,66,7,0,1,66,12,0,10,26,1,0,7,73,1,0,37,117,1,0,48,14,14,0,2,14,1,0,1,14,2,0,3,14,2,0,1,14,1,0,5,14,1,0,7,14,1,0,4,14,1,0,1,14,2,
0,1,15,14,0,2,15,1,0,1,15,2,0,3,15,2,0,1,15,1,0,5,15,1,0,7,15,1,0,4,15,1,0,1,15,2,0,2,17,1,0,11,120,1,0,50,121,1,0,36,66,1,0,11,66,1,122,1,66,1,0,1,66,1,0,3,66,2,0,2,66,7,0,1,66,12,0,22,123,1,0,10,53,1,54,1,55,
1,56,1,0,38,124,1,0,51,125,1,0,5,53,1,54,1,55,1,56,1,0,33,75,1,0,11,75,1,0,37,126,1,0,11,127,1,0,37,128,1,0,11,129,1,0,23,29,1,0,11,29,3,0,1,29,1,0,2,29,3,0,2,29,7,0,1,29,12,0,1,130,1,0,4,71,1,0,12,71,2,0,1,71,
1,0,3,71,2,0,4,71,5,0,1,71,12,0,6,39,1,0,12,39,2,0,1,39,1,0,4,39,1,0,4,39,5,0,1,62,1,63,1,64,1,65,1,66,1,67,1,0,12,56,1,0,12,56,2,0,1,56,1,0,3,56,2,0,4,56,5,0,1,56,12,0,6,57,1,0,12,57,2,0,1,57,1,0,3,57,2,0,4,57,
5,0,1,57,12,0,19,70,2,0,48,131,1,132,1,0,35,67,1,0,11,67,3,0,1,67,1,0,3,67,2,0,2,67,7,0,1,67,12,0,6,31,1,0,11,31,3,0,1,31,1,0,2,31,3,0,2,31,7,0,1,31,12,0,32,133,1,0,23,34,1,0,12,34,2,0,1,34,1,0,4,34,1,0,4,34,5,0,
19,35,1,0,12,35,2,0,1,35,1,0,4,35,1,0,4,35,5,0,19,36,1,0,12,36,2,0,1,36,1,0,4,36,1,0,4,36,5,0,19,37,1,0,12,37,2,0,1,37,1,0,4,37,1,0,4,37,5,0,19,48,1,0,12,48,2,0,1,48,1,0,3,48,2,0,4,48,5,0,1,48,10,68,1,69,1,0,6,49,
1,0,12,49,2,0,1,49,1,0,3,49,2,0,4,49,5,0,1,49,10,68,1,69,1,0,6,50,1,0,12,50,2,0,1,50,1,0,3,50,2,0,4,50,5,0,1,50,10,68,1,69,1,0,6,51,1,0,12,51,2,0,1,51,1,0,3,51,2,0,4,51,5,0,1,51,10,68,1,69,1,0,6,27,1,0,49,28,1,0,
26,53,1,54,1,55,1,56,1,0,19,41,1,0,12,41,2,0,1,41,1,0,4,41,1,0,4,41,5,0,1,41,6,57,1,58,1,59,1,60,1,0,8,42,1,0,12,42,2,0,1,42,1,0,4,42,1,0,4,42,5,0,1,42,6,57,1,58,1,59,1,60,1,0,8,43,1,0,12,43,2,0,1,43,1,0,4,43,1,0,
4,43,5,0,1,43,6,57,1,58,1,59,1,60,1,0,8,44,1,0,12,44,2,0,1,44,1,0,4,44,1,0,4,44,5,0,1,44,6,57,1,58,1,59,1,60,1,0,8,45,1,0,12,45,2,0,1,45,1,0,4,45,1,0,4,45,5,0,1,45,6,57,1,58,1,59,1,60,1,0,8,46,1,0,12,46,2,0,1,46,
1,0,4,46,1,0,4,46,5,0,1,46,6,57,1,58,1,59,1,60,1,0,8,53,1,0,12,53,2,0,1,53,1,0,3,53,2,0,4,53,5,0,1,53,12,0,6,54,1,0,12,54,2,0,1,54,1,0,3,54,2,0,4,54,5,0,1,54,12,0,5,4,14,0,2,4,1,0,1,4,2,0,3,4,2,0,1,4,1,0,5,4,1,0,
7,4,1,0,4,4,1,0,1,4,2,0,1,5,14,0,2,5,1,0,1,5,2,0,3,5,2,0,1,5,1,0,5,5,1,0,7,5,1,0,4,5,1,0,1,5,2,0,14,44,1,0,11,45,1,46,1,0,23,13,14,0,2,13,1,0,1,13,2,0,3,13,2,0,1,13,1,0,5,13,1,0,7,13,1,0,4,13,1,0,1,13,2,0,2,134,
1,0,103,136,1,0,8,137,1,0,36,59,1,0,12,59,2,0,1,59,1,0,3,59,2,0,4,59,5,0,1,59,12,0,6,30,1,0,11,30,3,0,1,30,1,0,2,30,3,0,2,30,7,0,1,30,12,0,6,138,1,0,50,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,
1,41,1,42,1,0,2,140,1,0,50,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,72,1,0,12,72,2,0,1,72,1,0,3,72,2,0,4,72,5,0,1,72,12,0,10,142,1,0,45,73,1,0,12,73,2,0,1,73,1,0,3,73,2,
0,4,73,5,0,1,73,12,0,7,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,68,1,0,11,68,3,0,1,68,1,0,3,68,2,0,2,68,7,0,1,68,12,0,7,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,
7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,32,1,0,11,32,3,0,1,32,1,0,2,32,3,0,2,32,7,0,1,32,12,0,5,3,1,22,1,23,1,24,1,25,1,26,1,27,1,28,1,29,1,30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,
7,39,1,0,4,40,1,0,1,41,1,42,1,0,15,146,1,147,1,0,48,20,2,0,35,18,1,0,48,3,1,22,1,23,1,24,1,25,1,26,1,27,1,28,1,29,1,30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,
1,41,1,42,1,0,22,149,1,0,17,57,1,58,1,59,1,60,1,0,7,3,1,22,1,23,1,24,1,25,1,26,1,27,1,28,1,29,1,30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,16,74,
1,0,11,74,1,0,68,151,1,0,18,77,1,0,11,77,1,0,36,69,2,0,46,152,1,0,38,19,1,0,53,153,1,0,56,155,1,0,5,156,1,0,33,23,1,24,1,25,1,26,1,0,7,33,1,0,12,37,1,0,13,39,1,0,4,40,1,0,1,41,1,42,1,0,24,158,1,0,28,23,1,24,1,25,1,
26,1,0,7,33,1,0,12,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,16,1,0,62,21,2,0,35,22,1,0,49,78,1,0,49,161,1,0,14,34,1,0,55,162,1,0,16,57,1,58,1,59,1,60,1,0,8,26,1,0,63,76,1,0,11,76,1,0,23,24,1,0,48,3,1,22,1,
23,1,24,1,25,1,26,1,27,1,28,1,29,1,30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,164,1,0,60,165,1,0,37,3,1,22,1,23,1,24,1,25,1,26,1,27,1,28,1,29,1,
30,1,31,1,32,1,3,1,33,1,0,2,34,1,0,1,3,1,35,1,0,3,3,1,36,1,0,1,37,1,0,5,38,1,0,7,39,1,0,4,40,1,0,1,41,1,42,1,0,2,23,1,0,71,167,1,0,27,25,1,0,47
};LR PT={50,26,168,0,0};db TF(V v){f(v.ty==VP){q (db)v.d.m;}ls f(v.ty==FP){q v.d.fp;}ls f(v.ty==VS){q atof(v.d.sr);}ls{q 0.0;}}u TI(V v){f(v.ty==VP){
q v.d.m;}ls f(v.ty==FP){q (u)v.d.fp;}ls f(v.ty==VS){q atoi(v.d.sr);}ls{q 0;}}u IsD(V ty){f(ty.ty==VP||ty.ty==FP||ty.ty==NIL){q 0;}ls f(ty.ty==VS||ty.ty==VF||ty.ty==VR||ty.ty==DT)
{q 1;}q 0;}V Cop(V sr){V copy;f(sr.ty==VS){copy.ty=VS;copy.d.sr=(char*)M(strlen(sr.d.sr)+1);sprintf((char*)copy.d.sr,"%s",sr.d.sr);copy.const_string=0;}
ls{copy.ty=NIL;copy.d.m=0;}q copy;}H* dpr=NULL;o PrintValue(V v){switch (v.ty)
{case VP:p("%i",v.d.m);K;case FP:p("%.16g",v.d.fp);K;case VS:p("%s",v.d.sr);K;case VR:PrintObject(v.d.ref);K;case VF:PrintFn(v.d.fn);K;case DT:PD(v.d.dt);K;default:
case NIL:p("[NIL]");K;}}o PrintFn(F* fn){f(fn){P* itr;f(fn->fm){p("%s(",fn->fm);}ls{p("fn(");}itr=fn->pm;w (itr){p("%s",itr->id);f(itr->x)p(",");itr=itr->x;}
p(")");}}o PrintObject(C* cx){V key;P* list;key.ty=VR;key.d.ref=cx;f(HG(key,dpr).ty != NIL){p("...");q;}ls{V v;v.ty=VP;v.d.m=1;hs(key,v,dpr);}
p("[");list=cx->list;w (list){PrintValue(list->v);f(list->x)p(",");list=list->x;}p("]");}o PD(H* dt){u i,sz;V key;sz=dt->sz;key.ty=DT;key.d.dt=dt;f(HG(key,dpr).ty != NIL){p("...");q;}ls{
V v;v.ty=VP;v.d.m=1;hs(key,v,dpr);}p("[");y(i=0;i < dt->cap;i++){f(dt->tb[i].key.ty != NIL){sz--;PrintValue(dt->tb[i].key);p(":");PrintValue(dt->tb[i].v);f(sz) p(",");}}p("]");}u Parses(u ac)
{L* ex;char* B;u e=0;V a=get("src",gCC);A.ty=VP;A.d.m=0;f(a.ty==VS){ex=LexSB(a.d.sr,&B,CFG);f(ex==NULL){A.d.m=-1;FL(ex,B);q 0;}A.d.m=PS(ex,PT,CFG);f(A.d.m==0){FL(ex,B);q 0;}FL(ex,B);}q e;}u Eval(u ac)
{L* ex;S* ast;char* B;u e=0;C* ccur;V a=get("src",gCC);u ple=le;S* pvp=pfp;A.ty=NIL;A.d.m=0;ccur=gCC;
gCC=gb;f(a.ty==VS){ex=LexSB(a.d.sr,&B,CFG);f(ex==NULL){p(es);FL(ex,B);q 1;}ast=PSc(ex,PT,CFG);f(ast==NULL){p(es);FL(ex,B);q 1;}RdPAST(&ast);e=IN(ast);f(e){p("%s\n",ErrorMessage(e));PST();FL(ex,B);FPT(ast);CCS(&gStackTrace);le=ple;pfp=pvp;q 1;}
GCAdL(ex,B);GCAdP(ast);}gCC=ccur;q e;}u DP(u ac){u e=0;V a=get("output",gCC);dpr=cht();PrintValue(a);fh(dpr);dpr=NULL;A.ty=NIL;A.d.m=0;q e;}
u DPLn(u ac){u e=0;V a=get("output",gCC);dpr=cht();PrintValue(a);fh(dpr);dpr=NULL;A.ty=NIL;A.d.m=0;p("\n");q e;}u Type(u ac){u e=0;V a=get("object",gCC);A.ty=VS;switch (a.ty)
{case NIL:A.d.sr="NIL";K;case VP:A.d.sr="INT";K;case FP:A.d.sr="FLOAT";K;case VS:A.d.sr="STRING";K;case VR:A.d.sr="REFERENCE";K;case VF:A.d.sr="FUNCTION";K;case DT:A.d.sr="DICTIONARY";K;default:A.d.sr="UNKNOWN";}
q e;}u DInt(u ac){u e=0;V a=get("v",gCC);A.ty=VP;A.d.m=TI(a);q e;}u DFl(u ac){u e=0;V a=get("v",gCC);A.ty=FP;A.d.fp=TF(a);q e;}u Len(u ac)
{u e=0;V a=get("array",gCC);f(a.ty==VR){u count=0;C* ref;P* itr;ref=a.d.ref;itr=ref->list;w (itr){count++;itr=itr->x;}
A.ty=VP;A.d.m=count;}ls f(a.ty==DT){ A.ty=VP;A.d.m=a.d.dt->sz;}ls {f(a.ty != NIL){A.ty=VP;A.d.m=1;}ls{A.ty=VP;A.d.m=0;}}q e;}u DQ(u ac){u e=0;hl=1;
A.ty=VP;A.ty=0;q e;}o BSL(){V dsl;V print;V println;V root;V int_c;V float_c;V length;V ty;V evalFunc;V parsesFunc;V duckQuit;
dsl=LN("duck");print=CrF(DP);AP(print,"output");LinkFn(dsl,"print",print);println=CrF(DPLn);AP(println,"output");LinkFn(dsl,"println",println);LCS(dsl,"newline","\n");root.ty=VR;root.d.ref=gb;int_c=CrF(DInt);AP(int_c,"v");
LinkFn(root,"u",int_c);float_c=CrF(DFl);AP(float_c,"v");LinkFn(root,"float",float_c);length=CrF(Len);AP(length,"array");LinkFn(root,"len",length);ty=CrF(Type);AP(ty,"object");LinkFn(root,"Type",ty);evalFunc=CrF(Eval);AP(evalFunc,"src");
LinkFn(root,"eval",evalFunc);parsesFunc=CrF(Parses);AP(parsesFunc,"src");LinkFn(dsl,"parses",parsesFunc);duckQuit=CrF(DQ);LinkFn(root,"quit",duckQuit);BindStringLibrary();}
u StringSplit(u ac){u e=0;V a=get("sr",gCC);f(a.ty==VS){cc* d=a.d.sr;u length=strlen(d);uint n;H* dt=cht();GCAdD(dt,&gGC);y(n=0;n < length;n++)
{char* char_string;V key;V xp;key.ty=VP;key.d.m=n;xp.ty=VS;xp.const_string=0;char_string=(char*)M(z(char)*2);GCAdS(char_string,&gGC);char_string[0]=d[n];char_string[1]='\0';xp.d.sr=char_string;hs(key,xp,dt);}
A.ty=DT;A.d.dt=dt;}ls{A.ty=NIL;A.d.m=0;}q e;} o BindStringLibrary(){V split;V stringLib=LN("sr");split=CrF(StringSplit);AP(split,"sr");LinkFn(stringLib,"split",split);}
cc* ErrorMessage(u e){q es;}u main(u argc,char* argv[]){cc* g;L* ex;S* ast;char* B;u e;uint i;uint k;uint j;u val;TAG=M(z(u)*4368);ACT=M(z(N)*8400);PT.actionTable=ACT;PT.gotoTable=TAG;
y(i=0,j=0;i<4368;){val=GOTO_TABLE[j];y(k=GOTO_TABLE[j+1];k;k--){TAG[i++]=val;}j+=2;}y(i=0,j=0;i<8400;){val=ACT_A[j];y(k=ACT_A[j+1];k;k--){ACT[i++].ty=val;}j+=2;}
y(i=0,j=0;i<8400;){val=ACT_B[j];y(k=ACT_B[j+1];k;k--){ACT[i++].v=val;}j+=2;}f(argc > 1){g=argv[1];f(strcmp(argv[1],"-h")==0||strcmp(argv[1],"-help")==0||strcmp(argv[1],"--help")==0)
{p("Usage:duck g.src\n");q 1;}}f(argc > 1){B=0;ex=LexSB(/*g*/demo,&B,CFG);f(ex==NULL){p(es);FL(ex,B);getchar();q 1;}}ls{B=0;ex=LexSB(demo,&B,CFG);f(ex==NULL){
p("Error 1.\n");FL(ex,B);getchar();q 1;}}ast=PSc(ex,PT,CFG);f(ast==NULL){p(es);FL(ex,B);getchar();q 1;}RdPAST(&ast);e=Inter(ast);f(e){p("%s\n",ErrorMessage(e));PST();FL(ex,B);FPT(ast);getchar();q 1;}
FE();FL(ex,B);FPT(ast);D(TAG);D(ACT);q 0;}