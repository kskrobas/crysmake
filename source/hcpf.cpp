//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "hcpf.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

#include <System.Character.hpp>
#include <limits.h>
#include <float.h>
#include <algorithm>
#include "genutil.h"


cdouble wi=0.5;//sin(30)
cdouble wj=sqrt(3)/2;//sin(60)
cdouble k2s3=2.0*sqrt(3);

//double StHcpPrm::tolPos=1.0/1024;


wofstream & operator<<(wofstream &o, const StHcpPrm &hcp)
{
		o<<"stack sequence: \r\n"<<(hcp.stackSeq.c_str())<<endl
		 <<"grain Parameters: "<<endl;
		o<<(*hcp.grainPrm)<<endl;

return o;
}


ostream & operator<<(ostream &o, const StHcpPrm &hcp)
{
	if(hcp.stackSeq.Length() && hcp.store){
	o<<"** hcp stack sequence: \r\n"<<(AnsiString(hcp.stackSeq).c_str())<<mendl
	 <<" cross section: "<<((hcp.ccross)?"circ":"hex")<<mendl
	 <<" hcp l.p. a: "<<hcp.a<<mendl
	 <<" hcp l.p. c: "<<hcp.c<<mendl
	 <<" hcp c/a: "<<hcp.c/hcp.a<<mendl
	 <<" hcp Eclp: "<<hcp.elp<<mendl
	 <<" hcp u: "<<hcp.u<<mendl;
	 //o<<mendl<<(*hcp.grainPrm);
	}
	else
	o<<"** hcp stack sequence is empty";

return o;
}
//---------------------------------------------------------------------------
__fastcall TFhcp::TFhcp(TComponent* Owner,StHcpPrm &hcpPrm_) : TForm(Owner),hcpPrm(hcpPrm_),grainPrm(*hcpPrm_.grainPrm)
{
		iniParameters();
		clac=caFree;

		refHcp=hcpPrm.refHcp;

		switch(refHcp){
		case 0 : RBforget->Checked=true;break;
		case 1 : RBab->Checked=true;break;
		case 2 : RBabc->Checked=true;break;
		case 3 : RBabac->Checked=true;break;
		}


//		Ethickness->Enabled=hcpPrm.coreth;

		ChBcutgrain->Checked=hcpPrm.cut;
//
		if(hcpPrm.coreThickness.IsEmpty()){

			Ethickness->Text=float2str(grainPrm.radius-grainPrm.lattConst*0.5);
		}
		else
			Ethickness->Text=hcpPrm.coreThickness;

//		if(hcpPrm.cut){
//
//
//
//		}


}
////---------------------------------------------------------------------------
//__fastcall TFhcp::TFhcp(TComponent* Owner,StGrainPrm *gp,String &dataDir,String *stackSeq): TForm(Owner)
//{
//        manual=false;
//        iniParameters(gp,dataDir,stackSeq);
//}
//---------------------------------------------------------------------------
void TFhcp::iniParameters()
{

		lhcp=grainPrm.lattConst/sqrt(2.0);  //latt. param. hcp latt
		Elhcpa->Text=FloatToStrF(lhcp,ffGeneral,5,5);

		Estacks->Text=hcpPrm.stackTemplate;
		Estacks123->Text=hcpPrm.stackTemplate123;

		if(hcpPrm.store) {
		cdouble ck0=hcpPrm.c/lhcp;

			Eprc->Text=float2str(hcpPrm.u);
			Elhcpc->Text=float2str(hcpPrm.c);
			Eca->Text=hcpPrm.ca;
			Ess->Text=hcpPrm.stackSeq;
			plusMinusSeq();

		}
		else {
			Elhcpc->Text=float2str(grainPrm.lattConst/sqrt(3));
			Eca->Text=float2str(sqrt(2.0/3.0));
			Ess->Text="";
		}


cdouble val=(k2s3*sqr(lhcp)*Elhcpc->Text.ToDouble());
		Eclp->Text=float2str(pow(val,1.0/3));


		if(grainPrm.shape==Eshape::sphere){
		cint maxN=(int)ceil(grainPrm.radius/lhcp/wj);
			 minStackN=2*(maxN+1);

		   Lminsizeval->Caption=IntToStr(minStackN);
		   ptrLatt=lattSphere;
		   ChBcut->Enabled=false;


		   Elayers->Enabled=false;
		}
		else {
			ptrLatt=lattCylinder;
			ChBcut->Checked=hcpPrm.ccross;
			Elayers->Text=IntToStr((int)hcpPrm.numOfLayers);
		}


}
//------------------------------------------------------------------------------
void TFhcp::buildHcp(VGrainNodes *hcpNodes)
{
const String seq=abcSeq.UpperCase();
cint seqSize=seq.Length();


		/// 1. ) program creates the base layer
		lattBase();

		if(grainPrm.shape!=Eshape::sphere)
			grainPrm.maxDist=sqrt(sqr(seqSize*d.z)+sqr(grainPrm.radius*2));

		hcpNodes->clear();

		/// 2. ) program creates whole grain based on ABC sequence and the base layer
		ptrLatt(hcpNodes); //pointer calls the proper method

		if(ChBsubnet->Checked && !ChBadlay->Checked){
		const int nsize=hcpNodes->size();
		const double prc=StrToFloat(Eprc->Text);
		const double dc=d.z*prc;
		const int la=(grainPrm.bilatt)?grainPrm.atomLB:grainPrm.atomLA;
		const String atomName=(grainPrm.bilatt)?grainPrm.atomNameB:grainPrm.atomNameA;

			 if(grainPrm.shape!=Eshape::sphere)
				grainPrm.maxDist=sqrt(sqr(seqSize*d.z+dc)+sqr(grainPrm.radius*2));

			 for(int i=0;i<nsize;i++){
				 grainNode=(*hcpNodes)[i];
				 grainNode.r=(*hcpNodes)[i].r;
                 grainNode.z+=dc;
				 grainNode.atomName=atomName;
				 grainNode.la=la;
				 hcpNodes->push_back(grainNode);
				 hcpNodes->back().layer=(*hcpNodes)[i].layer;
				 hcpNodes->back().centerLayer=false;
			 }
		  }

		clac=caFree;
}
//---------------------------------------------------------------------------
void TFhcp::lattCylinder(VGrainNodes *hcpNodes)
{
int i,j;
         if(ChBcut->Checked){
         const int bsize=base.size();
         cdouble gpR2=sqr(grainPrm.radius);
         double x,y,r2;
         vector<StGrainNode> cpbase;

                for(i=0;i<bsize;i++){
                    x=(base)[i].x;
                    y=(base)[i].y;
                    r2=sqr(x)+sqr(y);

                    if(r2<=gpR2)
                      cpbase.push_back((base)[i]);
                }
               base=cpbase;
         }

         /////////////////////////////////

cdouble ht=lhcp*sqrt(3)/2;

         a.x=0;
         a.y=0;
         a.z=0;

         b.x=ht/3; // shift 1/3
         b.y=lhcp/2.0;
         b.z=0;

         c.x=ht*2/3; //shift 2/3
         c.y=0;
         c.z=0;

const int bsize=base.size();
const String seq=abcSeq.UpperCase();
cint seqSize=seq.Length();
cint seqSize2=(seqSize+1)/2;
char lett;
StVector vd,vres;
StGrainNode tnode;
double r;

		 hcpNodes->reserve(bsize*seqSize);

		 for(i=0;i<seqSize;i++){
             lett=seq[i+1];

             vd=d*i;

			switch(lett){
			case 'A': vres=a+vd;break;
			case 'B': vres=b+vd;break;
			case 'C': vres=c+vd;break;
			default: msgDlgErr("Unknown letter. STOP!"); return;
            }

            for(j=0;j<bsize;j++){
				tnode=base[j];
				tnode+=vres;
				tnode.layer=i;
                r=sqrt( sqr(tnode.x)+sqr(tnode.y) );
				hcpNodes->push_back(tnode);
				hcpNodes->back().centerLayer=(seqSize2==i);
            }
         }
}
//---------------------------------------------------------------------------
void TFhcp::lattSphere(VGrainNodes *hcpNodes)
{
cdouble ht=lhcp*0.5*sqrt(3);

			 a.x=0;
			 a.y=0;
			 a.z=0;

			 b.x=ht/3; // shift 1/3
			 b.y=lhcp/2.0;
			 b.z=0;

			 c.x=ht*2/3; //shift 2/3
			 c.y=0;
			 c.z=0;

const int bsize=base.size();
const String seq=abcSeq.UpperCase();
cint seqSize=seq.Length();
cint seqSize2=(seqSize+2)/2;
cdouble gpR2=sqr(grainPrm.radius);
StVector vz,vres;

			if(ChBsubnet->Checked)
				hcpNodes->reserve(bsize*seqSize*2);
			else
				hcpNodes->reserve(bsize*seqSize);

StGrainNode tnodeP;
int i,j,k;
double r2,ii;

			if(ChBadlay->Checked){
			const double prc=StrToFloat(Eprc->Text);
			const double dc=d.z*prc;
			const String atomName=(grainPrm.bilatt)?grainPrm.atomNameB:grainPrm.atomNameA;

					for(ii=-seqSize>>1,k=1;k<=seqSize;ii++,k++){
						vz=d*ii;

						switch(seq[k]){
						case 'A': vres=a+vz;break;
						case 'B': vres=b+vz;break;
						case 'C': vres=c+vz;break;
						default: msgDlgErr("Unknown letter. STOP!"); return;
						}

						for(j=0;j<bsize;j++){
							tnodeP=base[j];
							tnodeP+=vres;
							r2=sqr(tnodeP.x)+sqr(tnodeP.y)+sqr(tnodeP.z);

							if(r2<=gpR2)  {
							   tnodeP.r=sqrt(r2);
							   hcpNodes->push_back(tnodeP);
							   hcpNodes->back().layer=k;
							}
						}
						/// sublattice

						vres.z+=dc;

						for(j=0;j<bsize;j++){
							tnodeP=base[j];
							tnodeP+=vres;
							r2=sqr(tnodeP.x)+sqr(tnodeP.y)+sqr(tnodeP.z);

							if(r2<=gpR2)  {
							   tnodeP.r=sqrt(r2);
							   hcpNodes->push_back(tnodeP);
							   hcpNodes->back().layer=k;
							   hcpNodes->back().atomName=atomName;
   							   hcpNodes->back().centerLayer=(seqSize2==k);
							}
						}

					}
			}
			else{
					for(ii=-seqSize>>1,k=1;k<=seqSize;ii++,k++){
						vz=d*ii;

						switch(seq[k]){
						case 'A': vres=a+vz;break;
						case 'B': vres=b+vz;break;
						case 'C': vres=c+vz;break;
						default: msgDlgErr("Unknown letter. STOP!"); return;
						}

						for(j=0;j<bsize;j++){
							tnodeP=base[j];
							tnodeP+=vres;
							r2=sqr(tnodeP.x)+sqr(tnodeP.y)+sqr(tnodeP.z);

							if(r2<=gpR2)  {
							   tnodeP.r=sqrt(r2);
							   //tnodeP.layer=k;
							   hcpNodes->push_back(tnodeP);
							   hcpNodes->back().layer=k;
							   hcpNodes->back().centerLayer=(seqSize2==k);

							 //  Memo->Lines->Add(String(hcpNodes->back().layer));
							}
						}
					}
			}

}
//---------------------------------------------------------------------------
void __fastcall TFhcp::lattBase()
{
		lhcp=Elhcpa->Text.ToDouble();

                 /*
		//k¹t nachylenia wersorów a, b +/-60 do osi X
        a.x=lhcp*wi;
		a.y=lhcp*wj;
        a.z=0;

		b.x=lhcp*wi;
        b.y=-lhcp*wj;
		b.z=0;
        */

		//k¹t nachylenia wersora a:  alfa=90, kat wersora b: beta=-30 do osi X
        a.x=0;
		a.y=lhcp;
		a.z=0;

		b.x=lhcp*wj;
		b.y=-lhcp*wi;
		b.z=0;

        c=a+b;

cdouble hcpC=Elhcpc->Text.ToDouble();

		d.x=0;
        d.y=0;
		d.z=hcpC;

		base.clear();

cint maxN=(int)ceil(grainPrm.radius/lhcp/wj);// a few atoms will be generated outside of circle
StVector va,vb,vd,vres;

		 if(ChBsubnet->Checked)
			base.reserve(6*maxN*(maxN-1)+2);
		 else
			base.reserve(3*maxN*(maxN-1)+1);

		grainNode.atomName=grainPrm.atomNameA;
		grainNode.la=grainPrm.atomLA;
		grainNode.rprev=0;

int l;
         //generuje bazowa siec : base
		 for(l=-maxN;l<=maxN;l++){   //wzdluz przekatnej
             va=a*l;

			 grainNode=va;
			 grainNode.r=va.getModule();

			 base.push_back(grainNode);
		 }

		 vd=a*maxN;

		 for(l=1;l<=maxN;l++){
			 va=vd*(-1)+c*l;
			 vb=vd-c*l;

			 for(int p=2*maxN-l+1;p>0;p--){

				 grainNode=va;
				 grainNode.r=va.getModule();

				 base.push_back(grainNode);

				 grainNode=vb;
				 grainNode.r=vb.getModule();

				 base.push_back(grainNode);
				 va=va+a;
				 vb=vb-a;
			  }
		 }
}

//---------------------------------------------------------------------------
void __fastcall TFhcp::EssChange(TObject *Sender)
{
cint tlength=Ess->Text.Length();
		Lstacksize->Caption=IntToStr(tlength);

		if(grainPrm.shape==Eshape::sphere){
		   if(Ess->Text.Length()>minStackN) {
			  Ess->Font->Color=clRed;
              Epm->Font->Color=clRed;
              Lstacksize->Font->Color=clRed;
           }

           if(tlength==minStackN){
			  Ess->Font->Color=clGreen;
              Epm->Font->Color=clGreen;
              Lstacksize->Font->Color=clGreen;
           }


           if(tlength<minStackN){
              Ess->Font->Color=clBlack;
              Epm->Font->Color=clBlack;
			  Lstacksize->Font->Color=clBlack;
           }
        }
        else{
        cdouble lpc=Elhcpc->Text.ToDouble();
        cdouble res=lpc*tlength;
                Lminsizeval->Caption=float2str(res);
        }

		plusMinusSeq();

		/// near center ABC seq
		//const int ncSeqLength=(tlength>5)? 4+tlength%2 : tlength;

		Labcnearcenter->Tag=(tlength>5)? tlength/2-1:0 ;
		String ncSeq=(tlength>5)?Ess->Text.SubString(Labcnearcenter->Tag,4+tlength%2):Ess->Text;

		Labcnearcenter->Caption=ncSeq;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::EssKeyPress(TObject *Sender, char &Key)
{
		if(Key=='a' || Key=='b' || Key=='c')
		   Key-=32;

		if(Key!='A' && Key !='B' && Key !='C' && Key!=8){
		   Beep(200,200);
		   Key=0;
		   return;
		}

		if(grainPrm.shape==Eshape::sphere){
        cint tlength=Ess->Text.Length();

           if(tlength==minStackN && Key!=8){
              Beep(500,200);
              Key=0;
           }
           if(tlength>minStackN)
              Key=8;  //backspace

        }
}
//---------------------------------------------------------------------------
void __fastcall TFhcp::EstacksKeyPress(TObject *Sender, System::WideChar &Key)
{
		if(Key=='a' || Key=='b' || Key=='c')
		   Key-=32;

		if(Key!='A' && Key !='B' && Key !='C' && Key!=8 && Key!='(' && Key!=')' && !IsDigit(Key)){
		   Beep(200,200);
		   Key=0;
		   return;
		}
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::EcaChange(TObject *Sender)
{
cint tlength=Eca->Text.Length();
		if(tlength){
		cdouble ck0=Eca->Text.ToDouble()*lhcp;
		cdouble val=ck0*k2s3*sqr(lhcp);

		   Eclp->Text=float2str(pow(val,1.0/3));

		   Elhcpc->OnChange=NULL;
		   Elhcpc->Text=float2str(ck0);
		   Elhcpc->OnChange=ElhcpcChange;
		}
}
//---------------------------------------------------------------------------
void __fastcall TFhcp::ElhcpcChange(TObject *Sender)
{
cint tlength=Elhcpc->Text.Length();
        if(tlength){
		cdouble ck0=Elhcpc->Text.ToDouble()/lhcp;

           Eca->OnChange=NULL;
           Eca->Text=float2str(ck0);
		   Eca->OnChange=EcaChange;

		cdouble val=(k2s3*sqr(lhcp)*Elhcpc->Text.ToDouble());
			Eclp->Text=float2str(pow(val,1.0/3));
		}
}
//---------------------------------------------------------------------------
void TFhcp::plusMinusSeq()
{
const String seq=Ess->Text.UpperCase();
cint seqLen=seq.Length();

        if(seqLen<2){
           Epm->Text="";
           return;
        }

String pmSeq;

int i;
char lettP,lettC;
int lp,lc;
int cond;

		lp=(int) seq[1];
        for(i=2;i<=seqLen;i++){

			lc=(int) seq[i];
			cond=lc-lp;

            if( cond ==1 || cond ==-2)
                pmSeq+="+";
            else
               if(cond==0)
                  pmSeq+="O";
               else
                  pmSeq+="-";

			lp=lc;
        }
        Epm->Text=pmSeq;
}
//---------------------------------------------------------------------------
void __fastcall TFhcp::ChBcutClick(TObject *Sender)
{
		hcpPrm.ccross=ChBcut->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::BokhexClick(TObject *Sender)
{

		if(grainPrm.shape==Eshape::sphere && Ess->Text.Length()>Lminsizeval->Caption.ToInt()){
			   //	sst.token=sst.token.SubString(1,Lminsizeval->Caption.ToInt());
			   Ess->Text=Ess->Text.SubString(1,Lminsizeval->Caption.ToInt());
		}

		hcpPrm.store=true;
		hcpPrm.a=StrToFloat(Elhcpa->Text);
		hcpPrm.c=StrToFloat(Elhcpc->Text);
		hcpPrm.elp=StrToFloat(Eclp->Text);
		hcpPrm.u=StrToFloat(Eprc->Text);
		hcpPrm.stackTemplate=Estacks->Text;
		hcpPrm.refHcp=refHcp;
		hcpPrm.ca=Eca->Text;
		hcpPrm.stackTemplate123=Estacks123->Text;
		hcpPrm.coreth=RBsurfatoms->Checked;
		RBcutbottomtop->Checked=!RBsurfatoms->Checked;//ChBcutbottomtop->Checked;
		hcpPrm.coreThickness=Ethickness->Text;
		hcpPrm.cut=ChBcutgrain->Checked;


		if(refHcp){
		String idealStr;

			switch(refHcp){
			case 1: idealSeq("AB",idealStr); break;
			case 2: idealSeq("ABC",idealStr);break;
			case 3: idealSeq("ABAC",idealStr);break;
			}
//
			if(Ess->Text.Length()>=idealStr.Length())
				lattSphereWithRefHcp(Ess->Text,idealStr);
			else{
				clac=caNone;
				msgDlgErr("Model sequence is too short. It should be longer than ideal structure sequence");
			}
		}
		else{
			abcSeq=Ess->Text;
			buildHcp(hcpPrm.nodesHcp);
		}

		hcpPrm.stackSeq=Ess->Text;

		if(ChBcutgrain->Checked){
		bool res=false;
				if(RBsurfatoms->Checked)
				   res=cutSphereGrain();
				else
                	res=cutSphereTopBottomLayers();


			 if(!res ) clac=caNone;

		}

		adjustCenterAtom();

}
//---------------------------------------------------------------------------

void __fastcall TFhcp::BabcClick(TObject *Sender)
{
TButton *lab=dynamic_cast<TButton *> (Sender);
String seq,result;

			switch(lab->Tag){
			case 0: seq="ABC";break;
			case 1: seq="AB";break;
			case 2: seq="ABAC";break;
			case 3: seq="AC";break;
			}

			idealSeq(seq,result);
			Ess->Text=result;
}
//---------------------------------------------------------------------------
void TFhcp::idealSeq(const String seq,String &result)
{
String ss;
unsigned seqLen,i;
int   minStackN=1;

			if(hcpPrm.grainPrm->shape==Eshape::sphere){
			cint maxN=(int)ceil(grainPrm.radius/lhcp/wj);
				minStackN=2*(maxN+1);
			}
			else{
				try{

				minStackN=Elayers->Text.ToInt();
				hcpPrm.numOfLayers=minStackN;
				}
				catch(...){
					msgDlgErr("Number of layers should be integer value greater than 0");
					return;
				}
			}

			for(i=0;i<minStackN;i++)
				ss+=seq[i%seq.Length()+1];

				result=ss;

}

//---------------------------------------------------------------------------
StExpr TFhcp::buildABC(const unsigned startPos=1)
{
const String &str=expression;
size_t pos=startPos;
StExpr expr;

		expr.start=pos;
		pos++;

		while( pos<exprLength && str[pos]!=')'){

			if(IsLetter(str[pos]))
				expr.token+=str[pos];
			else{
				if(IsDigit(str[pos])){
				StExpr exprUp=buildABCexpand(pos);
					pos+=exprUp.Length();
					expr.token+=exprUp.token;
				}
				else{
				StExpr exprUp=buildABC(pos);
					pos+=exprUp.Length();
					expr.token+=exprUp.token;
				}
			}
			pos++;
		}

//		if(pos==strLen && str[pos]!=')'  )
//		throw TExprException(TExprException::EXPRESSION,pos);
		if(pos>exprLength)
		throw TExprException(TExprException::EXPRESSION,pos);

		expr.stop=pos;

return	expr;
}
//---------------------------------------------------------------------------
StExpr TFhcp::buildABCexpand(const unsigned startPos)
{
const String &str=expression;
size_t pos=startPos;
StExpr expr;
String digits;

		expr.start=startPos;

		while(pos<exprLength && IsDigit(str[pos]) )
			digits+=str[pos++];

String letters;

		if(IsLetter(str[pos]))
			letters=str[pos];
		else
			if(str[pos]=='(') {
			StExpr exprUp=buildABC(pos);
					pos+=exprUp.Length();
					letters=exprUp.token;

			}
			else
			throw TExprException(TExprException::EXPRESSION,pos);

		expr.stop=pos;

const size_t intNum=StrToInt(digits);

			expr.token="";
			for(size_t i=0;i<intNum;i++)
				expr.token+=letters;

return expr;
}
//---------------------------------------------------------------------------
void __fastcall TFhcp::BtemplateClick(TObject *Sender)
{
		try{
				expression="("+Estacks->Text+")";
				exprLength=Estacks->Text.Length()+2;
		StExpr sst=buildABC();

				if( sst.stop != Estacks->Text.Length()+2)
					throw  TExprException(TExprException::POSITION,Estacks->Text.Length());

				if(grainPrm.shape==Eshape::sphere){  //obcina do odpowiedniej dlugosci
					sst.token=sst.token.SubString(1,Lminsizeval->Caption.ToInt());
				}
				Ess->Text=sst.token;

		}
		catch(TExprException &e){
				msgDlgErr("syntax error: "+String(e.what())+", at position: "+String(e.getPosition()));

		}
		catch(...){
				msgDlgErr("Unknown error");
		}
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::RBforgetClick(TObject *Sender)
{
TRadioButton *rb=dynamic_cast <TRadioButton *> (Sender);

		refHcp=rb->Tag;
}

//---------------------------------------------------------------------------
void __fastcall TFhcp::lattSphereWithRefHcp(const String &genHcpSeq,const String &refHcpSeq)
{
const size_t genHcpSize=genHcpSeq.Length();
const size_t refHcpSize=refHcpSeq.Length();
VGrainNodes * const genHcpNodes=hcpPrm.nodesHcp;
VGrainNodes * const refHcpNodes=hcpPrm.nodesHcpRef;

			lattBase();

const size_t bsize=base.size();
const String idealSeq=refHcpSeq.UpperCase();
const String genSeq=genHcpSeq.UpperCase();
int seqSize=idealSeq.Length();
cdouble gpR2=sqr(grainPrm.radius);

			genHcpNodes->clear();
			refHcpNodes->clear();

			genHcpNodes->reserve(bsize*seqSize);
			refHcpNodes->reserve(bsize*seqSize);


cdouble ht=lhcp*sqrt(3)/2;
            //// wektory translacji struktur: A, B, C
			a.x=0;
			a.y=0;
			a.z=0;

			b.x=ht/3; // shift 1/3
			b.y=lhcp/2.0;
			b.z=0;

			c.x=ht*2/3; //shift 2/3
			c.y=0;
		 	c.z=0;

int ii,j,k;
double r2;
StGrainNode tnodeP,tnodeI,tnodeG;
StVector vz,vres;


			for(ii=-seqSize/2.0,k=1;k<=seqSize;ii++,k++){

				vz=d*ii;

				switch(idealSeq[k]){
				case 'A': vres=a+vz;break;
				case 'B': vres=b+vz;break;
				case 'C': vres=c+vz;break;
				default: msgDlgErr("Unknown letter. STOP!"); return;
				}

				for(j=0;j<bsize;j++){
					tnodeP=base[j];
					tnodeI=tnodeP+vres;
					r2=sqr(tnodeI.x)+sqr(tnodeI.y)+sqr(tnodeI.z);


					if(r2<=gpR2)  { /// for sphere condition
						tnodeP.r=sqrt(r2);
						refHcpNodes->push_back(tnodeI);

						switch(genSeq[k]){
						case 'A': vres=a+vz;break;
						case 'B': vres=b+vz;break;
						case 'C': vres=c+vz;break;
						default: msgDlgErr("Unknown letter. STOP!"); return;
						}

						tnodeG=tnodeP+vres;
						r2=sqr(tnodeG.x)+sqr(tnodeG.y)+sqr(tnodeG.z);
						tnodeG.r=sqrt(r2);
						genHcpNodes->push_back(tnodeG);

					}
				}
			}


		//if(ChBsubnet->Checked){
		const int nsize=refHcpNodes->size();
		const double prc=StrToFloat(Eprc->Text);
        const double dc=d.z*prc;
		const int la=(grainPrm.bilatt)?grainPrm.atomLB:grainPrm.atomLA;
        const String atomName=(grainPrm.bilatt)?grainPrm.atomNameB:grainPrm.atomNameA;

			 if(grainPrm.shape!=Eshape::sphere)
				grainPrm.maxDist=sqrt(sqr(seqSize*d.z+dc)+sqr(grainPrm.radius*2));

			 for(int i=0;i<nsize;i++){
				 grainNode=(*refHcpNodes)[i];
				 grainNode.r=(*refHcpNodes)[i].r;
                 grainNode.z+=dc;
				 grainNode.atomName=atomName;
				 grainNode.la=la;
				 refHcpNodes->push_back(grainNode);

				 grainNode=(*genHcpNodes)[i];
				 grainNode.r=(*genHcpNodes)[i].r;
                 grainNode.z+=dc;
				 grainNode.atomName=atomName;
				 grainNode.la=la;
				 genHcpNodes->push_back(grainNode);


			 }
		//  }

}
//---------------------------------------------------------------------------
void TFhcp::loadFile(const String &fileName)
{
wfstream file(fileName.c_str(),ios::in);

			if(!file){
			   MessageDlg("Error, couldn't open file: "+SDial->FileName+" for writing",mtError, TMsgDlgButtons() << mbOK,0);
			return;
            }

const int oneLineSize=128;
wchar_t oneLine[oneLineSize-1];
wchar_t cmd[16];

        //czas
			file.getline(oneLine,oneLineSize,'\n');

unsigned size;
            file>>cmd>>size;
			file.ignore(1);

wchar_t *seq=new wchar_t[size+2];

			file.getline(seq,size+1);
            seq[size]=0;
			Ess->Text=String(seq);

			file.close();
delete [] seq;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::MIseqtofileClick(TObject *Sender)
{
		SDial->InitialDir=dataDir;

		if(SDial->Execute()){
		wfstream file(SDial->FileName.c_str(),ios::out);

                if(!file){
				   MessageDlg("Error, couldn't open file: "+SDial->FileName+" for writing",mtError, TMsgDlgButtons() << mbOK,0);
				return;
				}
		TDateTime today=Now();

				file<<"#date: "<<DateToStr(today).c_str()<<", "<<TimeToStr(today).c_str()<<endl;
				file<<"#size: "<<Ess->Text.Length()<<endl;
				file<<Ess->Text.c_str();

				file.close();
		}
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::ChBsubnetClick(TObject *Sender)
{
		Eprc->Enabled=ChBsubnet->Checked;
}

//---------------------------------------------------------------------------
void __fastcall TFhcp::FormClose(TObject *Sender, TCloseAction &Action)
{
		Action=clac;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::BitBtn3Click(TObject *Sender)
{
		this->ModalResult=mrCancel;
		clac=caFree;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::MIseqfromfileClick(TObject *Sender)
{
		ODial->InitialDir=dataDir;

		if(ODial->Execute()){
		   loadFile(ODial->FileName);
		}
}
//---------------------------------------------------------------------------

String getTokenPM(const String &str,int &start,const int &strLength,size_t &tokenLength)
{
int i;
String token;

		i=start;
		tokenLength=0;
		while(i<=strLength && IsDigit(str,i) ){
			i++;
			tokenLength++;
		}

	token=str.SubString(start,tokenLength);
	start+=tokenLength+1;

return token;
}

//---------------------------------------------------------------------------
void __fastcall TFhcp::Etemplate123Click(TObject *Sender)
{
//
int start=1;
String token,result("A");
const String seq=Estacks123->Text;

const size_t seqLength=seq.Length();
size_t tokenLength;
bool pm=true;

const String abc("ABC");
const String cba("CBA");
const String *ptrAbc=&abc;
size_t abcPos=1;

wchar_t lastLetter;

		while(start<=seqLength){

			lastLetter=(*result.LastChar());
			ptrAbc=(pm)?&abc:&cba;

			token=getTokenPM(seq,start,seqLength,tokenLength);

			if(tokenLength==0) continue;

			const int stackSize=token.ToInt();

			if(pm){ //plusy
					switch(lastLetter){
					case 'A' : abcPos=1;break;
					case 'B' : abcPos=2;break;
					case 'C' : abcPos=3;break;
					}
			}
			else{
					switch(lastLetter){
					case 'A' : abcPos=3;break;
					case 'B' : abcPos=2;break;
					case 'C' : abcPos=1;break;
					}

			}

				for(size_t i=0;i<stackSize;i++){
					abcPos=abcPos%3+1;
					result+=(*ptrAbc)[abcPos];
				}

			pm=!pm;
		}

		Ess->Text=result;

}

//---------------------------------------------------------------------------

void __fastcall TFhcp::Estacks123KeyPress(TObject *Sender, System::WideChar &Key)

{
	if(! IsDigit(Key) && (Key!=',') && Key!=8)
			Key=0;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::adjustCenterAtom()
{
VGrainNodes * const genHcpNodes=hcpPrm.nodesHcp;
const int hcpSize=genHcpNodes->size();
double rd=(*genHcpNodes)[0].r;


			if(rd==0)
				return;

int i;
int savePos=0;

			/// search for atom closest to the center 0,0,0 point
			for(i=0;i<hcpSize;i++)
				if( (*genHcpNodes)[i].r<rd){
					savePos=i;
					rd=(*genHcpNodes)[i].r;
				}

			StGrainNode::centerGNindex=savePos;


			if(rd==0)
				return;

StGrainNode dispVect=(*genHcpNodes)[savePos];

double x,y,z;
			for(i=hcpSize-1;i>=0;i--){
				(*genHcpNodes)[i]-=dispVect;

				x=(*genHcpNodes)[i].x;
				y=(*genHcpNodes)[i].y;
				z=(*genHcpNodes)[i].z;
				rd=x*x+y*y+z*z;
				(*genHcpNodes)[i].r=sqrt(rd);
			}

}
//---------------------------------------------------------------------------
void __fastcall TFhcp::cmSumXYZ(const iterVGN &start,const size_t N, double &sumX,double &sumY,double &sumZ)
{
const iterVGN stop=start+N;

		sumX=sumY=sumZ=0;

		for(auto iter=start;iter!=stop;iter++){
			sumX+=iter->x;
			sumY+=iter->y;
			sumZ+=iter->z;
		}


}
//---------------------------------------------------------------------------
bool __fastcall TFhcp::cutSphereGrain()
{
VGrainNodes * const genHcpNodes=hcpPrm.nodesHcp;

iterVGN iterLatt0=genHcpNodes->begin();//lattice ref
csize_t hcpSize=genHcpNodes->size();
csize_t hcpSize2=hcpSize>>1;
cdouble ihcpSize=1.0/hcpSize;
auto iterLatt1=iterLatt0+(hcpSize/2);  //subblatice
auto iterEnd=genHcpNodes->end();		// the end

double cmX,cmY,cmZ;
			// 1)  center of mass searching

			//IMPORTANT ASSUMPTION: all atoms have the same mass


		   cmSumXYZ(iterLatt0,hcpSize,cmX,cmY,cmZ);
		   cmX*=ihcpSize;
		   cmY*=ihcpSize;
		   cmZ*=ihcpSize;

cdouble r2core=sqr(Ethickness->Text.ToDouble());
//cdouble rcore=grainPrm
double r2;
size_t i,j;

VGrainNodes tmpnodesIN,tmpnodesOUT;

			tmpnodesIN.reserve(hcpSize);
			tmpnodesOUT.reserve(hcpSize);

			j=0;
			for(auto iter=iterLatt0;iter!=iterEnd;iter++) {
				r2=sqr(iter->x-cmX)+sqr(iter->y-cmY)+sqr(iter->z-cmZ);

				if(r2<=r2core){
					tmpnodesIN.push_back(*iter);
					j++;
				}
				else{
					tmpnodesOUT.push_back(*iter);
					tmpnodesOUT.back().r=r2;
				}
            }

			//Memo->Lines->Add("atoms out: "+String(tmpnodesOUT.size()));


			i=j=0;
			for(auto iter=tmpnodesOUT.begin();iter!=tmpnodesOUT.end();iter++){

					if(grainPrm.atomNameA==iter->atomName)
						i++;
					else
						j++;
			}


			if(i!=j){
				msgDlgErr(" i!=j   "+String(i)+", "+String(j));


			   if(i>j){
				   sort (tmpnodesOUT.begin(),tmpnodesOUT.begin()+i+1,sortByRadii);

			   auto iterOUT=tmpnodesOUT.begin();

				   while(i>j && i>0){
						tmpnodesIN.push_back((*iterOUT));
						iterOUT++;
						i--;
				   }
			   }
			   if(j>i){
					sort(tmpnodesOUT.begin()+i+1,tmpnodesOUT.end(),sortByRadii);

			   auto iterOUT=tmpnodesOUT.begin()+i+1;

					while(j>i && j>0){
						tmpnodesIN.push_back(*iterOUT);
						iterOUT++;
						j--;
					}
			   }
			//return false;
			}


			genHcpNodes->assign(tmpnodesIN.begin(),tmpnodesIN.end());
			genHcpNodes->resize(tmpnodesIN.size());

			sort(genHcpNodes->begin(),genHcpNodes->end(),sortByAtomName);

return true;
}

//---------------------------------------------------------------------------
bool __fastcall TFhcp::cutSphereTopBottomLayers()
{
VGrainNodes * const genHcpNodes=hcpPrm.nodesHcp;

int nodesToEraseBottom=0;
auto iter=genHcpNodes->begin();
auto iterEnd=genHcpNodes->end();
int refLayer=iter->layer;

		do{
			nodesToEraseBottom++;
			iter++;
		} while( iter->layer == refLayer && iter!=iterEnd);


		genHcpNodes->erase(genHcpNodes->begin(),genHcpNodes->begin()+nodesToEraseBottom);

auto iterBegin=genHcpNodes->begin();
int nodesToEraseTop=0;

		iter=genHcpNodes->end()-1;

		refLayer=iter->layer;
		nodesToEraseTop=0;

		do{
			nodesToEraseTop++;
			iter--;
		}while(iter->layer==refLayer && iter!=iterBegin);

		genHcpNodes->erase(genHcpNodes->end()-nodesToEraseTop, genHcpNodes->end());

		if(nodesToEraseBottom != nodesToEraseTop){
			msgDlgErr(" bottom != top   "+String(nodesToEraseBottom)+", "+String(nodesToEraseTop));
		}

return true;
}

//---------------------------------------------------------------------------
void __fastcall TFhcp::LabcnearcenterClick(TObject *Sender)
{
const int selStart=Labcnearcenter->Tag-1;
const int selLength=Labcnearcenter->Caption.Length();

			Ess->SetFocus();
			Ess->SelStart=selStart;
			Ess->SelLength=selLength;
}
//---------------------------------------------------------------------------
void __fastcall TFhcp::findCenterAtomsInLayers()
{
csize_t numOfLayers=Ess->Text.Length();
size_t i,j;



			for(i=0;i<numOfLayers;i++){




			}
}
//---------------------------------------------------------------------------
void __fastcall TFhcp::ChBcutgrainClick(TObject *Sender)
{
	 //	Ethickness->Enabled=ChBcutbottomtop->Checked;
}
//---------------------------------------------------------------------------


void __fastcall TFhcp::RBsurfatomsClick(TObject *Sender)
{
		Ethickness->Enabled=RBsurfatoms->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TFhcp::markCenterLayerAtoms()
{

}
