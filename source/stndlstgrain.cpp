

#include "stndlstgrain.h"
#include "elements.h"
#include <cmath>

const int headerVersion=0;
const int singleLineSize=63;

//#define debug
#include "clcmdlist.h"

//---------------------------------------------------------------------------
wifstream& operator >> (wifstream &in, StFileHeaderNDL &fileHeader)
{
wchar_t cmd[32],singleLine[singleLineSize+1];
std::wstring wcmd;
std::wstring wstr_date(L"#date:");
StGrainPrm &grainPrm=*(fileHeader.ptrGrainPrm);

        singleLine[singleLineSize]=0;
       ///czas

        do{  in>>wcmd; }
        while(wcmd.find(wstr_date)!=0 && !in.eof());

        in.getline(singleLine,singleLineSize,'\n');

        #ifdef debug
        wcout<<" header: time "<<singleLine<<endl;
        #endif

       ///rozmiar
std::wstring sizerc;
std::wstring fs(L"sizeRC");

        in>>sizerc;

        #ifdef debug
        wcout<<" header: sizerc "<<sizerc<<endl;
        #endif


int fff=sizerc.find(fs);

        if(fff>0)
            in>>fileHeader.dataRows>>fileHeader.dataCols;
        else
            in>>fileHeader.dataRows;

        ///#lattice parameter
        in>>cmd>>grainPrm.lattConst;

        ///#promien ziarna
        in>>cmd>>grainPrm.radius;

        ///struktura
        in>>cmd;
        in.ignore(1);
        in.getline(singleLine,singleLineSize,'\n');

        ///sphere/cylinder
int gshape;
        in>>cmd>>gshape;


        ///modyfikowana
        in>>cmd>>grainPrm.modified;

        ///nazwy atomów
wstring atomTypes;

        in>>cmd>>atomTypes;

vector<wstring> atypes(split<wstring> (atomTypes,L" "));


       grainPrm.atomsNameNumber.reserve(2);

       for(auto &atype: atypes)
            grainPrm.atomsNameNumber.push_back(atype);


        grainPrm.atomsNameNumber.shrink_to_fit();


        for(size_t i=0;i<grainPrm.atomsNameNumber.size();i++){
        const wstring amass(Elements::mass.at(grainPrm.atomsNameNumber[i].name));

                switch (i) {
                case 0: grainPrm.atomMassA=( amass.empty() ) ? L"?" : amass; grainPrm.bilatt=false; break;
                case 1: grainPrm.atomMassB=( amass.empty() ) ? L"?" : amass; grainPrm.bilatt=true; break;
                }
        }



        in>>cmd>>grainPrm.numOfatoms;

        #ifdef debug
        cout<<" header: numOfatoms"<<grainPrm.numOfatoms<<endl;
        #endif

return in;
}

//---------------------------------------------------------------------------

//StGrainPrm::eScFactors StGrainPrm::getAtomMass(const QString &aname, QString &amass)
//{
//const QString fileName("scFact.sft");

//            if(!QFile(fileName).exists()){
//                cerr<<" couldn't find scFact.sft"<<endl;
//                return SF_FILE_ERROR;
//            }

//QSettings fileSF(fileName,QSettings::IniFormat);
//const QString keyName(aname.toUpper()+"/mass"); //capital letters


//            amass=fileSF.value(keyName).toString();

//            if(amass.isEmpty())  {
//                cerr<<" empty mass value of : "<<aname.toStdString()<<endl;
//                return SF_NO_ITEM_EXISTS;
//            }

//            //cerr<<" mass "<< amass.toStdString()<<endl;

//return SF_OK;
//}


//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
coordinate sqrCoord(const coordinate &x)
{
return x*x;
}


/**
 * @brief StRotationMatrix::buildMatrix
 * @param axis_
 * ///https://en.wikipedia.org/wiki/Rotation_matrix
 */
void StRotationMatrix::buildMatrix(const StAxis &axis_)
{
            buildRotationAxis(axis_);

ccoord cosTheta=cos(theta);
ccoord sinTheta=sin(theta);
ccoord oneMinCos=1-cosTheta;

        m11=cosTheta+sqrCoord(ux)*oneMinCos; m12=ux*uy*oneMinCos-uz*sinTheta; m13=ux*uz*oneMinCos+uy*sinTheta;
        m21=uy*ux*oneMinCos+uz*sinTheta; m22=cosTheta+sqrCoord(uy)*oneMinCos; m23=uy*uz*oneMinCos-ux*sinTheta;
        m31=uz*ux*oneMinCos-uy*sinTheta; m32=uz*uy*oneMinCos+ux*sinTheta; m33=cosTheta+sqrCoord(uz)*oneMinCos;

        on=true;
}
//-----------------------------------------------------------------------------
void StRotationMatrix::buildMatrix(const StAxis &axis_, const coordinate &angle)
{
        theta=angle;
        buildRotationAxis(axis_);

ccoord cosTheta=cos(theta);
ccoord sinTheta=sin(theta);//sqrt(1-sqrCoord(cosTheta));
ccoord oneMinCos=1-cosTheta;

        m11=cosTheta+sqrCoord(ux)*oneMinCos; m12=ux*uy*oneMinCos-uz*sinTheta; m13=ux*uz*oneMinCos+uy*sinTheta;
        m21=uy*ux*oneMinCos+uz*sinTheta; m22=cosTheta+sqrCoord(uy)*oneMinCos; m23=uy*uz*oneMinCos-ux*sinTheta;
        m31=uz*ux*oneMinCos-uy*sinTheta; m32=uz*uy*oneMinCos+ux*sinTheta; m33=cosTheta+sqrCoord(uz)*oneMinCos;

        on=true;
}

//-----------------------------------------------------------------------------
void StRotationMatrix::buildRotationAxis(const StAxis &axis_)
{
        ux=axis_.a;
        uy=axis_.b;
        uz=axis_.c;

        normUxyz();
}

//-----------------------------------------------------------------------------
void StRotationMatrix::normUxyz()
{
ccoord sumSq=sqrCoord(ux)+sqrCoord(uy)+sqrCoord(uz);
ccoord norm=sqrt(1/sumSq);

            ux*=norm;
            uy*=norm;
            uz*=norm;

}
//-----------------------------------------------------------------------------
StVector3 StRotationMatrix::operator*(const StVector3 &a)
{
 coordinate bx,by,bz;

                bx=m11*a.x+m12*a.y+m13*a.z;
                by=m21*a.x+m22*a.y+m23*a.z;
                bz=m31*a.x+m32*a.y+m33*a.z;

return StVector3(bx,by,bz);
}

wostream & operator <<(wostream &o, StGrainAtom &atom)
{
    o<<atom.name<<"   "<<atom.x<<"    "<<atom.y<<"    "<<atom.z;
return o;
}
