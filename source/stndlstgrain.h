#ifndef STNDLSTGRAIN
#define STNDLSTGRAIN

//#include <QFile>
//#include <QString>
//#include <QSettings>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>


//#include "elements.h"

using namespace std;
typedef float position;
typedef const position cpos;
typedef position coordinate;
typedef const coordinate ccoord;

/*template<typename T>
vector<T>  split(const T & str, const T & delimiters)
{
vector<T> v;
typename T::size_type start = 0;
auto pos = str.find_first_of(delimiters, start);

        while(pos != T::npos) {
            if(pos != start) // ignore empty tokens
                v.emplace_back(str, start, pos - start);
            start = pos + 1;
            pos = str.find_first_of(delimiters, start);
        }
        if(start < str.length()) // ignore trailing delimiter
            v.emplace_back(str, start, str.length() - start); // add what's left of the string
return v;
}*/

//-----------------------------------------------------------------------------
struct StVector3{
    union{
    struct{coordinate x,y,z;};
    coordinate xyz[3];
    };

   StVector3(){ x=y=z=0;}
   StVector3(ccoord &x_,ccoord &y_,ccoord &z_)
   {x=x_;y=y_;z=z_;}

   void clear(){x=0;y=0;z=0;}

   void operator+=(const StVector3 &v);
   void operator*=(const coordinate &m);
   StVector3 operator*(const coordinate &m);

};

//---------------------------------------------------------------------------
struct StAxis{
union{
    struct{
    coordinate a,b,c,xo,yo,zo;
    };
    coordinate prm[6];
    };
bool on=false;
coordinate tmax,tmin;

//StVector3 begin,end;

    StAxis(const float &x, const float &y, const float &z){
        a=x;b=y;c=z;
    }

};
//---------------------------------------------------------------------------
///https://en.wikipedia.org/wiki/Rotation_matrix
///
class StRotationMatrix{
public:

union{
      struct{ coordinate m11,m12,m13,m21,m22,m23,m31,m32,m33;};
      coordinate m[9];
      coordinate mm[3][3];
};

coordinate ux,uy,uz;//axis of rotation
coordinate theta;
bool on=false;

    void buildMatrix(const StAxis &axis_);
    void buildMatrix(const StAxis &axis_,const coordinate &angle);
    StVector3 operator*(const StVector3 &a);
    //StAtom operator*=(StAtom &atom_);

private:
    void buildRotationAxis(const StAxis &axis_);
    void normUxyz();
};

//---------------------------------------------------------------------------
struct StGrainAtom{

wstring name;

union{
     StVector3 xyz;
     struct { coordinate x,y,z;};
};

    StGrainAtom(){ }
    StGrainAtom(const wstring &name_,position &x_, position &y_, position &z_){

        name=name_;
        x=x_;
        y=y_;
        z=z_;
        //r2=x*x+y*y+z*z;
    }
    void translate(cpos &x__, cpos &y__, cpos &z__)
    {
        x+=x__;
        y+=y__;
        z+=z__;
    }

friend wostream & operator << (wostream &o, StGrainAtom &atom);
};


typedef vector<StGrainAtom> vatoms;

//---------------------------------------------------------------------------

struct StGrainPrm{
float radiusI,radius,height,lattConst,T;
float radiusM,lattConstM,rel,std;
//position maxDist,maxX,maxY,maxZ;

union{
    struct {
    position minX,maxX;
    position minY,maxY;
    position minZ,maxZ;
    };
    position minmax[6];
};
enum Eminmax{MINX,MAXX,MINY,MAXY,MINZ,MAXZ};
enum Eshape{sphere=0,cyl,cube};
Eshape shape;

bool bT,bilatt,modified=false;
enum unitCell{SC,BCC,FCC,HCP,DIAM,XYZ,TEST};
unitCell ucType;

wstring fileName;
struct  StAtomNameNumber {
        wstring name; size_t count;
        StAtomNameNumber() {name=L""; count=0;}  ;
        StAtomNameNumber(const wstring name__) { name=std::move(name__); count=0;} ;


        } ;

vector<StAtomNameNumber> atomsNameNumber;
wstring atomMassA,atomMassB;


int atomLA,atomLB;
size_t numOfatoms,numOfAatoms;
unsigned number;

//bool remrand=false;



        StGrainPrm(){ }
        StGrainPrm(const double &radius,const double &lattConst){this->radius=radius;this->lattConst=lattConst;}
        StGrainPrm(const double &radius,const double &lattConst, const unitCell &uc)
                {this->radius=radius;this->lattConst=lattConst;ucType=uc;}

        bool renameAtomName(const wstring &from__, const wstring &to__)
        {
        bool retValue=false;

                for(auto &ann: atomsNameNumber){
                        if(ann.name==from__){
                            ann.name==to__;
                            retValue=true;
                        break;
                        }
                }


        return retValue;
        }



enum eScFactors{SF_OK,SF_FILE_ERROR,SF_NO_ITEM_EXISTS};


       static eScFactors getAtomMass(const string & aname,string &amass);

       friend ostream& operator << (ostream &os, const StGrainPrm &);
       friend wofstream& operator << (wofstream &os, const StGrainPrm &);
};
//------------------------------------------------------------------------------


struct StFileHeaderNDL{
StGrainPrm *ptrGrainPrm;
int dataRows,dataCols;
string crystStruct;
unsigned fpos_size,fwidth_size;
bool defaultSizeRC;
unsigned sizeRow,sizeCol;
//unsigned fpos_atoms,fwidth_atoms;

        StFileHeaderNDL() { }

        ~StFileHeaderNDL() { }

enum EnHeadErr{OK,WrongFormat};

//        void setSizeRC(const unsigned sr, const unsigned sc);
//        void updateSizeTag(wofstream &file,const unsigned &r, const unsigned &c);
//        void updateAtomsTag(wofstream &file,const QString &atoms);


//        StFileHeaderNDL::EnHeadErr removeExtHeader(wifstream &file);


//friend wofstream& operator << (wofstream &, StFileHeader &);
friend wifstream  & operator >> (wifstream &, StFileHeaderNDL &);
};



#endif // STNDLSTGRAIN

