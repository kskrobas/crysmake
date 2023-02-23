#ifndef DATABUFFER_H
#define DATABUFFER_H


#include <string>
#include <fstream>
#include <vector>
#include "stndlstgrain.h"


#define SET_BIT(val, bitIndex) val |= (1 << bitIndex)
#define CLEAR_BIT(val, bitIndex) val &= ~(1 << bitIndex)
#define TOGGLE_BIT(val, bitIndex) val ^= (1 << bitIndex)
#define BIT_IS_SET(val, bitIndex) (val & (1 << bitIndex))


#define XRANGE 1
#define YRANGE 2
#define ZRANGE 3
#define NUMATOMS 4

//typedef const unsigned cuns;
//extern cuns XRANGE;


struct dataBuffer;
typedef std::vector<dataBuffer> vdb;
typedef const float cfloat;
//-----------------------------------------------------------------------------

struct Stkeyvalue{
std::string key;
std::string valueA,valueB,valueC,valueD,valueE,valueF;

    Stkeyvalue(){ };

    Stkeyvalue(const std::string key__,const std::string &value__)
    {
        key=key__;
        valueA=value__;
    }

    Stkeyvalue(const std::string key__,
               const std::string &valueA__,const std::string &valueB__,
               const std::string &valueC__)
    {
        key=key__;
        valueA=valueA__;
        valueB=valueB__;
        valueC=valueC__;
    }

    Stkeyvalue(const std::string key__,
               const std::string &valueA__,const std::string &valueB__,
               const std::string &valueC__,const std::string &valueD__)
    {
        key=key__;
        valueA=valueA__;
        valueB=valueB__;
        valueC=valueC__;
        valueD=valueD__;
    }

    Stkeyvalue(const std::string key__,
               const std::string &valueA__,const std::string &valueB__,
               const std::string &valueC__,const std::string &valueD__,
               const std::string &valueE__,const std::string &valueF__)
    {

        //Stkeyvalue(key,valueA__,valueB__,valueC__,valueD__);
        key=key__;
        valueA=valueA__;
        valueB=valueB__;
        valueC=valueC__;
        valueD=valueD__;
        valueE=valueE__;
        valueF=valueF__;
    }
};

///----------------------------------------------


struct StRandomizePrm{
float dx,dy,dz;
};

//-----------------------------------------------------------------------------
struct dataBuffer{
StGrainPrm grainprm;
vatoms atoms;

static bool setbox;
static float xmargin,ymargin,zmargin;
enum EBox{xmin,xmax,ymin,ymax,zmin,zmax};
static float box[6];

static vdb *ptr_vdb;
float xt,yt,zt;
size_t removedAtoms;
std::string fileName,dataName;
unsigned char dispOptions=0;
bool sort;




        dataBuffer(){
            atoms.clear();
            fileName="?";
            dataName="?";
            removedAtoms=0;
            sort=true;
            //xshift=0;yshift=0;zshift=0;
        }
        dataBuffer(std::string &file__){
            fileName=file__;
            removedAtoms=0;

        }

vector<Stkeyvalue> cmdList;
static bool verbose;
        //----------------------------------------------------------------
        void doCmdList(const bool &verb=false);

        //----------------------------------------------------------------
        void translateAtoms(const float &x__, const float &y__, const float &z__);
        void rotateAtoms(const float &angle__,const float &x__, const float &y__, const float &z__);
        void renameAtoms(const wstring &from, const wstring &to);
        void removeAtomsCylXY(const float &rxy__,const bool &out=true,const float &h=1);
        void removeAtomsSph(const float &rxy__,const bool &out=true);
        void removeAtomsPlane(const float &a, const float &b, const float &c, const float &d,
                              const bool &out, const string &fileName);
        void removeAtomsPlaneRand(const float &a, const float &b, const float &c, const float &d,const float &e,
                                  const bool &out=true);
        void removeAtomsSlice(const float &a, const float &b, const float &c, const float &d,const float &e,
                              const bool &out=true);
        void removeAtomsSlice(const char &xyz,const bool &minmax,const float &thickness);
        void removeAtomsCuboid(const float &a, const float &b, const bool &out=true);
        void removeAtomsRandom(const float &a);

        void randDislocation(cfloat &dx__, cfloat &dy__,cfloat &dz__);
        void center();
        void markerLayer(const string &apos, const string &tol, const string &mode,
                         const string & sig);

        void sortNames(bool asc);


        //intspace {x|y|z} {pos} {space}
        void intspace(const string &xyz, const string &pos, const string &space);

        //void pushAtoms(const string &dataName__);
        void printStatistic(const string &option__);

        void rescalePosition(const float &x,const float &y,const float &z);

        void slices(const string &dirName,const float &start, const float & step, const float & stop);

        //----------------------------------------------------------------


};


//-----------------------------------------------------------------------------

#endif // DATABUFFER_H

