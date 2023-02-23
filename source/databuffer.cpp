#include "databuffer.h"

#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>

#include <cstdlib>


float dataBuffer::xmargin=0;
float dataBuffer::ymargin=0;
float dataBuffer::zmargin=0;
bool  dataBuffer::setbox=false;
float dataBuffer::box[]={0,0,0,0,0,0};


bool dataBuffer::verbose=false;
vdb  *dataBuffer::ptr_vdb=nullptr;


using namespace std;

//cuns XRANGE=1;
//cuns YRANGE=1<<1;
//cuns ZRANGE=1<<2;


void dataBuffer::doCmdList(const bool &verb)
{
        verbose=verb;

        for(Stkeyvalue &kv: cmdList){
            if(verb)
                cout<<"::: "<<kv.key<<endl;


            //-------------------------------------------------------
            if(kv.key.find("translate")!=string::npos){
            const float x=std::stof(kv.valueA);
            const float y=std::stof(kv.valueB);
            const float z=std::stof(kv.valueC);

                translateAtoms(x,y,z);

                grainprm.maxX+=x;
                grainprm.minX+=x;

                grainprm.maxY+=y;
                grainprm.minY+=y;


                grainprm.maxZ+=z;
                grainprm.minZ+=z;

                continue;
            }
            //-------------------------------------------------------
            if(kv.key.find("rotate")!=string::npos){
            const float angle=std::stof(kv.valueA);
            const float x=std::stof(kv.valueB);
            const float y=std::stof(kv.valueC);
            const float z=std::stof(kv.valueD);

                rotateAtoms(angle,x,y,z);
                continue;
            }
            //-------------------------------------------------------
            if(kv.key.find("remove")!=string::npos){

                if(kv.valueA.find("cyl")!=string::npos){
                        if(kv.valueB.find("xy")!=string::npos){
                        bool remout=(kv.valueC.find("out")!=string::npos)?true: false;
                        const float radii=std::stof(kv.valueD);
                            removeAtomsCylXY(radii,remout);
                        }
                }

                if(kv.valueA.find("plane")!=string::npos){
                bool remout=(kv.valueB.find("out")!=string::npos)?true: false;

                      //  removeAtomsPlane(std::stof(kv.valueC),std::stof(kv.valueD),
                                       //  std::stof(kv.valueE),std::stof(kv.valueF),remout);

                }
                continue;
            }
            //-------------------------------------------------------
            if(kv.key.find("randpos")!=string::npos){

            cfloat dx=std::stof(kv.valueA);
            cfloat dy=std::stof(kv.valueB);
            cfloat dz=std::stof(kv.valueC);

                    randDislocation(dx,dy,dz);

            }
            //-------------------------------------------------------

        }

}

//----------------------------------------------------------------
void dataBuffer::translateAtoms(const float &x__, const float &y__, const float &z__)
{

        for(StGrainAtom &atom: atoms){
            atom.x+=x__;
            atom.y+=y__;
            atom.z+=z__;
        }

}

//----------------------------------------------------------------
void dataBuffer::rotateAtoms(const float &angle__, const float &x__, const float &y__, const float &z__)
{
StAxis axis(x__,y__,z__);
StRotationMatrix rotMat;

        rotMat.buildMatrix(axis,angle__*M_PI/180);

        for(StGrainAtom &atom:atoms)
            atom.xyz=rotMat*atom.xyz;

}
//----------------------------------------------------------------
void dataBuffer::renameAtoms(const wstring &from, const wstring &to)
{

        if( grainprm.renameAtomName(from,to)){

                for(StGrainAtom &atom : atoms)
                        if(atom.name==from)
                            atom.name=to;
        }
        else
            wcerr<<"WARNING: renaming failure, unknown name: "<<from<<endl;
}

//----------------------------------------------------------------
void dataBuffer::removeAtomsCylXY(const float &rxy__, const bool & out, const float &h)
{
float xx,yy,hh;
float minX,maxX;
float minY,maxY;
float minZ,maxZ;
const float rxy2=rxy__*rxy__;
size_t  &aremoved=removedAtoms;

vector<StGrainAtom> tmpAtoms;

        tmpAtoms.reserve(atoms.size());

        minX=minY=minZ=std::numeric_limits<float>::max();
        maxX=maxY=maxZ=-minX;


        for(StGrainAtom &atom: atoms){
            xx=atom.x*atom.x;
            yy=atom.y*atom.y;
            hh=std::fabs(atom.z*0.5);

            if( (xx+yy<=rxy2 && hh<=h ) ^ !out ){
                tmpAtoms.push_back(atom);

                if(atom.x<minX) minX=atom.x;
                if(atom.x>maxX) maxX=atom.x;

                if(atom.y<minY) minY=atom.y;
                if(atom.y>maxY) maxY=atom.y;

                if(atom.z<minZ) minZ=atom.z;
                if(atom.z>maxZ) maxZ=atom.z;

            }
            else
                aremoved++;
       }

        //atoms.assign(tmpAtoms.begin(),tmpAtoms.end());
        atoms=std::move(tmpAtoms);
        grainprm.numOfatoms=atoms.size();
        grainprm.maxX=maxX;
        grainprm.minX=minX;
        grainprm.maxY=maxY;
        grainprm.minY=minY;
        grainprm.maxZ=maxZ;
        grainprm.minZ=minZ;

        if(verbose){
            cout<<":::    atom(s) removed "<<aremoved<<endl;
            cout<<":::    atom(s) in buffer "<<grainprm.numOfatoms<<endl;
        }
}
//----------------------------------------------------------------
void dataBuffer::removeAtomsSph(const float &r__, const bool &out)
{
float xx,yy,zz;
float minX,maxX;
float minY,maxY;
float minZ,maxZ;
const float r2=r__*r__;
size_t  &aremoved=removedAtoms;

vector<StGrainAtom> tmpAtoms;

        tmpAtoms.reserve(atoms.size());

        minX=minY=minZ=std::numeric_limits<float>::max();
        maxX=maxY=maxZ=-minX;


        for(StGrainAtom &atom: atoms){
            xx=atom.x*atom.x;
            yy=atom.y*atom.y;
            zz=atom.z*atom.z;

            if( (xx+yy+zz<=r2) ^ !out ){
                tmpAtoms.push_back(atom);

                if(atom.x<minX) minX=atom.x;
                if(atom.x>maxX) maxX=atom.x;

                if(atom.y<minY) minY=atom.y;
                if(atom.y>maxY) maxY=atom.y;

                if(atom.z<minZ) minZ=atom.z;
                if(atom.z>maxZ) maxZ=atom.z;

            }
            else
                aremoved++;
       }

        atoms.assign(tmpAtoms.begin(),tmpAtoms.end());
        grainprm.numOfatoms=atoms.size();
        grainprm.maxX=maxX;
        grainprm.minX=minX;
        grainprm.maxY=maxY;
        grainprm.minY=minY;
        grainprm.maxZ=maxZ;
        grainprm.minZ=minZ;

        if(verbose){
            cout<<":::    atom(s) removed "<<aremoved<<endl;
            cout<<":::    atom(s) in buffer "<<grainprm.numOfatoms<<endl;
        }
}
//----------------------------------------------------------------
void dataBuffer::removeAtomsPlane(const float &a, const float &b, const float &c, const float &d, const bool &out,const string &fileName)
{
size_t  &aremoved=removedAtoms;
vector<StGrainAtom> tmpAtoms;
float v;
float minX,maxX;
float minY,maxY;
float minZ,maxZ;
wfstream fileForRmAtoms;

            tmpAtoms.reserve(atoms.size());

            minX=minY=minZ=std::numeric_limits<float>::max();
            maxX=maxY=maxZ=-minX;


            if(!fileName.empty()){
                fileForRmAtoms.open(fileName,ios::out);
                if(!fileForRmAtoms){
                    cerr<<"ERROR: couldn't open a file "<<fileName<<endl;
                    cerr<<"       removed atoms will not be saved"<<endl;
                }

            }


            if(fileForRmAtoms.is_open()){
            size_t aremovedLoc=0;
                fileForRmAtoms<<"          "<<endl;
                fileForRmAtoms<<"removed atoms"<<endl;

                for(StGrainAtom &atom:atoms){
                    v=a*atom.x+b*atom.y+c*atom.z+d;

                    if( (v<=0) ^ !out){
                        tmpAtoms.push_back(atom);

                        if(atom.x<minX) minX=atom.x;
                        if(atom.x>maxX) maxX=atom.x;

                        if(atom.y<minY) minY=atom.y;
                        if(atom.y>maxY) maxY=atom.y;

                        if(atom.z<minZ) minZ=atom.z;
                        if(atom.z>maxZ) maxZ=atom.z;
                    }
                    else{
                        aremoved++;
                        aremovedLoc++;
                        fileForRmAtoms<<atom<<endl;
                    }
                }

                fileForRmAtoms.seekp(0,ios::beg);
                fileForRmAtoms<<aremovedLoc;
                fileForRmAtoms.close();

                if(!aremovedLoc)
                    cout<<"WARNING: file is empty"<<endl;
            }
            else{

                for(StGrainAtom &atom:atoms){
                    v=a*atom.x+b*atom.y+c*atom.z+d;

                    if( (v<=0) ^ !out){
                        tmpAtoms.push_back(atom);

                        if(atom.x<minX) minX=atom.x;
                        if(atom.x>maxX) maxX=atom.x;

                        if(atom.y<minY) minY=atom.y;
                        if(atom.y>maxY) maxY=atom.y;

                        if(atom.z<minZ) minZ=atom.z;
                        if(atom.z>maxZ) maxZ=atom.z;

                    }
                    else
                        aremoved++;
                }
            }

            atoms.assign(tmpAtoms.begin(),tmpAtoms.end());
            grainprm.numOfatoms=atoms.size();

            grainprm.maxX=maxX;
            grainprm.minX=minX;

            grainprm.maxY=maxY;
            grainprm.minY=minY;

            grainprm.maxZ=maxZ;
            grainprm.minZ=minZ;

            if(verbose){
                cout<<__LINE__<<" :::    atom(s) removed "<<aremoved<<endl;
                cout<<__LINE__<<" :::    atom(s) in buffer "<<grainprm.numOfatoms<<endl;
            }
}
//----------------------------------------------------------------
void dataBuffer::removeAtomsPlaneRand(const float &a, const float &b, const float &c, const float &d, const float &e, const bool &out)
{
default_random_engine generator;
uniform_real_distribution<float> distribution(d,e);
const string emptyFileName;


            generator.seed(chrono::system_clock::now().time_since_epoch().count());
            removeAtomsPlane(a,b,c,distribution(generator),out,emptyFileName);
}
//----------------------------------------------------------------
void dataBuffer::removeAtomsSlice(const float &a, const float &b, const float &c, const float &d, const float &e, const bool &out)
{
    size_t  &aremoved=removedAtoms;
    vector<StGrainAtom> tmpAtoms;
    float v0,v1;
    float minX,maxX;
    float minY,maxY;
    float minZ,maxZ;

                tmpAtoms.reserve(atoms.size());


                minX=minY=minZ=std::numeric_limits<float>::max();
                maxX=maxY=maxZ=-minX;


                for(StGrainAtom &atom:atoms){
                    v0=a*atom.x+b*atom.y+c*atom.z+d;
                    v1=v0-e;

                    if(  ( (v1<=0) ^ out) || ( (v0<=0) ^ !out)  ) {
                        tmpAtoms.push_back(atom);

                        if(atom.x<minX) minX=atom.x;
                        if(atom.x>maxX) maxX=atom.x;

                        if(atom.y<minY) minY=atom.y;
                        if(atom.y>maxY) maxY=atom.y;

                        if(atom.z<minZ) minZ=atom.z;
                        if(atom.z>maxZ) maxZ=atom.z;

                    }
                    else
                        aremoved++;
                }

                atoms.assign(tmpAtoms.begin(),tmpAtoms.end());
                grainprm.numOfatoms=atoms.size();

                grainprm.maxX=maxX;
                grainprm.minX=minX;
                grainprm.maxY=maxY;
                grainprm.minY=minY;
                grainprm.maxZ=maxZ;
                grainprm.minZ=minZ;

                if(verbose){
                    cout<<":::    atom(s) removed "<<aremoved<<endl;
                    cout<<":::    atom(s) in buffer "<<grainprm.numOfatoms<<endl;
                }

}
//----------------------------------------------------------------
bool fMin(const position &value,const position & xyzMinMax)
{
return value < xyzMinMax;
}
//----------------------------------------------------------------
bool fMax(const position &value,const position & xyzMinMax)
{
return value > xyzMinMax;
}

//----------------------------------------------------------------
void dataBuffer::removeAtomsSlice(const char &xyz, const bool &minmax, const float &thickness)
{
const size_t cmp=xyz-'x';
const size_t cMinMax=2*cmp+(size_t)minmax;
vector<StGrainAtom> tmpAtoms;
float xyzMinMax=grainprm.minmax[cMinMax]+thickness;   //
auto fcmp=(minmax) ? &fMax : &fMin ;

float min,max;


            tmpAtoms.reserve(atoms.size());

            min=INFINITY;
            max=-min;

            for(StGrainAtom &atom:atoms){

                if( ! fcmp(atom.xyz.xyz[cmp], xyzMinMax) ){
                    tmpAtoms.emplace_back(atom);


                    if(atom.xyz.xyz[cmp]>max) max=atom.xyz.xyz[cmp];
                    if(atom.xyz.xyz[cmp]<min) min=atom.xyz.xyz[cmp];

                }
                else
                    removedAtoms++;

            }

            atoms.assign(tmpAtoms.begin(),tmpAtoms.end());

            grainprm.numOfatoms=atoms.size();
            grainprm.minmax[2*cmp]=min;
            grainprm.minmax[2*cmp+1]=max;

}
//----------------------------------------------------------------
void dataBuffer::removeAtomsCuboid(const float &a, const float &b, const bool &out)
{
    size_t  &aremoved=removedAtoms;
    vector<StGrainAtom> tmpAtoms;

    float minX,maxX;
    float minY,maxY;
    float minZ,maxZ;

    const float ap2= a*0.5;
    const float am2=-a*0.5;

    const float bp2= b*0.5;
    const float bm2=-b*0.5;

                tmpAtoms.reserve(atoms.size());


                minX=minY=minZ=std::numeric_limits<float>::max();
                maxX=maxY=maxZ=-minX;


                for(StGrainAtom &atom:atoms){

                    if ( (atom.x<=ap2 && atom.x>=am2 && atom.y<=bp2 && atom.y>=bm2) ^ !out ) {
                       tmpAtoms.push_back(atom);


                       if(atom.x<minX) minX=atom.x;
                       if(atom.x>maxX) maxX=atom.x;

                       if(atom.y<minY) minY=atom.y;
                       if(atom.y>maxY) maxY=atom.y;

                       if(atom.z<minZ) minZ=atom.z;
                       if(atom.z>maxZ) maxZ=atom.z;


                    }
                    else
                        aremoved++;
                }


                atoms=std::move(tmpAtoms);
                grainprm.numOfatoms=atoms.size();


                grainprm.maxX=maxX;
                grainprm.minX=minX;
                grainprm.maxY=maxY;
                grainprm.minY=minY;
                grainprm.maxZ=maxZ;
                grainprm.minZ=minZ;


}
//----------------------------------------------------------------
void dataBuffer::removeAtomsRandom(const float &a)
{
size_t  &aremoved=removedAtoms;
vector<StGrainAtom> tmpAtoms;

default_random_engine generator;
uniform_real_distribution<double> distribution(0,100);

size_t removedRandAtoms=0;
const double k0=100.0/atoms.size();

float minX,maxX;
float minY,maxY;
float minZ,maxZ;
double numrandom;



            generator.seed(chrono::system_clock::now().time_since_epoch().count());

            tmpAtoms.reserve(atoms.size());

            minX=minY=minZ=std::numeric_limits<float>::max();
            maxX=maxY=maxZ=-minX;


             for(StGrainAtom &atom:atoms){
                 numrandom=distribution(generator);

                 if(numrandom>a){
                     tmpAtoms.emplace_back(atom);


                     if(atom.x<minX) minX=atom.x;
                     if(atom.x>maxX) maxX=atom.x;

                     if(atom.y<minY) minY=atom.y;
                     if(atom.y>maxY) maxY=atom.y;

                     if(atom.z<minZ) minZ=atom.z;
                     if(atom.z>maxZ) maxZ=atom.z;
                 }
                 else{
                     aremoved++;
                     removedRandAtoms++;
                 }

             }


             atoms=std::move(tmpAtoms);
             grainprm.numOfatoms=atoms.size();


             grainprm.maxX=maxX;
             grainprm.minX=minX;
             grainprm.maxY=maxY;
             grainprm.minY=minY;
             grainprm.maxZ=maxZ;
             grainprm.minZ=minZ;


             if(verbose){
                 cout<<"::: random atoms removed "<<removedRandAtoms<<",  "<<removedRandAtoms*k0<<"%"<<endl;
                 cout<<":::    atom(s) in buffer "<<grainprm.numOfatoms<<endl;
             }


}
//----------------------------------------------------------------
void dataBuffer::randDislocation(cfloat &dx__, cfloat &dy__, cfloat &dz__)
{
std::srand ( unsigned ( std::time(0) ) );
float rdx,rdy,rdz;
cfloat invRM=1.0f/RAND_MAX;

            for(StGrainAtom &atom:atoms){
                rdx=std::rand()*invRM-0.5;
                rdy=std::rand()*invRM-0.5;
                rdz=std::rand()*invRM-0.5;

                rdx*=dx__;
                rdy*=dy__;
                rdz*=dz__;

                atom.x+=rdx;
                atom.y+=rdy;
                atom.z+=rdz;
            }


}
//----------------------------------------------------------------
void dataBuffer::center()
{
float cx,cy,cz;
const size_t numOfatoms=atoms.size();

            if(!numOfatoms){
                cerr<<"::: WARNING: centering is impossible, number of atoms is equal 0\n";
            return;
            }

const float inumOfatoms=-1.0f/numOfatoms;


            cx=cy=cz=0;

            for(StGrainAtom &atom:atoms){
                cx+=atom.x;
                cy+=atom.y;
                cz+=atom.z;
            }

            cx*=inumOfatoms;
            cy*=inumOfatoms;
            cz*=inumOfatoms;

            if(verbose){
                cout<<"::: translation vector of centering: "<<cx<<", "<<cy<<", "<<cz<<endl;
            }

            for(StGrainAtom &atom:atoms){
                atom.x+=cx;
                atom.y+=cy;
                atom.z+=cz;
            }
}

//----------------------------------------------------------------
float & aretX( StGrainAtom &atom) { return atom.x; }
float & aretY( StGrainAtom &atom) { return atom.y; }
float & aretZ( StGrainAtom &atom) { return atom.z; }
//----------------------------------------------------------------

void dataBuffer::markerLayer(const string &apos, const string &tol, const string &mode, const string &sig)
{
float & (*fptr)(StGrainAtom &)=(mode=="x") ? &aretX :(mode=="y") ? &aretY : &aretZ;
float pos;
const float aposH=std::stof(apos)+std::stof(tol);
const float aposL=std::stof(apos)-std::stof(tol);
wstring wsig;
size_t mcatoms=0;

            wsig.assign(sig.begin(),sig.end());

            for(StGrainAtom &atom:atoms){

                pos=fptr(atom);

                if( pos>= aposL && pos<=aposH){
                    atom.name+=wsig;
                    mcatoms++;
                }
            }


            if(verbose)
                cout<<"::: marked atoms "<<mcatoms<<endl;
}
//----------------------------------------------------------------

bool sortByAtomNameA (const StGrainAtom & lg, const StGrainAtom & rg)
{	return lg.name < rg.name; }

bool sortByAtomNameD (const StGrainAtom & lg, const StGrainAtom & rg)
{	return lg.name > rg.name; }




void dataBuffer::sortNames(bool asc)
{

//bool (*fsort)(const StGrainAtom &, const StGrainAtom &)=(asc)? &sortByAtomNameA: &sortByAtomNameD;

//        sort(atoms.begin(),atoms.end(),fsort);

}
//----------------------------------------------------------------
void dataBuffer::intspace(const string &xyz, const string &position, const string &space)
{
float & (*fptr)( StGrainAtom &)=(xyz=="x") ? &aretX :(xyz=="y") ? &aretY : &aretZ;
const float fposition=std::stof(position);
const float fspace=std::stof(space);

                for(StGrainAtom &atom:atoms)
                    fptr(atom)+=(fptr(atom)>fposition)? fspace : -fspace;

}
//----------------------------------------------------------------
void dataBuffer::printStatistic(const string &option__)
{

            if(option__=="xrange")
                cout<<"xmin, xmax: "<<grainprm.minX<<"\t"<<grainprm.maxX<<"\n";
            else
                if(option__=="yrange")
                    cout<<"ymin, ymax: "<<grainprm.minY<<"\t"<<grainprm.maxY<<"\n";
                else
                    if(option__=="zrange")
                        cout<<"zmin, zmax: "<<grainprm.minZ<<"\t"<<grainprm.maxZ<<"\n";                                                        
                    else
                        if(option__=="name")
                            cout<<"dataName "<<dataName<<"\n";
                        else
                            if(option__=="rematoms")
                                cout<<"removed "<<removedAtoms<<"\n";
                            else
                                if(option__=="file")
                                    cout<<"file "<<fileName<<"\n";
                                else
                                    if(option__=="numatoms"){
                                        cout<<"num of atoms: ";
                                        if(grainprm.bilatt)
                                            cout<<"A "<<grainprm.numOfAatoms<<" B "<<grainprm.numOfatoms-grainprm.numOfAatoms<<"\n";
                                        else
                                            cout<<"A "<<grainprm.numOfatoms<<"\n";
                                    }
                                    else{
                                        //separator
                                        if(option__.length()==9)
                                            cout<<"\n";
                                        else{
                                        char sepchar=option__[9];
                                        int len=1;
                                            if(option__.length()>10){
                                            string slen=option__.substr(10) ;
                                                    len=std::stoi(slen);
                                            }

                                            std::cout<<setfill(sepchar)<<setw(len)<<sepchar<<"\n";
                                            std::cout<<setfill(' ');
                                        }
                                    }            
}
//----------------------------------------------------------------
void dataBuffer::rescalePosition(const float &x, const float &y, const float &z)
{

        for(StGrainAtom &atom: atoms){
                    atom.x*=x;
                    atom.y*=y;
                    atom.z*=z;
        }
}
//----------------------------------------------------------------
void dataBuffer::slices(const string &dirName, const float &start, const float &step, const float &stop)
{

const int N=(stop-start)/step;

            if(N<=0){
                cerr<<"ERROR: wrong numerical values of 'slices' command"<<endl;
            return;
            }


const int errDir = system(std::string("mkdir -p "+dirName).c_str());

            if (errDir==-1)
            {
                cerr<<"ERROR: coulnd't create a directory"<<endl;
            return;
            }


            // remove content if any
            system(std::string("rm -f  ./"+dirName+"/*").c_str());

//StGrainAtom &atom:atoms
//vatoms sliceAtoms;

const float a=(N-1.0)/(stop-start);
const float b=-start*a;

vector<vatoms> vslices;
int sliceNum;
size_t appAtomsPerSlice=1+atoms.size()/N;


            vslices.resize(N);

            for(auto &slice: vslices)
                slice.reserve(2*appAtomsPerSlice);


            for(const StGrainAtom &atom : atoms ){
                sliceNum=(int)(a*atom.z+b);

                if(sliceNum>=0 && sliceNum<N){
                    vslices[sliceNum].push_back(atom);
                }
            }


            for(int i=0;i<N;i++){
            string fileName(dirName+"/"+std::to_string(i)+".xyz");
            wofstream fout(fileName,ios::out);
                    if(!fout){
                        cerr<<"ERROR: couldn't save  slice "<<i<<endl;
                    break;
                    }

            vatoms &atoms=vslices[i];
            const size_t sN=atoms.size();

                    fout<<sN<<endl;
                    fout<<"slice "<<i<<" "<<(i-b)/a<<endl;

                    for(size_t j=0;j<sN;j++)
                        fout<<atoms[j]<<endl;

                    fout.close();

            }


}

