#include "iofile.h"
#include "elements.h"
#include <limits>
#include <map>

#include<string>
#include <codecvt>
#include <locale>
#include <algorithm>
#include <chrono>


// conversion
using convert_t = std::codecvt_utf8<wchar_t>;
std::wstring_convert<convert_t, wchar_t> strconverter;

std::string to_string(std::wstring wstr)
{
    return strconverter.to_bytes(wstr);
}

std::wstring to_wstring(std::string str)
{
    return strconverter.from_bytes(str);
}

using namespace std;
//-----------------------------------------------------------------------------
bool sortByAtomName (const StGrainAtom & lg, const StGrainAtom & rg)
{	return lg.name < rg.name; }




//-----------------------------------------------------------------------------
bool loadAtomsFromNDL(dataBuffer &db,bool verbose)
{
const string fileName=db.fileName;
wifstream file(fileName.c_str(),ios::in);

                if(!file){
                    cerr<<endl<<"ERROR: couldn't find the input file "<<endl;
                    file.close();
                return false;
                }
                else
                    if(verbose) cout<<"::: file reading: "<<fileName<<endl;

StFileHeaderNDL fndl;
                fndl.ptrGrainPrm=&db.grainprm;

                file.exceptions(ios::badbit );
                try{
                    file>>fndl;
                }
                catch(...){
                    cerr<<endl<<" header error "<<endl;
                    file.close();
                return false;
                }

wstring atname;
position x,y,z;
position minZ,maxZ,minX,maxX,minY,maxY;
//position x2max,y2max,z2max;
size_t i;

                db.atoms.reserve(db.grainprm.numOfatoms);

                try{

                        //x2max=y2max=z2max=0;

                        file>>atname>>x>>y>>z;
                        db.atoms.push_back(StGrainAtom(atname,x,y,z));
                        minX=x;
                        maxX=x;

                        minY=y;
                        maxY=y;

                        minZ=z;
                        maxZ=z;

                        for(i=1;i<db.grainprm.numOfatoms;i++){
                            file>>atname>>x>>y>>z;

                            db.atoms.emplace_back(StGrainAtom(atname,x,y,z));

                            if(x>maxX) maxX=x;
                            if(x<minX) minX=x;

                            if(y>maxY) maxY=y;
                            if(y<minY) minY=y;

                            if(z>maxZ) maxZ=z;
                            if(z<minZ) minZ=z;
                        }


                        db.grainprm.minX=minX;
                        db.grainprm.maxX=maxX;

                        db.grainprm.minY=minY;
                        db.grainprm.maxY=maxY;

                        db.grainprm.minZ=minZ;
                        db.grainprm.maxZ=maxZ;

                        if(verbose){

                            cout<<"::: "<<setw(2)<<std::right<<"*"<<setw(10)<<std::right<<"min"<<setw(10)<<std::right<<"max"<<endl;
                            cout<<"::: "<<setw(2)<<std::right<<"X"<<setw(10)<<std::right<<minX<<setw(10)<<std::right<<maxX<<endl;
                            cout<<"::: "<<setw(2)<<std::right<<"Y"<<setw(10)<<std::right<<minY<<setw(10)<<std::right<<maxY<<endl;
                            cout<<"::: "<<setw(2)<<std::right<<"Z"<<setw(10)<<std::right<<minZ<<setw(10)<<std::right<<maxZ<<endl;
                        }


                }
                catch(...){
                    cerr<<endl<<" error during data reading, row:  "<<i<<endl;
                    file.close();
                    return false;
                }

            file.close();

            db.atoms.shrink_to_fit();

            if(db.grainprm.bilatt){


                //sort (db.atoms.begin(),db.atoms.end(),sortByAtomName);

                for(db.grainprm.numOfAatoms=1;
                    db.grainprm.numOfAatoms<db.grainprm.numOfatoms;
                    db.grainprm.numOfAatoms++)
                {
                        if(db.atoms[0].name!=db.atoms[db.grainprm.numOfAatoms].name){
                            break;
                         }
                }


                    if(verbose){
                        wcout<<"\n\n::: atom types listing\n::: type number\n";

                        for(auto &ann: db.grainprm.atomsNameNumber)
                            wcout<<"   "<<ann.name<<"\t"<<ann.count<<endl;


                        wcout<<"::: reading SUCCESS\n"<<endl;
                    }


            }
            else
                db.grainprm.numOfAatoms=db.grainprm.numOfatoms;

return true;
}
/////---------------------------------------------------------------------------
//bool saveDataToConfig(std::string &filename,dataBuffer *db, const size_t dbsize)
//{

//}
///---------------------------------------------------------------------------
bool saveDataToXYZ(const std::string &filename, vdb *ptr_vdb)
{
wfstream file(filename,ios::out);

size_t numOfatoms=0;
const size_t dbsize=ptr_vdb->size();


                for(dataBuffer &db: *ptr_vdb)
                        numOfatoms+=db.atoms.size();


                if(!numOfatoms){

                    //if(verbose){
                        cout<<"::: WARNING, no atoms to be saved\n";
                    // }

                        file<<numOfatoms<<endl;
                        file<<"empty set"<<endl;


                    return true;
                }


                file<<numOfatoms<<endl;
                file<<"num of molecules: "<<dbsize<<endl;
				
				constexpr double mnoz=10;
				
				//cerr<<"  UWAGA: pozycje przemnozona przez "<<mnoz<<endl;


                for(dataBuffer &db: *ptr_vdb){


                    for(StGrainAtom &atom: db.atoms)
                        file<<atom.name<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;
                }


                file.close();
return true;
}
///---------------------------------------------------------------------------
bool saveDataToCONFIG(const std::string &filename, vdb *ptr_vdb, bool verbose)
{
wfstream file(filename,ios::out);
size_t numOfatoms=0;
float xmin,xmax,ymin,ymax,zmin,zmax;


                xmin=std::numeric_limits<float>::max();//(*ptr_vdb)[0].atoms[0].x;
                xmax=-xmin;
                ymin= xmin;
                ymax=-xmin;
                zmin= xmin;
                zmax=-xmin;


                for(dataBuffer &db: *ptr_vdb){

                        if(db.atoms.empty()) {

                            if(verbose){
                                cout<<"::: WARNING, current buffer is empty \n";
                                cout<<":::   fileName: "<<db.fileName<<"\n";
                                cout<<":::   dataName: "<<db.dataName<<"\n";
                             }

                            continue;
                        }

                        numOfatoms+=db.atoms.size();


                       for(StGrainAtom &atom: db.atoms){

                            if(atom.x<xmin) xmin=atom.x;
                            if(atom.x>xmax) xmax=atom.x;

                            if(atom.y<ymin) ymin=atom.y;
                            if(atom.y>ymax) ymax=atom.y;

                            if(atom.z<zmin) zmin=atom.z;
                            if(atom.z>zmax) zmax=atom.z;
                       }


                 }


                 if(!numOfatoms){

                     if(verbose){
                         cout<<"::: WARNING, no atoms to be saved\n";
                      }

                     file.close();
                     return true;
                 }



const float w=(xmax-xmin)+2*dataBuffer::xmargin;
const float l=(ymax-ymin)+2*dataBuffer::ymargin;
const float h=(zmax-zmin)+2*dataBuffer::zmargin;

                if(verbose){
                    cout<<"::: "<<setw(2)<<std::right<<"*"<<setw(10)<<std::right<<"min"<<setw(10)<<std::right<<"max"<<endl;

                    cout<<"::: "<<setw(2)<<std::right<<"X"<<setw(10)<<std::right<<xmin<<setw(10)<<std::right<<xmax<<endl;
                    cout<<"::: "<<setw(2)<<std::right<<"Y"<<setw(10)<<std::right<<ymin<<setw(10)<<std::right<<ymax<<endl;
                    cout<<"::: "<<setw(2)<<std::right<<"Z"<<setw(10)<<std::right<<zmin<<setw(10)<<std::right<<zmax<<endl;
                }

                /// 1 row
                file<<"*** merged slices , total number of atoms "<<numOfatoms<<" ***"<<endl;
                /// 2 row
                file<<setw(10)<<std::right<<0<<setw(10)<<std::right<<2<<setw(10)<<std::right<<numOfatoms<<setw(20)<<std::right<<1<<endl;
                /// 3, 4, 5 row
                file<<setw(20)<<std::right<< w <<setw(20)<<std::right<<0<<setw(20)<<std::right<<0<<endl;
                file<<setw(20)<<std::right<<0<<setw(20)<<std::right<<  l <<setw(20)<<std::right<<0<<endl;
                file<<setw(20)<<std::right<<0<<setw(20)<<std::right<<0<<setw(20)<<std::right<<  h  <<endl;                       \


                /// atom positions
                size_t k=0;
                    for(dataBuffer &db : *ptr_vdb){
                        for(StGrainAtom &atom: db.atoms){
                            file<<atom.name<<"        "<<k++<<endl;
                            file<<setw(16)<<right<<atom.x<<"    ";
                            file<<setw(16)<<right<<atom.y<<"    ";
                            file<<setw(16)<<right<<atom.z<<"    ";
                            file<<endl;
                        }
                    }

                file.close();
return true;
}
//-----------------------------------------------------------------------------

bool saveDataToCOORD(const std::string &filename, vdb *ptr_vdb,bool verbose)
{
wfstream file(filename,ios::out);
size_t numOfatoms=0;
float xmin,xmax,ymin,ymax,zmin,zmax;


                xmin=std::numeric_limits<float>::max();//(*ptr_vdb)[0].atoms[0].x;
                xmax=-xmin;
                ymin= xmin;
                ymax=-xmin;
                zmin= xmin;
                zmax=-xmin;


                for(dataBuffer &db: *ptr_vdb){

                        if(db.atoms.empty()) {

                            if(verbose){
                                cout<<"::: WARNING, current buffer is empty \n";
                                cout<<":::   fileName: "<<db.fileName<<"\n";
                                cout<<":::   dataName: "<<db.dataName<<"\n";
                             }

                            continue;
                        }

                        numOfatoms+=db.atoms.size();


                       for(StGrainAtom &atom: db.atoms){

                            if(atom.x<xmin) xmin=atom.x;
                            if(atom.x>xmax) xmax=atom.x;

                            if(atom.y<ymin) ymin=atom.y;
                            if(atom.y>ymax) ymax=atom.y;

                            if(atom.z<zmin) zmin=atom.z;
                            if(atom.z>zmax) zmax=atom.z;
                       }


                        if(db.grainprm.bilatt){

                            sort(db.atoms.begin(),db.atoms.end(),sortByAtomName);
                        const wstring lastAtomName=db.atoms.back().name;

                            db.grainprm.numOfAatoms=0;

                            for(StGrainAtom &atom: db.atoms){
                                if(atom.name!=lastAtomName){
                                    db.grainprm.numOfAatoms++;
                                }
                                else
                                    break;
                            }

                            if(verbose){
                                wcout<<"::: num of "<<db.atoms[0].name<<": "<<db.grainprm.numOfAatoms<<endl;
                                wcout<<"::: num of "<<lastAtomName<<": "<<db.atoms.size()-db.grainprm.numOfAatoms<<endl;
                            }
                        }
                 }


                 if(!numOfatoms){

                     if(verbose){
                         cout<<"::: WARNING, no atoms to be saved\n";
                      }

                     file.close();
                     return true;
                 }

//const float h=(zmax-zmin)+dataBuffer::zmargin;
return true;
}
//-----------------------------------------------------------------------------
bool saveDataToLammps(const std::string &filename, vdb *ptr_vdb,bool verbose)
{
////wfstream file(filename,ios::out);
wfstream file(filename,ios::out);

            if(verbose)
                cout<<" reading data file"<<endl;



            if(!file){
                cerr<<endl<<"ERROR: couldn't open file for saving"<<endl;
                file.close();
            return false;
            }


size_t numOfatoms=0,i;
float xmin,xmax,ymin,ymax,zmin,zmax;
float *minmax[6];
std::map<std::wstring, std::wstring> elementMass;
std::map<std::wstring, std::wstring>::iterator it_atom;
std::wstring headerComment(L" lammps data file created by crysmake");

                xmin=std::numeric_limits<float>::max();//(*ptr_vdb)[0].atoms[0].x;
                xmax=-xmin;
                ymin= xmin;
                ymax=-xmin;
                zmin= xmin;
                zmax=-xmin;

                minmax[0]=&xmin;
                minmax[1]=&xmax;
                minmax[2]=&ymin;
                minmax[3]=&ymax;
                minmax[4]=&zmin;
                minmax[5]=&zmax;


                for(dataBuffer &db: *ptr_vdb){
                        numOfatoms+=db.atoms.size();
                        for(StGrainAtom &atom: db.atoms){

                            if(atom.x<xmin) xmin=atom.x;
                            if(atom.x>xmax) xmax=atom.x;

                            if(atom.y<ymin) ymin=atom.y;
                            if(atom.y>ymax) ymax=atom.y;

                            if(atom.z<zmin) zmin=atom.z;
                            if(atom.z>zmax) zmax=atom.z;

                            it_atom=elementMass.find(atom.name);
                            if(it_atom==elementMass.end()){
                            auto it_elemMass=Elements::mass.find(atom.name);
                            wstring elemMass(L"?");

                                if(it_elemMass==Elements::mass.end()){
                                    wcerr<<" the mass of "<<atom.name<<" is unknown"<<endl;
                                    headerComment+=L"   ERROR: the mass of "+atom.name+L" is unknown";
                                }
                                else
                                    elemMass=it_elemMass->second;

                                 elementMass.insert(std::pair<wstring,wstring>(atom.name,elemMass));
                            }
                       }
                }


                if(dataBuffer::setbox){

                    for(int i=dataBuffer::xmin;i<=dataBuffer::zmax;i++){

                        if( std::fabs(dataBuffer::box[i]) < std::fabs(*minmax[i]) )
                            cout<<"WARNING: box size is less than grain size, component #"<<i
                                <<", box limit: "<<dataBuffer::box[i]<<", grain limit: "<<(*minmax[i])<<endl;
                        else
                            *minmax[i]=dataBuffer::box[i];
                    }

                }
                else{
                    xmin-=dataBuffer::xmargin;
                    xmax+=dataBuffer::xmargin;
                    ymin-=dataBuffer::ymargin;
                    ymax+=dataBuffer::ymargin;
                    zmin-=dataBuffer::zmargin;
                    zmax+=dataBuffer::zmargin;
                }


                /// 1. lammps file start writing

                file<<headerComment<<endl<<endl;

                file<<"\t"<<numOfatoms<<"\tatoms"<<endl;
                file<<"\t"<<elementMass.size()<<"\t atom types"<<endl;

                file<<"    "<<xmin<<"    "<<xmax<<"    xlo xhi"<<endl;
                file<<"    "<<ymin<<"    "<<ymax<<"    ylo yhi"<<endl;
                file<<"    "<<zmin<<"    "<<zmax<<"    zlo zhi"<<endl;

                file<<"  Masses"<<endl<<endl;


std::map<std::wstring, size_t>  nameID; /// mapa skojarzen:  nazwaPierwiastka,  pozycja w elementMass
std::map<std::wstring, size_t>::iterator iter_nameID;

                for(i=1,it_atom=elementMass.begin();it_atom!=elementMass.end();i++,it_atom++){
                        file<<"\t"<<i<<"\t"<<it_atom->second<<"    # "<<it_atom->first<<endl;
                        nameID.insert(std::pair<wstring,size_t>(it_atom->first,i));
                }


                file<<"  Atoms"<<endl<<endl;
                i=1;  // an atom's number

                for(dataBuffer &db : *ptr_vdb)
                    for(StGrainAtom &atom: db.atoms){
                        iter_nameID=nameID.find(atom.name);
                        file<<i++<<" "<<iter_nameID->second<<" "<<atom.x<<" "<<atom.y<<" "<<atom.z<<endl;
                    }

                file.close();
return true;
}
//-----------------------------------------------------------------------------
bool loadAtomsFromXYZ(dataBuffer &db, bool verbose)
{
const string fileName=db.fileName;
wifstream file(fileName.c_str(),ios::in);

                if(!file){
                    cerr<<endl<<"ERROR: couldn't find the input file "<<endl;
                    file.close();
                return false;
                }
                else
                    if(verbose) cout<<"::: file reading: "<<fileName<<endl;

wstring atname;
position x,y,z;
position minZ,maxZ,minX,maxX,minY,maxY;
//position x2max,y2max,z2max;
size_t i,size=0;
//wchar_t cmd[128];

                file.exceptions(ios::badbit );

                try{

                    file>>size;
                    while (file.get()!='\n' && !file.eof())
                    ;

                    while(file.get()!='\n' && !file.eof())
                    ;


                    if(!size){
                        cerr<<" ERROR, empty set"<<endl;
                        throw 0;
                    }

                    if(!file.good()){
                        cerr<<" ERROR, file's header is not good"<<endl;
                    throw 0;
                    }


                    db.atoms.reserve(size);
                    db.grainprm.atomsNameNumber.reserve(size);

                    file>>atname>>x>>y>>z;
                    db.atoms.push_back(StGrainAtom(atname,x,y,z));
                    db.grainprm.atomsNameNumber.push_back(atname);
                    auto annRef=db.grainprm.atomsNameNumber.back();

                    minX=x;
                    maxX=x;

                    minY=y;
                    maxY=y;

                    minZ=z;
                    maxZ=z;

                    file.clear();
                    file.exceptions(ios::failbit);

                    for(i=1;i<size;i++){
                        file>>atname>>x>>y>>z;

                       // if(!file.good()){
                        //    cerr<<"ERROR, file is not good (improper number of atoms, conversion, ?), line: "<<i+3<<endl;
                       // throw 0;
                       // }


                        if(atname!=annRef.name){
                            db.grainprm.atomsNameNumber.push_back(atname);
                            annRef=db.grainprm.atomsNameNumber.back();
                        }


                        db.atoms.emplace_back(StGrainAtom(atname,x,y,z));


                        if(x>maxX) maxX=x;
                        if(x<minX) minX=x;

                        if(y>maxY) maxY=y;
                        if(y<minY) minY=y;

                        if(z>maxZ) maxZ=z;
                        if(z<minZ) minZ=z;

                    }


                    db.grainprm.numOfatoms=size;

                    db.grainprm.minX=minX;
                    db.grainprm.maxX=maxX;

                    db.grainprm.minY=minY;
                    db.grainprm.maxY=maxY;

                    db.grainprm.minZ=minZ;
                    db.grainprm.maxZ=maxZ;


                    /// UWAGA: ponisze procedury zapamietuja nazwy typow atomow

                    auto sortWstr=[](StGrainPrm::StAtomNameNumber  &a, StGrainPrm::StAtomNameNumber &b){return a.name<(b.name);};
                    sort(db.grainprm.atomsNameNumber.begin(), db.grainprm.atomsNameNumber.end(),sortWstr);

                    auto cmpWstr=[](StGrainPrm::StAtomNameNumber  &a, StGrainPrm::StAtomNameNumber &b){return !a.name.compare(b.name);};
                    std::vector<StGrainPrm::StAtomNameNumber>::iterator it;
                    it = std::unique (db.grainprm.atomsNameNumber.begin(), db.grainprm.atomsNameNumber.end(),cmpWstr);

                    db.grainprm.atomsNameNumber.resize( std::distance(db.grainprm.atomsNameNumber.begin(),it) );

                    ///


                    if(verbose){
                        cout<<"::: "<<setw(2)<<std::right<<"*"<<setw(10)<<std::right<<"min"<<setw(10)<<std::right<<"max"<<endl;
                        cout<<"::: "<<setw(2)<<std::right<<"X"<<setw(10)<<std::right<<minX<<setw(10)<<std::right<<maxX<<endl;
                        cout<<"::: "<<setw(2)<<std::right<<"Y"<<setw(10)<<std::right<<minY<<setw(10)<<std::right<<maxY<<endl;
                        cout<<"::: "<<setw(2)<<std::right<<"Z"<<setw(10)<<std::right<<minZ<<setw(10)<<std::right<<maxZ<<endl;
                    }


                }
                catch(...){
                        file.close();

                return false;
                }

                file.close();

                //--------------------------------------------
                db.atoms.shrink_to_fit();
                db.grainprm.atomsNameNumber.shrink_to_fit();

                    for(auto &atom: db.atoms)
                    {
                            for(auto &ann: db.grainprm.atomsNameNumber)
                                if(ann.name==atom.name)
                                    ann.count++;

                    }



                    if(verbose){
                        wcout<<"\n\n::: atom types listing\n::: type number\n";

                        for(auto &ann: db.grainprm.atomsNameNumber)
                            wcout<<"   "<<ann.name<<"\t"<<ann.count<<endl;


                        wcout<<"::: reading SUCCESS\n"<<endl;
                    }

return true;
}
//-----------------------------------------------------------------------------
bool loadAtomsFromLammps(dataBuffer &db, bool verbose)
{
const string fileName=db.fileName;
fstream file(fileName.c_str(),ios::in);

                if(!file){
                    cerr<<endl<<"ERROR: couldn't find the input file "<<endl;
                    file.close();
                return false;
                }
                else
                    if(verbose) cout<<"::: file reading: "<<fileName<<endl;

string line;
string tmp,atype,amass;
position minZ,maxZ,minX,maxX,minY,maxY;
size_t numOfatoms=0,numOftypes=0;


                minX=minY=minZ=INFINITY;
                maxX=maxY=maxZ=-minX;


                try{

                    file.exceptions(ios::badbit | ios::failbit );

                    while (file.get()!='\n' && !file.eof())  ;
                    while (file.get()!='\n' && !file.eof())  ;

                    ///  1. num of atoms, types
                    file>>numOfatoms; while (file.get()!='\n' && !file.eof())  ;
                    file>>numOftypes; while (file.get()!='\n' && !file.eof())  ;

                    if(!numOfatoms || !numOftypes) {
                        cout<<" number of atoms/types error \n";
                    throw 0;
                    }


                    /// 2. ignore box size, save types

                    do{
                        std::getline(file,line);
                    }
                    while (line.find("Masses")==string::npos) ;

                    while (file.get()!='\n' && !file.eof())  ;

                    for(size_t i=0;i<numOftypes;i++){

                        file>>atype>>amass;

                        if(amass.find('?')!=string::npos)
                            cout<<"WARNING: unspecified mass of atom type: "<<atype<<endl;

                    wstring watype(to_wstring(atype));
                    wstring wamass(to_wstring(atype));

                        Elements::addElement(watype,wamass);

                    }


                    /// 3. find "Atoms" section

                    do{
                        std::getline(file,line);
                    }
                    while (line.find("Atoms")==string::npos) ;

                    while (file.get()!='\n' && !file.eof())  ;

                    db.atoms.resize(numOfatoms);


                    /// 4. read atom positions


                    for(StGrainAtom & atom: db.atoms){


                        file>>tmp>>tmp;
                        atom.name=to_wstring(tmp);

                        file>>atom.x>>atom.y>>atom.z;

                        while (file.get()!='\n' && !file.eof())  ;

                        if(atom.x>maxX) maxX=atom.x;
                        if(atom.x<minX) minX=atom.x;

                        if(atom.y>maxY) maxY=atom.y;
                        if(atom.y<minY) minY=atom.y;

                        if(atom.z>maxZ) maxZ=atom.z;
                        if(atom.z<minZ) minZ=atom.z;

                    }


                    db.grainprm.numOfatoms=numOfatoms;

                    db.grainprm.minX=minX;
                    db.grainprm.maxX=maxX;

                    db.grainprm.minY=minY;
                    db.grainprm.maxY=maxY;

                    db.grainprm.minZ=minZ;
                    db.grainprm.maxZ=maxZ;


                    file.close();
                }

                catch(...){

                    file.close();
                return false;
                }



            if(verbose){
                wcout<<"\n\n::: atom types listing\n::: type number\n";

                for(auto &ann: db.grainprm.atomsNameNumber)
                    wcout<<"   "<<ann.name<<"\t"<<ann.count<<endl;


                wcout<<"::: reading SUCCESS\n"<<endl;
            }




return true;
}

//-----------------------------------------------------------------------------
bool saveDataToDMP(const string &filename, vdb *ptr_vdb)
{
fstream foutDMP(filename,ios::out);


                if(!foutDMP){
                    cerr<<endl<<"ERROR: couldn't open file for saving"<<endl;
                    foutDMP.close();
                return false;
                }


float xmin,xmax,ymin,ymax,zmin,zmax;
streampos spos;
size_t numOfatoms=0;



                for(dataBuffer &db: *ptr_vdb)
                        numOfatoms+=db.atoms.size();



                xmin=std::numeric_limits<float>::max();
                xmax=-xmin;
                ymin= xmin;
                ymax=-xmin;
                zmin= xmin;
                zmax=-xmin;

                foutDMP<<"ITEM: TIMESTEP"<<endl;
                foutDMP<<"0"<<endl;
                foutDMP<<"ITEM: NUMBER OF ATOMS"<<endl;
                foutDMP<<numOfatoms<<endl;
                foutDMP<<"ITEM: BOX BOUNDS pp pp pp"<<endl;

                spos=foutDMP.tellp();

                foutDMP<<setw(75)<<" "<<endl; // 3x foutDMP<< StAtomAveStd::minX <<" "<<setprecision(6)<<setw(9)<<StAtomAveStd::maxX<<endl;
                foutDMP<<"ITEM: ATOMS id type x y z"<<endl;

                for(dataBuffer &db: *ptr_vdb){
                const size_t dbsize=db.atoms.size();


                    for(size_t i=0;i<dbsize;i++){
                    const StGrainAtom &atom=db.atoms[i];

                        foutDMP<<setw(3)<<i<<" 1 "<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;

                        if(atom.x<xmin) xmin=atom.x;
                        if(atom.x>xmax) xmax=atom.x;

                        if(atom.y<ymin) ymin=atom.y;
                        if(atom.y>ymax) ymax=atom.y;

                        if(atom.z<zmin) zmin=atom.z;
                        if(atom.z>zmax) zmax=atom.z;

                    }

                }


                foutDMP.seekp(spos,ios::beg);
                foutDMP.precision(6);
                foutDMP.width(9);

                foutDMP<< xmin <<" "<< xmax <<endl<<" ";
                foutDMP<< ymin <<" "<< ymax <<endl<<" ";
                foutDMP<< zmin <<" "<< zmax;




                foutDMP.close();

return true;
}

//-----------------------------------------------------------------------------

bool saveDataToD2S(const string &filename, vdb *ptr_vdb, bool verbose)
{
                if(verbose)
                    cout<<"save data as: "<<filename<<endl;

wfstream fout(filename,ios::out);

                if(!fout){
                    cerr<<endl<<"ERROR: couldn't open file for saving"<<endl;
                    fout.close();
                return false;
                }


size_t numOfatoms=0;
float xmin,ymin,xmax,ymax,zmin,zmax,bwidth;
vector<wstring> anames;
auto tnow = std::chrono::system_clock::now();
auto tdate = std::chrono::system_clock::to_time_t(tnow);

streampos widthnm,backnm;


                xmin=ymin=INFINITY;
                xmax=ymax=-xmin;
                zmin= xmin;
                zmax=-xmin;

                anames.push_back(ptr_vdb->data()->atoms.data()->name);



                for(dataBuffer &db: *ptr_vdb){
                        numOfatoms+=db.atoms.size();
                        for(StGrainAtom &atom: db.atoms){

                            if(atom.x<xmin) xmin=atom.x;
                            if(atom.y<ymin) ymin=atom.y;

                            if(atom.x>xmax) xmax=atom.x;
                            if(atom.y>ymax) ymax=atom.y;

                            if(atom.z<zmin) zmin=atom.z;
                            if(atom.z>zmax) zmax=atom.z;

                            auto iter=std::find(anames.begin(),anames.end(),atom.name);

                            if(iter==anames.end())
                                anames.push_back(atom.name);
                       }
                }

                bwidth=xmax-xmin;

                if(float yminmax=ymax-ymin; yminmax>bwidth)
                    bwidth=yminmax;

                fout<<"#dateTime: "<<std::ctime(&tdate);
                fout<<"#target: "; for(auto &a: anames) fout<<"  "<<a;  fout<<endl;
                fout<<"#shape: brick"<<endl;
                fout<<"#width_nm: "<<(bwidth)*0.1f<<endl;
                fout<<"#front_nm: 0"<<endl;
                fout<<"#back_nm: 0"<<endl; //backnm=fout.tellp();  fout<<"          "<<endl;
                fout<<"#thickness_nm: "<<(zmax-zmin)*0.1f<<endl;
                fout<<"#klm: 001"<<endl;
                fout<<"#numberOfatoms: "<<numOfatoms<<endl;
                fout<<"#\n";
                fout<<"#ele.\tx(pm)"<<"\t"<<setw(9)<<" y(pm)"<<"\t"<<setw(9)<<" z(pm)"<<endl;



                for(dataBuffer &db : *ptr_vdb)
                    for(StGrainAtom &atom: db.atoms){
                        fout<<atom.name<<"\t"<<
                              std::setprecision(6)<<std::setw(9)<<100*(atom.x-xmin)<<"\t"<<
                              std::setprecision(6)<<std::setw(9)<<100*(atom.y-ymin)<<"\t"<<
                              std::setprecision(6)<<std::setw(9)<<100*(atom.z-zmin)<<endl;
                    }

                fout.close();

return true;
}
