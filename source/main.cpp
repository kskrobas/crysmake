//#include <QCoreApplication>

#include <regex>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <string>



#include "clcmdlist.h"

//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{

            if(argc==1){

                cerr<<" usage:   crysmake {configuration file} [-v -h] [-if {arg}] \n";
                cerr<<" -v  - verbose mode"<<endl;
                cerr<<" -h  - help, script commands listing"<<endl;
                cerr<<" -if {arg} - conditional processing argument\n";
                cerr<<" date: "<<__DATE__<<endl;

            return RETVALUE::NARGS_ERR;
            }


std::string configFileName;
std::vector<std::string> cparg;
vfunct funcList;

int i=1;
size_t pos;
bool verbose=false;


                while(i< argc) {
                std::string sarg(argv[i]);

                if(sarg=="-v"){
                    verbose=true;
                    dataBuffer::verbose=true;
                    i++;
                    continue;
                }

                if(argv[i][0]!='-'){
                    configFileName=std::string(argv[i]);
                    pos=configFileName.find_last_of('/');
                std::string dirName(configFileName.substr(0,pos));

                        if(chdir(dirName.c_str())==-1 && pos!=string::npos){
                            cerr<<" ERROR , specified dirName doesn't exist: "<<dirName<<"\n";
                        return  RETVALUE::CFILE_ERR;
                        }

                        if(verbose){
                            cout<<"current directory ";
                            system("pwd");
                        }

                        configFileName=configFileName.substr(pos+1);
                fstream file(configFileName);


                        if(!file){
                            cerr<<" ERROR , specified script file doesn't exist: "<<configFileName<<endl;
                        return  RETVALUE::CFILE_ERR;
                        }

                        file.close();

                        i++;
                        continue;
                 }

                 if(sarg=="-l" || sarg=="-h"){
                     cout<<"\n list of script file commands accepted by the program\n";
                     cout<<"\n -------------------------------------------\n";
                     for(size_t i=0;i<commandListHelp::clhSize;i++)
                         cout<<commandListHelp::clh[i]<<"\n";

                     cout<<"\n============================================"<<endl;
                     i++;

                     if(argc==2)
                         return RETVALUE::OK;
                     continue;
                 }


                 if(sarg=="-if"){
                    i++;

                    if(i<argc){
                        cparg.push_back(std::string(argv[i]));
                        i++;
                        continue;
                    }
                    else{
                        cerr<<"ERROR, \'if\' option needs an argument\n";
                    return RETVALUE::UNKOPT;
                    }


                 }

                         cerr<<" unknown option "<<endl;

                return RETVALUE::UNKOPT;
                }

                //-----------------------------------------------
                //config file processing

                if(verbose){
                    cout<<"::: working directory: ";cout.flush();
                    system("pwd"); cout<<endl;
                    cout<<"::: == reading script file: "<<configFileName<<" =="<<endl;
                }

                if(configFileName.empty()){
                    cerr<<"ERROR, script file is not given"<<endl;
                return RETVALUE::CFILE_ERR;
                }


vector<ClKeyValues> cmdlist;
vstruct structList;
fstream fileConfig(configFileName,ios::in);
size_t cmdLineCounter=0,sizeDataBuffer=0,sizeStack=0,popNumber=0;


                cmdlist.reserve(32);
                fileConfig.exceptions(ios::badbit);


StInputParametrs ip;
                ip.cmdlist=&cmdlist;
                ip.cmdLineCounter=&cmdLineCounter;
                ip.sizeDataBuffer=&sizeDataBuffer;
                ip.sizeStack=&sizeStack;
                ip.popNumber=&popNumber;
                ip.cparg=&cparg;
                ip.verbose=verbose;
                ip.processing=true;
                ip.funct=&funcList;
                ip.structPrm=&structList;


RETVALUE::eRetValue erv=cmdListProcessing(fileConfig,ip);

                if(erv!=RETVALUE::OK){

                    if(erv==RETVALUE::C_END)
                        cerr<<"ERROR: misplaced (orphaned) end stateman"<<endl;
                    if(erv==RETVALUE::C_ENDF)
                        cerr<<"ERROR: misplaced (orphaned) endf stateman"<<endl;


                    cerr<<"ERROR, code="<<erv<<endl;

                    fileConfig.close();
                    return erv;
                 }

                if(verbose)
                    cout<<"::: == scirpt file finished reading ==\n\n";


                //-------------------------------------------------------------


vdb vdataBuffer,vstackBuffer;

                 vdataBuffer.resize(sizeDataBuffer+popNumber);
                 vstackBuffer.resize(sizeStack);


auto ptr_dataBuffer=vdataBuffer.data();
auto ptr_stackBuffer=vstackBuffer.data();
string key;


                if(verbose)
                    cout<<"+ == commands processing \n\n";



                for(ClKeyValues &keyvalue : cmdlist){
                    keyvalue>>key;

                    if(verbose)
                        cout<<"+ "<<key<<endl;


                    if(key=="center"){
                    string mode;

                            keyvalue>>mode;

                            if(mode=="loc")
                                (ptr_dataBuffer-1)->center();
                            else{
                            float cx,cy,cz;
                                    cx=cy=cz=0;
                            size_t numOfatoms=0;

                                    for(dataBuffer &db : vdataBuffer)
                                        for(StGrainAtom &atom:db.atoms){
                                            cx+=atom.x;
                                            cy+=atom.y;
                                            cz+=atom.z;
                                            numOfatoms++;
                                        }


                                    if(!numOfatoms){
                                        cerr<<"+ WARNING: centering is impossible, number of atoms is equal 0\n";
                                    continue;
                                    }


                            const float inumOfatoms=-1.0f/numOfatoms;

                                    cx*=inumOfatoms;
                                    cy*=inumOfatoms;
                                    cz*=inumOfatoms;

                                    if(verbose){
                                        cout<<"+ translation vector of centering: "<<cx<<", "<<cy<<", "<<cz<<endl;
                                    }

                                    for(dataBuffer &db : vdataBuffer)
                                        for(StGrainAtom &atom:db.atoms){
                                            atom.x+=cx;
                                            atom.y+=cy;
                                            atom.z+=cz;
                                        }
                            }
                    }


                    if(key=="create"){
                    string structname;
                            keyvalue>>structname;

                       continue;
                    }


                    if(key=="infile"){
                    string filename;
                    //string ksort;

                        keyvalue>>filename;

                        ptr_dataBuffer->fileName=filename;


                        if(verbose)
                            cout<<"::: "<<ptr_dataBuffer->fileName<<endl;


                        if(filename.find(".xyz")!=string::npos){
                            if(!loadAtomsFromXYZ(*ptr_dataBuffer,verbose))
                            return RETVALUE::IFILE_ERR;
                        }
                        else{

                            if(filename.find(".lmp")!=string::npos){
                                if(!loadAtomsFromLammps(*ptr_dataBuffer,verbose))
                                    return RETVALUE::INC_FILE_ERR;
                            }
                            else
                                if(!loadAtomsFromNDL(*ptr_dataBuffer,verbose))
                                return RETVALUE::IFILE_ERR;
                        }


                        ptr_dataBuffer++;
                        continue;
                    }


                    if(key=="intspace"){
                    string xyz,pos,space;

                            keyvalue>>xyz>>pos>>space;
                            (ptr_dataBuffer-1)->intspace(xyz,pos,space);

                            continue;
                    }


                    if(key=="marker"){
                    string apos,tol,mode,sig,sname;

                            keyvalue>>apos>>tol>>mode>>sig;
                            (ptr_dataBuffer-1)->markerLayer(apos,tol,mode,sig);
                            continue;
                    }

                    if(key=="mass"){
                    string mode,mass;
                    auto ptr_currData=ptr_dataBuffer-1;

                            keyvalue>>mode>>mass;

                    wstring wmass(mass.begin(),mass.end());

                            if(mode=="A")
                                ptr_currData->grainprm.atomMassA=wmass;
                            else
                                ptr_currData->grainprm.atomMassB=wmass;


                            continue;
                    }

                    if(key=="merge"){
                    //no arguments

                    const size_t dbSize=(ptr_dataBuffer-vdataBuffer.data());

                            if(dbSize>1){
                            vdb m_dataBuffer(sizeDataBuffer-dbSize+1);
                            vatoms::iterator m_dataPtr;
                            size_t numOfatoms=0;


                                        sizeDataBuffer=m_dataBuffer.size();

                                        for(auto ptr=vdataBuffer.data();ptr!=ptr_dataBuffer;ptr++)
                                            numOfatoms+=ptr->atoms.size();

                                        m_dataBuffer.data()->atoms.resize(numOfatoms);
                                        m_dataPtr= m_dataBuffer.data()->atoms.begin();

                                        for(auto ptr=vdataBuffer.data();ptr!=ptr_dataBuffer;ptr++){
                                                std::copy(ptr->atoms.begin(),ptr->atoms.end(),m_dataPtr);
                                                m_dataPtr+=ptr->atoms.size();
                                        }


                                        vdataBuffer=std::move(m_dataBuffer);
                                        ptr_dataBuffer=vdataBuffer.data();



                            coordinate minX,maxX,minY,maxY,minZ,maxZ;

                                        minX=minY=minZ=std::numeric_limits<float>::max();
                                        maxX=maxY=maxZ=-minX;


                                        for(auto & data : vdataBuffer)
                                            for(StGrainAtom & atom: data.atoms){

                                                if(atom.x>maxX) maxX=atom.x;
                                                if(atom.x<minX) minX=atom.x;

                                                if(atom.y>maxY) maxY=atom.y;
                                                if(atom.y<minY) minY=atom.y;

                                                if(atom.z>maxZ) maxZ=atom.z;
                                                if(atom.z<minZ) minZ=atom.z;

                                            }


                                        ptr_dataBuffer->grainprm.minX=minX;
                                        ptr_dataBuffer->grainprm.maxX=maxX;

                                        ptr_dataBuffer->grainprm.minY=minY;
                                        ptr_dataBuffer->grainprm.maxY=maxY;

                                        ptr_dataBuffer->grainprm.minZ=minZ;
                                        ptr_dataBuffer->grainprm.maxZ=maxZ;

                                        ptr_dataBuffer=vdataBuffer.data()+1;

                            }


                    continue;
                    }


                    /*
                    if(key=="minmax"){
                    coordinate minX,maxX,minY,maxY,minZ,maxZ;

                                //minX=minY=minZ=std::numeric_limits<float>::max();
                                maxX=maxY=maxZ=-minX;

                                for(auto & data : vstackBuffer)
                                    for(StGrainAtom & atom: data.atoms){

                                        if(atom.x>maxX) maxX=atom.x;
                                        if(atom.x<minX) minX=atom.x;

                                        if(atom.y>maxY) maxY=atom.y;
                                        if(atom.y<minY) minY=atom.y;

                                        if(atom.z>maxZ) maxZ=atom.z;
                                        if(atom.z<minZ) minZ=atom.z;

                                    }


                                minX=ptr_dataBuffer->grainprm.minX;
                                maxX=ptr_dataBuffer->grainprm.maxX;

                                minY=ptr_dataBuffer->grainprm.minY;
                                maxY=ptr_dataBuffer->grainprm.maxY;

                                minZ=ptr_dataBuffer->grainprm.minZ;
                                maxZ=ptr_dataBuffer->grainprm.maxZ;


                                cout<<"minX, maxX: "<<minX<<", "<<maxX<<"\n";
                                cout<<"minY, maxY: "<<minY<<", "<<maxY<<"\n";
                                cout<<"minZ, maxZ: "<<minZ<<", "<<maxZ<<"\n";

                                continue;
                    }
                    */


                    if(key=="outfile"){
                    string outputfile;
                            keyvalue>>outputfile;

                    vector<string> fnTokens(split<string>(outputfile,"."));
                    const string fileExt(fnTokens.back());


                                try{
                                    if(fileExt.find("xyz")!=string::npos){
                                        throw saveDataToXYZ(outputfile,&vdataBuffer) ;
                                     }

                                    if(fileExt.find("poly")!=string::npos )
                                        throw saveDataToCONFIG(outputfile,&vdataBuffer);

                                    if(fileExt.find("lmp")!=string::npos )
                                        throw saveDataToLammps(outputfile,&vdataBuffer);

                                    if(fileExt.find("d2s")!=string::npos)
                                        throw saveDataToD2S(outputfile,&vdataBuffer,verbose);

                                    if(fileExt.find("dmp")!=string::npos )
                                        throw saveDataToDMP(outputfile,&vdataBuffer);

                                    cerr<<" unrecognized file format "<<endl;
                                    throw false;

                                }
                                catch(bool b){
                                    if(verbose){
                                        cout<<"::: file: "<<outputfile;
                                        cout<<((b)? "  success":"  failure")<<endl;
                                    }
                                }

                        continue;
                    }

                    if(key=="pause"){
                    //no arguments

                    string cmd("read -p \"pause, press any key to continue ...\"");
                    const int retval=system(cmd.c_str());
                            if(verbose)
                                cout<<"retval "<<retval<<endl;

                            continue;
                    }


                    if(key=="pop"){
                    string dataName;
                        keyvalue>>dataName;

                    auto currPtr=vstackBuffer.data();


                            for(auto & data : vstackBuffer)
                                if(data.dataName==dataName)
                                    break;
                                else currPtr++;

                            (*ptr_dataBuffer)=(*currPtr);
                            ptr_dataBuffer++;
                            continue;
                    }

                    if(key=="push"){
                    string dataName;

                        (*ptr_stackBuffer)=(*(ptr_dataBuffer-1));
                        keyvalue>>ptr_stackBuffer->dataName;
                        ptr_stackBuffer++;
                        continue;
                    }


                    if(key=="print"){
                    vector<string> vkeyvalues;

                            keyvalue>>vkeyvalues;

                    const size_t vsize=vkeyvalues.size();
                            for(size_t i=1;i<vsize;i++)
                                (ptr_dataBuffer-1)->printStatistic(vkeyvalues[i]);

                        continue;
                    }

                    if(key=="printf"){
                    string argument;

                            keyvalue>>argument;
                    const string cmd("printf "+argument);

                                system(cmd.c_str());
                        continue;
                    }

                    if(key=="randpos"){
                    string mode,sxr,syr,szr;
                        keyvalue>>mode>>sxr>>syr>>szr;

                    cfloat rx=std::stof(sxr);
                    cfloat ry=std::stof(syr);
                    cfloat rz=std::stof(szr);

                        if(mode=="loc")
                            (ptr_dataBuffer-1)->randDislocation(rx,ry,rz);
                        else{
                            for(dataBuffer &db : vdataBuffer)
                                db.randDislocation(rx,ry,rz);
                        }

                        continue;
                    }


                    if(key=="remove"){//remove plane out 0 0 1 -10
                    string rtype,mode,sa,sb,sc,sd,se,sf;

                            keyvalue>>rtype>>mode>>sa>>sb>>sc>>sd>>se>>sf;

                            if(rtype=="plane"){
                            const bool bm=(mode=="out");
                            const float a=std::stof(sa);
                            const float b=std::stof(sb);
                            const float c=std::stof(sc);
                            string fileName(se);


                                    if(sd=="rand"){
                                    const float d=std::stof(se);
                                    const float e=std::stof(sf);

                                        (ptr_dataBuffer-1)->removeAtomsPlaneRand(a,b,c,d,e,bm);
                                    }
                                    else{                                                                        
                                    const float d=std::stof(sd);

                                        (ptr_dataBuffer-1)->removeAtomsPlane(a,b,c,d,bm,fileName);
                                    }
                            }
                            else{

                                if(rtype=="cyl"){
                                const bool bm=(sa=="out");
                                const float r=std::stof(sb);
                                const float h=(sc.empty())? std::stof("inf"): std::stof(sc);

                                    (ptr_dataBuffer-1)->removeAtomsCylXY(r,bm,h);
                                }
                                else{

                                    if(rtype=="cuboid"){
                                    const bool bm=(mode=="out");
                                    const float a=std::stof(sa);
                                    const float b=std::stof(sb);

                                        (ptr_dataBuffer-1)->removeAtomsCuboid(a,b,bm);

                                    }
                                    else
                                            if(rtype=="slice"){//mode.find("out")!=string::npos  || mode.find("in")!=string::npos
                                                if(regex_match(mode,regex ("(out|in)"))){
                                                const bool bm=(mode=="out");
                                                const float a=std::stof(sa);
                                                const float b=std::stof(sb);
                                                const float c=std::stof(sc);
                                                const float d=std::stof(sd);
                                                const float e=std::stof(se);

                                                    (ptr_dataBuffer-1)->removeAtomsSlice(a,b,c,d,e,bm);
                                                }
                                                else{//[xyz](min|max)
                                                const char xyz=mode[0];
                                                const bool bmax=(mode.find("max")!=string::npos);
                                                const float a=std::stof(sa);


                                                    (ptr_dataBuffer-1)->removeAtomsSlice(xyz,bmax,a);

                                                }
                                            }
                                            else
                                                if(rtype=="random"){
                                                const float a=std::stof(sa);

                                                    (ptr_dataBuffer-1)->removeAtomsRandom(a);

                                                }
                                                else{




                                                    //sphere
                                                const bool bm=(mode=="out");
                                                const float r=std::stof(sa);

                                                    (ptr_dataBuffer-1)->removeAtomsSph(r,bm);
                                                }
                                }
                            }


                    continue;
                    }


                    if(key=="rename"){
                    string sfrom,sto;

                            keyvalue>>sfrom>>sto;

                    wstring wfrom(sfrom.begin(),sfrom.end());
                    wstring wto(sto.begin(),sto.end());

                            (ptr_dataBuffer-1)->renameAtoms(wfrom,wto);

                    continue;
                    }



                    if(key=="rescale"){
                    string value;
                    float fx,fy,fz;


                            if(keyvalue.numOfvalues()==1){
                                keyvalue>>fz;
                                fx=fy=fz;
                            }
                            else
                                keyvalue>>fx>>fy>>fz;

                            (ptr_dataBuffer-1)->rescalePosition(fx,fy,fz);

                            continue;
                    }


                    if(key=="rotate"){
                    string  sa,sx,sy,sz;
                            keyvalue>>sa>>sx>>sy>>sz;
                    const float a=std::stof(sa);
                    const float x=std::stof(sx);
                    const float y=std::stof(sy);
                    const float z=std::stof(sz);

                        (ptr_dataBuffer-1)->rotateAtoms(a,x,y,z);
                        continue;
                    }

                    if(key=="slices"){
                    string dirName,sstart,sstep,sstop;
                            keyvalue>>dirName>>sstart>>sstep>>sstop;


                    const float start=std::stof(sstart);
                    const float step=std::stof(sstep);
                    const float stop=std::stof(sstop);




                        (ptr_dataBuffer-1)->slices(dirName,start,step,stop);

                    continue;
                    }


                    if(key=="sort"){
                    string sorder;

                            keyvalue>>sorder;
                    bool asc=(sorder =="a" || sorder=="A");

                            (ptr_dataBuffer-1)->sortNames(asc);
                            continue;
                    }


                    if(key =="stop")
                        break;

                    if(key=="system"){
                    string arguments;
                            keyvalue>>arguments;

                    const int retval=system(arguments.c_str());
                            if(verbose)
                                cout<<"retval "<<retval<<endl;

                            continue;
                    }

                    if(key=="translate"){
                    string  sx,sy,sz;
                            keyvalue>>sx>>sy>>sz;
                    const float x=std::stof(sx);
                    const float y=std::stof(sy);
                    const float z=std::stof(sz);

                        (ptr_dataBuffer-1)->translateAtoms(x,y,z);
                        continue;
                    }

                    if(key.find("space")==1){
                    string value,valueInd;
                        keyvalue>>value>>valueInd;

                    auto ptr_prevData=(valueInd.empty()) ? ptr_dataBuffer-2 : vdataBuffer.data() + std::stoi(valueInd) ;

                                if(ptr_prevData<vdataBuffer.data()){
                                    cerr<<" IGNORED zspace "<<value<<" operation \n";
                                    cerr<<" zspace needs at least 2 elements in data buffer"<<endl;
                                continue;
                                }

                    auto ptr_currData=ptr_dataBuffer-1;

                    const int mode= (key[0]-'x');
                    const int emin= mode*2;
                    const int emax= mode*2+1;

                    const float space=std::stof(value);
                    const float pminmax=(space>0) ? ptr_prevData->grainprm.minmax[emax] : ptr_prevData->grainprm.minmax[emin];
                    const float cminmax=(space>0) ? ptr_currData->grainprm.minmax[emin] : ptr_currData->grainprm.minmax[emax];

                    const float deltaShift=(space>0) ? space+pminmax-cminmax : space+pminmax-cminmax;


                            for (StGrainAtom &atom: ptr_currData->atoms)
                                atom.xyz.xyz[mode]+=deltaShift;

                            ptr_currData->grainprm.minmax[emin]+=deltaShift;
                            ptr_currData->grainprm.minmax[emax]+=deltaShift;

                        continue;
                    }
                }


return RETVALUE::OK;
}
//


//-------------------------------------------------------------------------------------------------


