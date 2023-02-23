#include "clcmdlist.h"
#include <map>


typedef std::map<std::string, std::string> varmap;


namespace  commandListHelp {
const string clh[]={


                   "*  center {loc|tot}",
                   "*  end - end of conditional processing of current nesting",
                   "*  endf - end of function ",
                   "*  exit - stop script processing",
                   "*  fcall {name} -  call a function",
                   "*  fixbox {xmin xmax ymin ymax zmin zmax} - cuboid box sizes",
                   "*  function {name} - function",
                   "*  if {arg} - conditional processing of following lines, arg is passed by command line options",
                   "*  infile  {the name of input data file | $var}",
                   "*  include {the name of input crysmake file}",
                   "*  intspace {x|y|z} {pos} {space} - creates two subgrains separated by {space} begining at {pos} ",
                   "*  marker {a. pos} {tol} {x|y|z} {signature} - a.pos stands for approximate position",
                   "*  mass {A|B} <value> - set atom A or B mass",
                   "*  merge",
                   "*  minZzero {no|yes} - tranlate atoms along z-axis to fulfill a condition  (x,y,z>=0)",
                   "*  {b|x|y|z}margin <value> - b parameter sets margin for all directions; x, y, z for selected direction only",
                   "*  outfile {the name of output file | $var} - acceptable formats: xyz poly lmp dmp",
                   "*  pause",
                   "*  push {the name of set}",
                   "*  pop  {the name of set}",
                   "*  print <any number, any order of arguments>",
                   "      arguments: [xyz]range , numatoms, rematoms, separator([printableChar]{1}[digit]{0,2})? ",
                   "*  printf \"<string>\"",
                   "*  quit - stop script processing ",
                   "*  randpos (loc|tot) {x|$var y|$var z|$var}",
                   "*  remove plane {out|in} {A B C D} (fileName)?- plane params Ax+By+Cz+D=0, out/in - vector direction",
                   "*  remove slice [xyz](max|min)  {value} - remove a slice relative the most deviated atoms along [xyz] axis",
                   "*  remove slice {out|in} {A B C D thickness} - as above",
                   "*  remove cuboid {out|in} {a b} - a, b lenght of sides along x, y axis respectively",
                   "*  remove cyl {xy} {out|in} {r}",
                   "*  remove sph  {out|in} {r}",
                   "*  remove random {value} - remove atom with propability equal to value given in percent",
                   "*  rename {from} {to} - renaming element",
                   "*  rescale {value}{1,3}",
                   "*  rotate {angle x y z} - angle in deg units",
                   "*  slices {outDirName} {zStart zStep zStop} - save a set of layers in a separate files of selected directory",
                   "*  sort {[aAdD]} - [aA] ascending order",
                   "*  stop - stop command executing",
                   "*  system {command} - start shell program (in background with  &) ",
                   "*  translate {x y z}",
                   "*  [xyz]space  {value} [id]- relative space between the current and the previous(default) or [id] grain",
                   "*  $var=value - variable definition (digits and letters are allowed)",
                   "*  /*    */ - block comment",
                   "*  # or  %   - comments"

                    };


const size_t clhSize=sizeof(clh)/sizeof(string);
}

const std::string strNUMBER(RE_NUMBER);
const std::string strPNUMBER(PRE_NUMBER);

//-----------------------------------------------------------------------------
void ClKeyValues::fini()
{
        ptr_frontValue=&keyvalues.back();
        ptr_lastValue=ptr_frontValue;
        ptr_f=&ClKeyValues::fcon;
        cout<<" f ini "<<endl;
}
//-----------------------------------------------------------------------------
void ClKeyValues::fcon()
{
        ptr_lastValue=&keyvalues.back();

        cout<<" f con "<<endl;
}
//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator <<(const char *ckv)
{
const string kv(ckv);

        if(keyvalues.empty()){
            keyvalues.push_back(kv);
            ptr_frontValue=&keyvalues.back();
         }
        else
            keyvalues.push_back(kv);

        ptr_lastValue=&keyvalues.back();

return *this;
}

ClKeyValues &ClKeyValues::operator <<(const  vector<string>  & vs)
{
            keyvalues=std::move(vs);
            ptr_frontValue=& keyvalues.front();
            ptr_lastValue= & keyvalues.back();

return *this;
}

ClKeyValues &ClKeyValues::operator <<(string & kv)
{
        if(keyvalues.empty()){
            keyvalues.push_back(kv);
            ptr_frontValue=&keyvalues.back();
         }
        else
            keyvalues.push_back(kv);

        ptr_lastValue=&keyvalues.back();


return *this;
}
//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator >>(string &kv)
{

    if(ptr_frontValue <= ptr_lastValue ){
        kv=*ptr_frontValue;
        ptr_frontValue++;
     }
    else
        kv.clear();

    return *this;
}

ClKeyValues &ClKeyValues::operator >>(float &f)
{
    if(ptr_frontValue <= ptr_lastValue ){
        //kv=*ptr_frontValue;
        f=std::stof(*ptr_frontValue);
        ptr_frontValue++;
     }
    else
        f=0;

return *this;
}
//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator >>(vector<string> &vs)
{
        vs=keyvalues;
return *this;
}

ClKeyValues::ClKeyValues(const ClKeyValues &v)
{
       //cout<<"copy constr"<<endl;
        keyvalues=v.keyvalues;
        ptr_frontValue= & keyvalues.front();
        ptr_lastValue= & keyvalues.back();

}

ClKeyValues::ClKeyValues(ClKeyValues &&v)
{
       //cout<<"copy constr"<<endl;
        keyvalues=std::move(v.keyvalues);
        ptr_frontValue= & keyvalues.front();
        ptr_lastValue= & keyvalues.back();

}

//ClKeyValues ClKeyValues::operator =(const ClKeyValues &c)
//{
//    cout<<"copy "<<endl;
//}


//ClKeyValues ClKeyValues::operator =( ClKeyValues &&m)
//{
//    cout<<"move "<<endl;
//}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

std::string encdollar( varmap & vm, std::string &cmd , size_t pos=0);
std::string encdollarY( varmap & vm, std::string &cmd , size_t pos)
{
std::string var("$"),zzz;
//int iend=pos;

            while(  isalnum(cmd[pos]) ){
                var+=cmd[pos++];
            }

std::map<string,string>::iterator it=vm.find(var);

            if(it!=vm.end())
                zzz=it->second;
            else{
                cerr<<"ERROR: unknown variable "<<var<<endl;
                throw RETVALUE::C_UNKVAR;
            }

            if(pos==cmd.length())
            return zzz;

            if(cmd[pos]!='$')
                zzz+=encdollar(vm,cmd,pos);
            else
                zzz+=encdollarY(vm,cmd,pos+1);
return zzz;
}

std::string encdollar(varmap &vm, std::string &cmd , size_t pos)
{
std::string zzz;

            try{

                if(cmd[pos]=='$')
                    zzz+=encdollarY(vm,cmd,pos+1);
                else{

                    while(cmd[pos]!='$' && pos!=cmd.length() ){
                        zzz+=cmd[pos];
                        pos++;
                    }

                    if(cmd[pos]=='$')
                        zzz+=encdollarY(vm,cmd,pos+1);

                }
            }
            catch(RETVALUE::eRetValue e){

                throw e;
            }

return zzz;
}

void replaceVariable(varmap &vm,vector<string> &tokens, const size_t from, const size_t to)
{
std::map<string,string>::iterator it;

            /// if exists replace $variable with proper value
            for(size_t i=from;i<to;i++){

                if(tokens[i][0]=='$'){
                    it=vm.find(tokens[i]);

                    if(it!=vm.end()){
                        if(std::regex_match(it->second,std::regex("[[:s:]]*[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?")))
                            tokens[i]=it->second;
                        else{
                            cerr<<"ERROR, variable is not a number: "<<tokens[i]<<endl;
                            throw RETVALUE::C_NONUM;
                        }
                    }
                    else{
                        cerr<<"ERROR, unknown variable: "<<tokens[i]<<endl;
                        throw RETVALUE::C_UNKVAR;
                    }
                }
            }
}



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//
RETVALUE::eRetValue cmdListProcessing(fstream &file,StInputParametrs &prm)
{

const int recnest=prm.recnest+1;
const bool verbose=prm.verbose;
std::string cmdline;
vector<ClKeyValues> &cmdlist=*prm.cmdlist;
fstream &fileConfig=file;
size_t &cmdLineCounter=*prm.cmdLineCounter;
size_t &sizeDataBuffer=*prm.sizeDataBuffer;
size_t &sizeStack=*prm.sizeStack;
size_t &popNumber=*prm.popNumber;
size_t pos;
varmap variables; /// 'name', 'value'


                    try{

                        while(!fileConfig.eof()){

                            std::getline(fileConfig,cmdline);
                            cmdLineCounter++;

                            if(verbose) cout<<"::: " <<cmdline<<endl;

                            trim(cmdline);//usun biale znaki

                            //--
                            if(cmdline.empty()) continue;
                            //--

                            //--
                            //komentarz
                            if(cmdline[0]=='#' || cmdline[0]=='%') continue;

                            if(cmdline.find("/*")!=string::npos){

                                while(cmdline.find("*/")==string::npos && !fileConfig.eof() ){
                                    std::getline(fileConfig,cmdline);
                                    cmdLineCounter++;
                                }


                            continue;
                            }

                            //--


                            //-----------------------------------------------------------
                            if(std::regex_match(cmdline,std::regex("end"))){
                                return RETVALUE::C_END;
                            }
                            //---

                            if(std::regex_match(cmdline,std::regex("else")))
                            return RETVALUE::C_ELSE;

                            //---

                            if(std::regex_match(cmdline,std::regex("endf"))){
                                return RETVALUE::C_ENDF;
                            }

                            //-----------------------------------------------------------


                            //---
                            if(std::regex_match(cmdline,std::regex("if[[:s:]]+\\w+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            const bool sprocessing=prm.processing;
                            const bool argExists=prm.argumentExists(tokens[1]) ;
                            const bool aprocessing=argExists & sprocessing;

                                        prm.processing=aprocessing;

                            RETVALUE::eRetValue erv=cmdListProcessing(file,prm);

                                        if(erv==RETVALUE::C_END){ /// OK . no errors
                                           prm.processing=sprocessing;
                                        }
                                        else{
                                            if(erv==RETVALUE::C_ELSE){
                                            const bool else_sprocessing= sprocessing;
                                            const bool else_aprocessing= !argExists & else_sprocessing;

                                                        prm.processing=else_aprocessing;

                                            RETVALUE::eRetValue erv=cmdListProcessing(file,prm);

                                                    if(erv==RETVALUE::C_END)
                                                        prm.processing=sprocessing; /// OK . no errors
                                                    else{
                                                        cerr<<"ERROR called by else preprocessing"<<endl;
                                                        if(erv==RETVALUE::PROC_ERR)
                                                            cerr<<" missing \'end\' statement"<<endl;

                                                        return RETVALUE::C_ELSE_ERR;
                                                    }
                                            }
                                            else{
                                                if(erv==RETVALUE::PROC_ERR)
                                                    cerr<<"ERROR, missing \'end\' statement"<<endl;

                                            return RETVALUE::C_END_ERR;
                                            }
                                        }

                                continue;
                            }

                            //---
                            if(std::regex_match(cmdline,std::regex("function[[:s:]]+\\w+[[:s:]]*"))){

                                        if(recnest>0){
                                            cerr<<"ERROR, functions can't be nested"<<endl;
                                        return RETVALUE::C_NEST_FUNC_ERR;
                                        }

                            vector<string> tokens(split<string> (cmdline," "));
                            const string &funcname=tokens[1];
                            const bool sprocessing=prm.processing;

                                        if(prm.functionExists(funcname) && sprocessing){
                                            cerr<<"ERROR, function names \'"<<funcname<<"\' collision"<<endl;
                                            cerr<<"  first defined here: "<<(*prm.funct)[prm.fpos].cmdPos<<endl;
                                        return RETVALUE::ERR_FNAMECOLL;
                                        }


                             const int funcpos=file.tellg();
                             const int tmpCmdLine=cmdLineCounter; //save value of line counter

                                        prm.processing=false;

                             RETVALUE::eRetValue erv=cmdListProcessing(file,prm);

                                        if(erv==RETVALUE::C_ENDF)
                                            prm.processing=sprocessing;
                                         else{

                                            if(erv==RETVALUE::PROC_ERR)
                                                 cerr<<"ERROR, missing \'endf\' statement  "<<endl;

                                         return RETVALUE::C_ENDF_ERR;
                                        }

                              ClKeyValues keyvalues;
                                          keyvalues<<tokens;

                                          cmdlist.emplace_back(keyvalues);
                                          prm.funct->emplace_back(StFunct(funcname,funcpos,tmpCmdLine));

                                          continue;
                            }

                            //---
                            if(std::regex_match(cmdline,std::regex("struct[[:s:]]+\\w+[[:s:]]*"))){

                                            if(recnest>0){
                                                cerr<<"ERROR, struct procedures can't be nested"<<endl;
                                            return RETVALUE::C_NEST_FUNC_ERR;
                                            }


                            vector<string> tokens(split<string> (cmdline," "));
                            const string &funcname=tokens[1];
                            const bool sprocessing=prm.processing;

                                        if(prm.structNameExists(funcname) && sprocessing){
                                            cerr<<"ERROR, create names \'"<<funcname<<"\' collision"<<endl;
                                            //cerr<<"  first defined here: "<<(*prm.funct)[prm.fpos].cmdPos<<endl;
                                        return RETVALUE::ERR_FNAMECOLL;
                                        }

                                        //prm.processing=false;

                             RETVALUE::eRetValue erv=cmdListCreateProcessing(file,prm);


                                        if(erv==RETVALUE::C_ENDC){
                                             prm.processing=sprocessing;
                                             prm.structPrm->back().name=funcname;

                                             sizeDataBuffer++;
                                        }
                                        else{
                                            if(erv==RETVALUE::PROC_ERR)
                                                 cerr<<"ERROR, missing \'endc\' statement  "<<endl;

                                         return RETVALUE::C_ENDC_ERR;
                                        }


                                        continue;
                            }

                            //-----------------------------------------------------------

                            if(!prm.processing){

                                if(verbose)
                                    cout<<"---> the previous command will be not processed !"<<endl;

                                continue;
                            }

                            //-----------------------------------------------------------



                            //--                            
                            if(std::regex_match(cmdline,std::regex("stop")))
                            {
                                ClKeyValues keyvalues;
                                            keyvalues<<"stop";
                                            cmdlist.emplace_back(keyvalues);
                            }



                            //--
                            if(std::regex_match(cmdline,std::regex("fcall[[:s:]]+\\w+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            const string &funcname=tokens[1];

                                        if(!prm.functionExists(funcname)){
                                            cerr<<"ERROR, function named \'"<<funcname<<"\' doesn't exist"<<endl;
                                            cerr<<" line: "<<cmdLineCounter<<endl;
                                        return RETVALUE::C_FUNC_ERR;
                                        }

                             const int &fpos=prm.fpos;
                             const size_t currCmdLineCounter=cmdLineCounter;
                             streampos currFilePos=-1;

                                        if(!file.eof())
                                            currFilePos=file.tellg();
                                        else
                                            file.clear();

                                        file.seekg( (*prm.funct)[fpos].filePos);
                                        cmdLineCounter=(*prm.funct)[fpos].cmdPos;

                             RETVALUE::eRetValue erv=cmdListProcessing(file,prm);

                                         if(erv!=RETVALUE::C_ENDF){
                                             return erv;
                                         }

                                         if(currFilePos<0)
                                            file.setstate(ios::eofbit); // if end of file
                                          else
                                            file.seekg(currFilePos);// if not end of file

                                         cmdLineCounter=currCmdLineCounter;
                                        continue;
                            }
                            //--
                            if(std::regex_match(cmdline,std::regex("create[[:s:]]+\\w+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));

                                            if(!prm.structNameExists(tokens[1])) {
                                                cerr<<"ERROR, structure named \'"<<tokens[1]<<"\' doesn't exist"<<endl;
                                                cerr<<" line: "<<cmdLineCounter<<endl;
                                            return RETVALUE::C_CREATE_ERR;
                                            }

                           ClKeyValues keyvalues;

                                           keyvalues<<tokens;
                                           cmdlist.emplace_back(keyvalues);

                                    continue;
                             }

                            //--
                            if(std::regex_match(cmdline,std::regex("center[[:s:]]+(loc|tot)[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                           keyvalues<<tokens;
                                           cmdlist.emplace_back(keyvalues);

                                    continue;
                             }                                                                                                             
                            //--
                            if(std::regex_match(cmdline,std::regex("exit[[:s:]]*"))){
                            //vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                           keyvalues<<"exit";
                                           cmdlist.emplace_back(keyvalues);

                                    break;
                            }
                            //--
                            if(std::regex_match(cmdline,std::regex("include[[:s:]]+[[:print:]]+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            string fileName;

                                    if(tokens[1][0]=='~'){
                                    const char* pPath= getenv ("HOME");
                                            fileName=string(pPath)+"//"+tokens[1].substr(1);
                                    }
                                    else{
                                        if(tokens[1][0]=='$'){
                                        int pos=tokens[1].find('/');
                                        const char* pPath= getenv (tokens[1].substr(1,pos-1).c_str());

                                                if(!pPath){
                                                    cerr<<"ERROR, unknown variable"<<endl;
                                                return RETVALUE::UNK_LOC_VAR;
                                                }

                                                fileName=string(pPath)+tokens[1].substr(pos);

                                                if(verbose)
                                                    cout<<tokens[1].substr(1,pos-1)<<"="<<pPath<<endl;
                                        }
                                        else
                                            fileName=std::move((tokens[1]));
                                    }

                                    if(verbose)
                                        cout<<"include file: "<<fileName<<endl;


                            fstream file(fileName,ios::in);

                                    if(!file){
                                        cerr<<"ERROR, couldn't open file: "<<tokens[1]<<endl;
                                    return RETVALUE::INC_FILE_ERR;
                                    }


                                    fileConfig.exceptions(ios::badbit);

                                    if(cmdListProcessing(file,prm)!=RETVALUE::OK){
                                        cerr<<"ERROR, file: "<<tokens[1]<<endl;
                                        file.close();
                                    return RETVALUE::INC_FILE_ERR;
                                    }

                                    file.close();
                                    continue;
                            }


                            //--
                            pos=cmdline.find("infile");
                            if(pos!=string::npos){
                            //std::string fileInName(cmdline.substr(7));
                              //      trim(fileInName);

                            vector<string> tokens(split<string> (cmdline," "));
                            std::string fileInName(tokens[1]);

                                    fileInName=encdollar(variables,fileInName);

                            std::string syscmd("ls "+fileInName+">___fileInName.tmp");

                            const int ret=system(syscmd.c_str());

                            std::ifstream fileTmp("___fileInName.tmp");

                                    if(fileTmp){
                                    int numOffiles=0;

                                            while (!fileTmp.eof())
                                                if(fileTmp.get()=='\n')
                                                    numOffiles++;

                                             fileTmp.clear();
                                             fileTmp.seekg(0,ios::beg);


                                             if(numOffiles==0){
                                                cerr<<"ERROR file doesn't exist"<<fileInName<<endl;
                                                throw RETVALUE::IFILE_ERR;
                                             }

                                            fileInName="";

                                            if(numOffiles==1)
                                                fileTmp>>fileInName;
                                            else{
                                            int i=0;
                                            int sel;
                                                cout<<"WARNING: Many file names matching the pattern.\n Select the file: "<<endl;

                                                while(i<numOffiles){
                                                    fileTmp>>fileInName;
                                                    cout<<i++<<"\t"<<fileInName<<endl;
                                                }

                                                cout<<">> ";
                                                cin>>sel;

                                                if(sel>=numOffiles){
                                                    cerr<<"ERROR, out of range"<<endl;
                                                    throw RETVALUE::IFILE_ERR;
                                                }

                                                file.clear();
                                                fileTmp.seekg(0);

                                                i=0;
                                                do{
                                                    fileTmp>>fileInName;
                                                }while(i!=sel);

                                            }

                                       fileTmp.close();
                                       system("rm -f ___fileInName.tmp");

                                    }
                                    else{
                                        cerr<<"ERROR file doesn't exist"<<fileInName<<endl;
                                        throw RETVALUE::IFILE_ERR;
                                    }



                            fstream file(fileInName);

                                    if(!file){
                                        cerr<<"ERROR file doesn't exist: "<<fileInName<<endl;
                                        throw RETVALUE::IFILE_ERR;
                                    }
                                    file.close();


                            ClKeyValues keyvalues;

                                        keyvalues<<("infile")<<fileInName;
                                        cmdlist.emplace_back(keyvalues);
                                        sizeDataBuffer++;

                                continue;
                            }
                            //--
                            if(std::regex_match(cmdline,std::regex("intspace[[:s:]]+(x|y|z)("+strNUMBER+"){2}[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                            keyvalues<<tokens;
                                            cmdlist.emplace_back(keyvalues);

                                continue;
                            }



                            //
                            if(std::regex_match(cmdline,std::regex("marker("+strNUMBER+"){2}[[:s:]]+(x|y|z)[[:s:]]+\\w+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                            keyvalues<<tokens;
                                            cmdlist.emplace_back(keyvalues);

                                continue;
                            }

                            //
                            const string strre("mass[[:s:]]+[aAbB]"+strNUMBER+"[[:s:]]*");
                            if(std::regex_match(cmdline,std::regex(strre))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                            keyvalues<<tokens;
                                            cmdlist.emplace_back(keyvalues);

                                continue;
                            }

                            if(std::regex_match(cmdline,std::regex("merge"))){
                            //vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                            keyvalues<<"merge";
                                            cmdlist.emplace_back(keyvalues);

                                continue;
                            }



                            //--
                            /*
                            if(std::regex_match(cmdline,std::regex("minmax[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                            keyvalues<<tokens;
                                            cmdlist.emplace_back(keyvalues);

                                continue;
                            }*/

                            //--

                            if(std::regex_match(cmdline,std::regex("[bxyz]margin"+strNUMBER+"[[:space:]]*"))){
                            char type=cmdline[0];
                            float value=std::stof(cmdline.substr(1+7));

                                        switch (type){
                                        case 'b' : dataBuffer::xmargin=value;
                                                   dataBuffer::ymargin=value;
                                                   dataBuffer::zmargin=value;break;

                                        case 'x' : dataBuffer::xmargin=value;break;
                                        case 'y' : dataBuffer::ymargin=value;break;
                                        case 'z' : dataBuffer::zmargin=value;break;

                                        default :
                                                    cerr<<" ERROR, invalid command: "<<cmdline<<endl;
                                                    throw RETVALUE::CFILE_ERR;
                                        }

                                        dataBuffer::setbox=false;


                            continue;
                            }
                            //--



                            //--
                            if(std::regex_match(cmdline,std::regex("fixbox("+strNUMBER+"){6}"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            const string *values=&tokens[1];

                                        for(int i=dataBuffer::xmin;i<=dataBuffer::zmax;i++,values++)
                                            dataBuffer::box[i]=std::stof( *values);

                                        dataBuffer::setbox=true;
                            continue;
                            }
                            //--



                            //--
                            pos=cmdline.find("outfile");
                            if(pos!=string::npos){

                            std::string fileOutName=(cmdline.substr(pos+8));

                                    trim(fileOutName);

                                    fileOutName=encdollar(variables,fileOutName);

                            vector<string> fnTokens(split<string>(fileOutName,"."));
                            const string fileFormats("CONFIG xyz poly lmp dmp d2s");

                                    if(fileFormats.find(fnTokens.back())==string::npos){
                                        cerr<<"ERROR, unrecognized format of output file: "<<fnTokens.back()<<endl;
                                        throw RETVALUE::OFILE_ERR;
                                    }




                            fstream file(fileOutName,ios::out);
                                    if(!file){
                                        cerr<<"ERROR, output file couldn't be created: "<<fileOutName<<endl;
                                        throw RETVALUE::OFILE_ERR;
                                    }
                                    file.close();

                             ClKeyValues keyvalues;

                                     keyvalues<<"outfile"<<fileOutName;
                                     cmdlist.emplace_back(keyvalues);

                                    continue;
                            }
                            //--

                            //--
                            if(std::regex_match(cmdline,std::regex("pause[[:s:]]*"))){
                            ClKeyValues keyvalues;
                                          keyvalues<<"pause";
                                          cmdlist.emplace_back(keyvalues);
                            continue;

                            }
                            //--

                            if(std::regex_match(cmdline,std::regex("print([[:space:]]+([xyz]range|numatoms|rematoms|name|file|separator([[:print:]][[:digit:]]{0,2})?))+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                       keyvalues<<tokens;
                                       cmdlist.emplace_back(keyvalues);

                                 continue;
                            }
                            //--

                            //if(std::regex_match(cmdline,std::regex("printf[[:s:]]+\\\"(\\w*[[:s:]]*)*\\\""))) {
                            pos=cmdline.find("printf");
                            if(pos!=string::npos && pos==0){
                            ClKeyValues keyvalues;
                            string token(cmdline.substr(pos+1+6));

                                       keyvalues<<"printf"<<token;
                                       cmdlist.emplace_back(keyvalues);


                                continue;
                            }

                            //--

                            //--
                            if(std::regex_match(cmdline,std::regex("(push|pop)[[:space:]]+\\w+[[:space:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);

                                        if(tokens[0]=="push")
                                            sizeStack++;
                                        else{
                                            popNumber++;
                                        }

                                    continue;
                            }
                            //--
                            if(std::regex_match(cmdline,std::regex("quit[[:s:]]*"))){
                            //vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                           keyvalues<<"quit";
                                           cmdlist.emplace_back(keyvalues);

                                    break;
                            }
                            //--
                            if(regex_match(cmdline,
                               regex("rescale([[:space:]]+[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?){1,3}"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                    keyvalues<<tokens;
                                    cmdlist.emplace_back(keyvalues);

                                    continue;
                            }

                            //--
                            if(regex_match(cmdline,
                               regex("rotate([[:space:]]+[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?){4}[[:space:]]*"))){

                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                       keyvalues<<tokens;
                                       cmdlist.emplace_back(keyvalues);

                                    continue;
                            }
                            //--

                            if(std::regex_match(cmdline,std::regex("randpos[[:s:]]+(loc|tot)+("+strNUMBER+"|[[:s:]]+\\$\\w+){3}[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));                            
                                       replaceVariable(variables,tokens,2,5);

                            ClKeyValues keyvalues;

                                       keyvalues<<tokens;
                                       cmdlist.emplace_back(keyvalues);
                            continue;
                            }

                            //--
                            if(std::regex_match(cmdline,regex("rename[[:s:]]+[[:print:]]+[[:s:]]+[[:print:]]+"))){
                            vector<string> tokens(split<string> (cmdline," "));
                            ClKeyValues keyvalues;

                                       keyvalues<<tokens;
                                       cmdlist.emplace_back(keyvalues);

                            continue;
                            }


                            //--
                             if(regex_match(cmdline,
                                regex ("remove[[:space:]]+plane[[:space:]]+(out|in)("+strNUMBER+"|[[:s:]]+\\$\\w+){3}[[:space:]]+rand("+strNUMBER+"){2}"))){

                                    if(!sizeDataBuffer){
                                        cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                    throw  RETVALUE::CFILE_ERR;
                                    }

                             vector<string> tokens(split<string> (cmdline," "));
                                        replaceVariable(variables,tokens,3,7);

                             ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);
                             continue;
                             }


                            //--
                             if(regex_match(cmdline,
                                regex ("remove[[:space:]]+plane[[:space:]]+(out|in)("+strNUMBER+"|[[:s:]]+\\$\\w+){4}([[:s:]]+[[:print:]]+)?"))){

                                    if(!sizeDataBuffer){
                                        cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                    throw  RETVALUE::CFILE_ERR;
                                    }

                             vector<string> tokens(split<string> (cmdline," "));
                                        replaceVariable(variables,tokens,3,7);

                             ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);
                             continue;
                             }





                             //--remove plane [xyz](max|min)  {value}   (out|in)("+strNUMBER+"|[[:s:]]+\\$\\w+){4}[[:space:]]*
                             if(regex_match(cmdline,regex ("remove[[:space:]]+slice[[:space:]]+[xyz](max|min)"+strNUMBER))){

                                    if(!sizeDataBuffer){
                                        cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                    throw  RETVALUE::CFILE_ERR;
                                    }

                             vector<string> tokens(split<string> (cmdline," "));
                                        replaceVariable(variables,tokens,3,3);

                             ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);
                             continue;
                             }


                              //--
                              if(regex_match(cmdline,
                                 regex ("remove[[:s:]]+slice[[:s:]]+(out|in)("+strNUMBER+"|[[:s:]]+\\$\\w+){5}[[:space:]]*"))){

                                     if(!sizeDataBuffer){
                                         cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                     throw  RETVALUE::CFILE_ERR;
                                     }

                              vector<string> tokens(split<string> (cmdline," "));
                                         replaceVariable(variables,tokens,3,8);// ignore first 3 tokens replace if variable from 3 to 8

                              ClKeyValues keyvalues;

                                         keyvalues<<tokens;
                                         cmdlist.emplace_back(keyvalues);
                              continue;
                              }

                              //--
                               if(regex_match(cmdline,
                                  regex ("remove[[:s:]]+cuboid[[:s:]]+(out|in)("+strNUMBER+"|[[:s:]]+\\$\\w+){2}[[:space:]]*"))){

                                      if(!sizeDataBuffer){
                                          cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                      throw  RETVALUE::CFILE_ERR;
                                      }

                               vector<string> tokens(split<string> (cmdline," "));
                                          replaceVariable(variables,tokens,3,5);// ignore first 3 tokens replace if variable from 3 to 5

                               ClKeyValues keyvalues;

                                          keyvalues<<tokens;
                                          cmdlist.emplace_back(keyvalues);
                               continue;
                               }


                            //--                                                                                                                                                 
                             if(regex_match(cmdline,
                                regex ("remove[[:space:]]+cyl[[:space:]]+xy[[:space:]]+(out|in)"+strNUMBER+"("+strNUMBER+")?"))){

                                    if(!sizeDataBuffer){
                                        cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                        throw  RETVALUE::CFILE_ERR;
                                    }

                             vector<string> tokens(split<string> (cmdline," "));
                             ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);


                                    continue;
                                }

                             //--
                             if(regex_match(cmdline,
                                regex ("remove[[:space:]]+sph[[:space:]]+(out|in)"+strNUMBER))){

                                    if(!sizeDataBuffer){
                                        cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                        throw  RETVALUE::CFILE_ERR;
                                    }

                             vector<string> tokens(split<string> (cmdline," "));
                             ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);


                                    continue;
                                }

                             //--

                             if(regex_match(cmdline,
                                regex ("remove[[:s:]]+random"+strPNUMBER))){

                                    if(!sizeDataBuffer){
                                        cerr<<"ERROR, the "<<cmdline<<" needs at least one input file"<<endl;
                                        throw  RETVALUE::CFILE_ERR;
                                    }

                             vector<string> tokens(split<string> (cmdline," "));
                             ClKeyValues keyvalues;
                             double value=std::stod(tokens[2]);

                                        if(value>100){
                                            cerr<<"ERROR, value should be ranging 0..100"<<endl;
                                        throw RETVALUE:: C_VALRANGE_ERR;
                                        }

                                        keyvalues<<tokens[0]<<tokens[1]<<"void"<<tokens[2];
                                        cmdlist.emplace_back(keyvalues);


                                    continue;
                                }

                             //--
                             //"*  slices {outDirName} {zStart zStep zStop} - save a set of layers in a separate files of selected directory"
                             if(regex_match(cmdline,regex("slices[[:s:]]+[[:print:]]+("+strNUMBER+"){3}"))){
                             vector<string> tokens(split<string> (cmdline," "));
                             ClKeyValues keyvalues;

                                         if(std::stod(tokens[3])<1){
                                             cerr<<"ERROR: step must be the number greater than 1"<<endl;
                                            throw RETVALUE::C_VALRANGE_ERR;
                                         }

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);

                             continue;
                             }




                             //--

                             if(regex_match(cmdline,regex("sort[[:space:]]+[aAdD][[:s:]]*"))){
                              vector<string> tokens(split<string> (cmdline," "));
                              ClKeyValues keyvalues;

                                         keyvalues<<tokens;
                                         cmdlist.emplace_back(keyvalues);

                            continue;
                             }

                             //--

                             if(std::regex_match(cmdline,std::regex("[xyz]space"+strNUMBER+"([[:s:]]+[0-9]+)?"))){
                             vector<string> tokens(split<string> (cmdline," "));
                             ClKeyValues keyvalues;

                                        keyvalues<<tokens;
                                        cmdlist.emplace_back(keyvalues);

                             continue;
                             }


                            //--
                              if(std::regex_search(cmdline,std::regex("system[[:s:]]+"))){
                              string cmdvalues=cmdline.substr(7);
                              ClKeyValues keyvalues;

                                            keyvalues<<"system"<<cmdvalues;
                                            cmdlist.emplace_back(keyvalues);
                                continue;

                              }

                              //--
                              if(regex_match(cmdline,regex("translate("+strNUMBER+"){3}[[:space:]]*"))){
                               vector<string> tokens(split<string> (cmdline," "));
                               ClKeyValues keyvalues;

                                          keyvalues<<tokens;
                                          cmdlist.emplace_back(keyvalues);

                                      continue;
                              }
                              //--

                              if(std::regex_match(cmdline, std::regex("\\$\\w+[[:s:]]*=[[:s:]]*[-._*[:alnum:]]+[[:s:]]*"))){
                              vector<string> tokens(split<string> (cmdline,"="));

                                  variables.emplace(rtrim(tokens[0]),trim(tokens[1]));

                                  continue;
                              }

                            cerr<<"ERROR unknown or invalid command: "<<cmdline<<endl;
                            throw RETVALUE::UNKOPT;

                        }


                        fileConfig.close();
                    }

                    catch (RETVALUE::eRetValue erv){
                        cerr<<"line "<<cmdLineCounter<<endl;

                        fileConfig.close();
                        return erv;

                    }
                    catch (std::exception & e){
                        cerr<<" exception during reading the script file: "<<e.what()<<endl;
                        cerr<<"line: "<<cmdLineCounter<<endl;

                        fileConfig.close();
                        return RETVALUE::CFILE_ERR;
                    }

                    catch(...){
                        cerr<<" ERROR, unknown exception:  "<<endl;
                        cerr<<"line: "<<cmdLineCounter<<endl;

                        fileConfig.close();
                        return RETVALUE::CONV_ERR;
                    }


                    if(verbose){
                        cout<<"::: variables:\n";
                        for (auto& x: variables)
                           std::cout << " [" << x.first << '=' << x.second << ']'<<endl;
                    }



return    (prm.processing) ? RETVALUE::OK : RETVALUE::PROC_ERR;  //END missing after 'if' processing
}
//-----------------------------------------------------------------------------
RETVALUE::eRetValue cmdListCreateProcessing(fstream &file, StInputParametrs &prm)
{

//const int recnest=prm.recnest+1;

fstream &fileConfig=file;
size_t &cmdLineCounter=*prm.cmdLineCounter;
std::string cmdline;
const bool verbose=prm.verbose;
StStructCmd scmd;
int numOfprm=0;
                    try{

                        while(!fileConfig.eof()){

                            std::getline(fileConfig,cmdline);
                            cmdLineCounter++;

                            if(verbose) cout<<"::: " <<cmdline<<endl;

                            ltrim(cmdline);//usun biale znaki na poczatku

                            //--
                            if(cmdline.empty()) continue;
                            //--

                            //--
                            //komentarz
                            if(cmdline[0]=='#' || cmdline[0]=='%') continue;
                            //--

                            //-----------------------------------------------------------
                            if(std::regex_match(cmdline,std::regex("endc[[:s:]]*"))){

                                if(prm.processing){

                                    if(numOfprm!=6){
                                            cerr<<"ERROR, one or more parameters are not set"<<endl;
                                            cerr<<"line: "<<cmdline;
                                    return RETVALUE::C_CREATE_NP_ERR;
                                    }

                                    prm.structPrm->emplace_back(scmd);

                                }

                                return RETVALUE::C_ENDC;
                            }
                            //---




                            //-----------------------------------------------------------

                            if(!prm.processing){

                                if(verbose)
                                    cout<<"---> the previous command will be not processed !"<<endl;

                                continue;
                            }

                            //-----------------------------------------------------------


                            if(std::regex_match(cmdline,std::regex("mode[[:s:]]+(hcp|sg)[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));

                                    scmd.mode=tokens[1];
                                    numOfprm++;
                                    continue;
                            }


                            if(std::regex_match(cmdline,std::regex("hcpa"+strNUMBER+"[[:space:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                                    scmd.hcpa=tokens[1];
                                    numOfprm++;
                                    continue;
                            }


                            if(std::regex_match(cmdline,std::regex("hcpc"+strNUMBER+"[[:space:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));

                                    scmd.hcpc=tokens[1];
                                    numOfprm++;
                                    continue;
                            }


                            if(std::regex_match(cmdline,std::regex("Na"+strNUMBER+"[[:space:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                                    //Na=tokens[1];

                                    scmd.Na=std::stof(tokens[1]);

                                    if(scmd.Na<0){
                                        cerr<<"ERROR, Na<0"<<endl;
                                        throw RETVALUE::WR_FORMAT;
                                    }
                                    numOfprm++;

                                    continue;
                            }


                            if(std::regex_match(cmdline,std::regex("Nc"+strNUMBER+"[[:space:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                                    //Nc=tokens[1];

                                    scmd.Nc=std::stof(tokens[1]);

                                    if(scmd.Nc<0){
                                        cerr<<"ERROR, Nc<0"<<endl;
                                        throw RETVALUE::WR_FORMAT;
                                    }
                                    numOfprm++;

                                    continue;
                            }


                            if(std::regex_match(cmdline,std::regex("atom[[:s:]]+\\w+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                                    //aname=tokens[1];
                                    scmd.aname=tokens[1];
                                    continue;
                            }

                            if(std::regex_match(cmdline,std::regex("layers[[:s:]]+[ABC\\d()]+[[:s:]]*"))){
                            vector<string> tokens(split<string> (cmdline," "));
                                    //seq=tokens[1];
                                    scmd.seq=tokens[1];

                                    numOfprm++;
                                    continue;
                            }


                            cerr<<"ERROR unknown or invalid command: "<<cmdline<<endl;
                            throw RETVALUE::UNKOPT;
                        }//while's end


                    }
                    catch (RETVALUE::eRetValue erv){
                        cerr<<"line "<<cmdLineCounter<<endl;

                        fileConfig.close();
                        return erv;

                    }
                    catch (std::exception e){
                        cerr<<" exception during reading the script file: "<<e.what()<<endl;
                        cerr<<"line: "<<cmdLineCounter<<endl;

                        fileConfig.close();
                        return RETVALUE::CFILE_ERR;
                    }

                    catch(...){
                        cerr<<" ERROR, unknown exception:  "<<endl;
                        cerr<<"line: "<<cmdLineCounter<<endl;

                        fileConfig.close();
                        return RETVALUE::CONV_ERR;
                    }

return   RETVALUE::PROC_ERR;  //END missing after if processing
}


