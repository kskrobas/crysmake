#ifndef CLCMDLIST_H
#define CLCMDLIST_H


#include <vector>
#include <regex>
#include <iostream>
#include <fstream>
#include <string>

//#include <QString>
//#include <QFileInfo>

using namespace std;


#include "iofile.h"
#include "databuffer.h"


namespace  commandListHelp {
extern const std::string clh[];
extern const size_t clhSize;
}

//http://stackoverflow.com/questions/1657883/variable-number-of-arguments-in-c
//#include <iostream>
//#include <string>
//#include <initializer_list>

//template <typename T>
//void func(T t)
//{
//    std::cout << t << std::endl ;
//}

//template<typename T, typename... Args>
//void func(T t, Args... args) // recursive variadic function
//{
//    std::cout << t <<std::endl ;

//    func(args...) ;
//}

//template <class T>
//void func2( std::initializer_list<T> list )
//{
//    for( auto elem : list )
//    {
//        std::cout << elem << std::endl ;
//    }
//}

//int main()
//{
//    std::string
//        str1( "Hello" ),
//        str2( "world" );

//    func(1,2.5,'a',str1);

//    func2( {10, 20, 30, 40 }) ;
//    func2( {str1, str2 } ) ;
//}


#define RE_NUMBER "[[:s:]]+[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"   //signed real value
#define PRE_NUMBER "[[:s:]]+[+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"   //positive real value

//-----------------------------------------------------------------------------
namespace RETVALUE{
enum eRetValue{OK,NARGS_ERR,UNKOPT,CFILE_ERR,CFILE_STOP,C_END,C_END_ERR,C_END_MIS_ERR,
               IFILE_ERR,OFILE_ERR,CONV_ERR,NTODO,CNARGS_ERR,WR_FORMAT,INC_FILE_ERR,UNK_LOC_VAR,
               C_ENDF,ERR_FNAMECOLL,C_ENDF_ERR,C_FUNC_ERR,C_NEST_FUNC_ERR,
               PROC_ERR,C_ELSE,C_ELSE_ERR,STREAM_ERR,
               C_ENDC,C_CREATE_ERR,C_ENDC_ERR,C_CREATE_NP_ERR,C_UNKVAR,C_NONUM,C_VALRANGE_ERR
              };
}
//--------------------------------------------------------------------------
namespace ASTRUCT{

enum mode{UNK,HCP,SG};


}



class ClKeyValues;


class ClKeyValues{
private:

        vector<string> keyvalues;
        string *ptr_frontValue;
        string *ptr_lastValue;

        void fini();
        void fcon();
        void (ClKeyValues::*ptr_f)();

public:

        ClKeyValues(){ keyvalues.reserve(8);}
        ClKeyValues(const char *key){keyvalues.reserve(8);keyvalues.emplace_back(string(key));}

        ClKeyValues & operator << (string &kv);
        ClKeyValues & operator << (const char *);
        ClKeyValues & operator << (const vector<string> &vs ) ;
        //ClKeyValues & operator << (string kv);
        ClKeyValues & operator >> (string &kv);
        ClKeyValues & operator >> (float &f);
        ClKeyValues & operator >> (vector<string> &vs ) ;

       // ClKeyValues operator = (const ClKeyValues &);
       // ClKeyValues operator = (ClKeyValues &&);
        ClKeyValues(const ClKeyValues &);
        ClKeyValues( ClKeyValues &&);

        size_t numOfvalues(){return keyvalues.size()-1;}

};


//-----------------------------------------------------------------------------
// trim from start
static inline std::string &ltrim(std::string &s)
{
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s)
{
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s)
{
        return ltrim(rtrim(s));
}

//-----------------------------------------------------------------------------
template<typename T>
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
}



//-----------------------------------------------------------------------------
struct StAtomsToRem{
bool on=false;
float randRatio;

bool onAB=false;
float randRatioA,randRatioB;
string fileName;

        StAtomsToRem() { }
        StAtomsToRem(const bool on_, const float randRatio_)
        { on=on_; randRatio=randRatio_;}

};
//-----------------------------------------------------------------------------
struct StFunct{
const string name;
const int filePos;
const int cmdPos;




    StFunct():name(""),filePos(-1),cmdPos(-1)
    {

    }


    StFunct(string name__,int filePos__):name(name__),filePos(filePos__),cmdPos(-1)
    {

    }


    StFunct(string name__,int filePos__, int cmdPos__):name(name__),filePos(filePos__),cmdPos(cmdPos__)
    {

    }


};

typedef vector<StFunct> vfunct;


//-----------------------------------------------------------------------------
struct StStructCmd{

const int filePos;
const int cmdPos;

string name;
string mode,aname;
string hcpa,hcpc,seq;
float Na,Nc;

            StStructCmd():filePos(-1),cmdPos(-1)
            {
                aname="X";
            }

//            StStructCmd(string name__,int filePos__, int cmdPos__):name(name__),filePos(filePos__),cmdPos(cmdPos__)
//            {

//            }
};

typedef vector<StStructCmd> vstruct;
//-----------------------------------------------------------------------------
struct StInputParametrs {
size_t *cmdLineCounter, *sizeDataBuffer, *sizeStack, *popNumber;
vector<ClKeyValues> *cmdlist;
vector<string> *cparg;
vfunct *funct;
vstruct *structPrm;
int recnest,fpos;
bool verbose,processing;

        StInputParametrs()
        {
            recnest=-1;
        }


        bool argumentExists(const std::string &arg){

                    for(std::string &str: *cparg)
                        if(str==arg)
                            return true;
        return false;
        }


        bool functionExists(const std::string &fname){
                    fpos=-1;
                    for(StFunct &fun: *funct){
                        fpos++;
                        if(fun.name==fname)
                            return true;
                    }
        return false;
        }


        bool structNameExists(const std::string &structName){
                    for(StStructCmd & crcmd: *structPrm){
                        if(crcmd.name==structName)
                            return true;
                    }
        return false;
        }


};

//-----------------------------------------------------------------------------
RETVALUE::eRetValue cmdListProcessing(fstream &file,StInputParametrs &prm);
RETVALUE::eRetValue cmdListCreateProcessing(fstream &file, StInputParametrs &prm);


typedef std::string String;

//---------------------------------------------------------------------------
struct StExpr{
int start=0,stop=0;
String token;

    StExpr(){ }

    int Length(){return stop-start;}
};
//---------------------------------------------------------------------------
class TExprException: public exception
{
private:
        const char *descr[4]={"unknown","unbalanced brackets","invalid expression","invalid pointer position"};
        unsigned ce=0;
        unsigned pos=0;
public:
        enum ErrorCode{UNKNOWN,UNBRACKETS,EXPRESSION,POSITION};
        TExprException(){ }
        TExprException(const ErrorCode ce_){ce=ce_;}
        TExprException(const ErrorCode ce_,const unsigned pos_){ce=ce_;pos=pos_;}

        virtual const char* what() const throw()
        {
        return descr[ce];
        }

        const	unsigned & getPosition() {return pos;}
} ;




////-----------------------------------------------------------------------------//---------------------------------------------------------------------------
//StExpr buildABC(const unsigned startPos=1)
//{
////const String &str=expression;
////size_t pos=startPos;
////StExpr expr;

////        expr.start=pos;
////        pos++;

////        while( pos<exprLength && str[pos]!=')'){

////            if(IsLetter(str[pos]))
////                expr.token+=str[pos];
////            else{
////                if(IsDigit(str[pos])){
////                StExpr exprUp=buildABCexpand(pos);
////                    pos+=exprUp.Length();
////                    expr.token+=exprUp.token;
////                }
////                else{
////                StExpr exprUp=buildABC(pos);
////                    pos+=exprUp.Length();
////                    expr.token+=exprUp.token;
////                }
////            }
////            pos++;
////        }

//////		if(pos==strLen && str[pos]!=')'  )
//////		throw TExprException(TExprException::EXPRESSION,pos);
////        if(pos>exprLength)
////        throw TExprException(TExprException::EXPRESSION,pos);

////        expr.stop=pos;

////return	expr;
//}
////---------------------------------------------------------------------------
//StExpr TFhcp::buildABCexpand(const unsigned startPos)
//{
//const String &str=expression;
//size_t pos=startPos;
//StExpr expr;
//String digits;

//		expr.start=startPos;

//		while(pos<exprLength && IsDigit(str[pos]) )
//			digits+=str[pos++];

//String letters;

//		if(IsLetter(str[pos]))
//			letters=str[pos];
//		else
//			if(str[pos]=='(') {
//			StExpr exprUp=buildABC(pos);
//					pos+=exprUp.Length();
//					letters=exprUp.token;

//			}
//			else
//			throw TExprException(TExprException::EXPRESSION,pos);

//		expr.stop=pos;

//const size_t intNum=StrToInt(digits);

//			expr.token="";
//			for(size_t i=0;i<intNum;i++)
//				expr.token+=letters;

//return expr;
//}


#endif // CLCMDLIST_H
