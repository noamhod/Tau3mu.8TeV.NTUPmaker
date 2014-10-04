// http://stackoverflow.com/questions/2670816/how-can-i-use-the-compile-time-constant-line-in-a-string

// #ifndef __CINT__
// #include <boost/lexical_cast.hpp>
// #else
// // If needed forward declare the boost entity as 'simply' as possible
// template <typename T> class boost_ptr;
// #endif
// // #include <boost/lexical_cast.hpp>

#ifndef TauLFVCommonTools_logs_H
#define TauLFVCommonTools_logs_H

namespace
{

#include "rawStd.h"

enum DBG
{
	SILENT,
	VISUAL
};

enum MSG
{
	DBG,
	INF,
	WRN,
	ERR,
	FAT
};

static vector<int> msglvl (5,SILENT);

static inline void setMSGlevel(int dbg, int inf, int wrn)
{
	msglvl[DBG] = dbg;
	msglvl[INF] = inf;
	msglvl[WRN] = wrn;
	msglvl[ERR] = VISUAL;
	msglvl[FAT] = VISUAL;
}

string log(const char* file, int line, int type, int level, string message)
{
	ostringstream os;
	
	string stype = "";
	switch(type)
	{
		case DBG: stype = "[DEBUG]   "; break;
		case INF: stype = "[INFO]    "; break;
		case WRN: stype = "[WARNING] "; break;
		case ERR: stype = "[ERROR]   "; break;
		case FAT: stype = "[FATAL]   "; break;
		default:  stype = "[???]     "; break;
	}
	
	os << stype << file << " +" << line << ": " << message;
	
	if(type==ERR) // ERRORS are always shown !
	{
		cout << os.str() << endl;
	}
	else if(type==FAT) // FATALS are always shown !
	{
		os << "\n   exitting now !\n";
		cout << os.str() << endl;
	}
	else if(msglvl[type]==VISUAL)
	{
		if(level==VISUAL) cout << os.str() << endl;
	}
	
	return os.str();
}

#define LOG(x,y,z)  log(__FILE__, __LINE__,  (x), (y),   (z)) // general log
#define _DEBUG(x)   log(__FILE__, __LINE__, DBG, VISUAL, (x)) // for DEBUG VISUAL
#define _INFO(x)    log(__FILE__, __LINE__, INF, VISUAL, (x)) // for INFO VISUAL
#define _WARNING(x) log(__FILE__, __LINE__, WRN, VISUAL, (x)) // for WARNING VISUAL
#define _ERROR(x)   log(__FILE__, __LINE__, ERR, VISUAL, (x)) // for ERROR VISUAL
#define _FATAL(x)   {log(__FILE__, __LINE__, FAT, VISUAL, (x)); exit(-1);} // for FATAL VISUAL


inline bool isnaninf(double x)
{
	if(std::isinf(x)) { _WARNING("value is infinity");     return true; }
	if(std::isnan(x)) { _WARNING("value is not a number"); return true; }
	return false;
}


// static inline double validate_double(string str)
// {
// 	if(str.substr(0,2)=="00")
// 	{
// 		size_t pos = str.find('.');
// 		if(pos==string::npos)
// 		{
// 			LOG(ERR,VISUAL,"cannot validate string=" + str + " (didn't find a decimal period). exit now.");
// 			exit(-1);
// 		}
// 	}
// 
// 	double x;
// 	try
// 		{
// 			x = boost::lexical_cast<double>(str);
// 		}
// 		catch(boost::bad_lexical_cast&)
// 		{
// 			LOG(ERR,VISUAL,"cannot validate string=" + str + " (cannot cast it to double). exit now.");
// 			exit(-1);
// 		}
// 	return x;
// }

static inline bool validate_bool(string str)
{
	if(str=="false" || str=="0") return false;
	else if(str=="true" || str=="1") return true;
	else
	{
		cout << "ERROR: cannot validate string=" << str << " (cannot be read as a boolean). exit now." << endl;
		exit(-1);
	}
	return true;
}

static inline string tostring(short x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return str;
}

static inline string tostring(unsigned short x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return str;
}

static inline string tostring(int x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return str;
}

static inline string tostring(unsigned int x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return str;
}

static inline string tostring(double x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return str;
}

static inline string tostring(float x)
{
	stringstream strm;
	string str;
	strm << x;
	strm >> str;
	return str;
}

static inline string tostring(double x, int prcn)
{
	stringstream strm;
	string str;
	strm << setprecision(prcn) << fixed << x; 
	strm >> str;
	return str;
}

static inline string tostring(float x, int prcn)
{
	stringstream strm;
	string str;
	strm << setprecision(prcn) << fixed << x; 
	strm >> str;
	return str;
}

static inline int toint(string str)
{
	stringstream strm;
	int x;
	strm << str;
	strm >> x;
	return x;
}

static inline double todouble(string str)
{
	stringstream strm;
	double x;
	strm << str;
	strm >> x;
	return x;
}

// static inline double       validate_double(const char* cc) { return validate_double((string)cc); }
// static inline unsigned int validate_uint(string str)       { return (unsigned int)validate_double(str); }
// static inline unsigned int validate_uint(const char* cc)   { return validate_uint((string)cc); }
// static inline int          validate_int(string str)        { return (int)validate_double(str); }
// static inline int          validate_int(const char* cc)    { return validate_int((string)cc); }
// static inline float        validate_float(string str)      { return (float)validate_double(str); }
// static inline float        validate_float(const char* cc)  { return validate_float((string)cc); }
static inline bool         validate_bool(const char* cc)   { return validate_bool((string)cc); }		
static inline int          toint(const char* cc)           { return toint((string)cc); }
static inline double       todouble(const char* cc)        { return todouble((string)cc); }

static inline string _s(unsigned int x)     {return tostring(x);}
static inline string _s(int x)              {return tostring(x);}
static inline string _s(double x)           {return tostring(x);}
static inline string _s(float x)            {return tostring(x);}
static inline string _s(short x)            {return tostring(x);}
static inline string _s(unsigned short x)   {return tostring(x);}
static inline string _s(double x, int prcn) {return tostring(x,prcn);}
static inline string _s(float x, int prcn)  {return tostring(x,prcn);}

}

#endif

