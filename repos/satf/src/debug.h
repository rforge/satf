#ifndef __SATF_DEBUG_H__
#define __SATF_DEBUG_H__

#include <stdio.h>
#include <vector>
#include <iostream>

class DebugMethod {
  public:
      DebugMethod(const char* class_name, const char* fn_name, int dbg_level, int at_level);
      ~DebugMethod();
      
      void log(int at_level, const char* format, ...);
      
  public:
      std::string mClassName;
      std::string mFunctionName;
      int mDbgLevel;
      int mAtLevel;
      
      static int mIndentationLevel;
};

int DebugMethod::mIndentationLevel = 0;

inline DebugMethod::DebugMethod(const char* class_name, const char* fn_name, int dbg_level, int at_level) {
    mDbgLevel = dbg_level;
    mAtLevel = at_level;
    mClassName = class_name;
    mFunctionName = fn_name;
    log(mAtLevel, fn_name);
    mIndentationLevel++;
}

inline DebugMethod::~DebugMethod() {
    mIndentationLevel--;
    log(mAtLevel, "--");
}

inline void DebugMethod::log(int at_level, const char * format, ...) 
{
    if(at_level+mAtLevel > mDbgLevel)
      return;
      
    for(int i=0; i < mIndentationLevel; i++)
      putc('\t', stdout);
    
    char buffer[1000];
    va_list args;
    va_start (args, format);
    vsprintf (buffer, format, args);
    printf(" %s", buffer);
    va_end (args);
    printf(" (%s)\n", mClassName.c_str());
}

#if 0
  #define _dbg_class_init                   static const char* m_dbg_class; static const int m_dbg_level
  #define _dbg_class_set(class_name, name, dbg_level) const char* class_name::m_dbg_class = name; const int class_name::m_dbg_level = dbg_level
  #define _dbg_function(fn_name, at_level)  static DebugMethod dbg_function(m_dbg_class, fn_name, m_dbg_level, at_level)
  #define _dbg(param)                       dbg_function.log param
#else
  #define _dbg_class_init                   
  #define _dbg_class_set(class_name, name, dbg_level) 
  #define _dbg_function(fn_name, at_level)  
  #define _dbg(param)                              
#endif

#endif  //__SATF_DEBUG_H__
