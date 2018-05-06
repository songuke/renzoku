#ifndef _LOG_H_
#define _LOG_H_

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#include "vec3.h"
#include "time.h"

namespace Renzoku {

// Custom endline character. 
// This is a work around as currently endl does not work with Log stream.
const char endn = '\n';

struct ExitCode {
    enum kExitCode {
        EXIT_CODE_OK,
        EXIT_CODE_LIBRARY_FAILED_TO_LOAD,
        EXIT_CODE_SHADER_FAILED_TO_COMPILE,
		EXIT_CODE_FILE_FAILED_TO_LOAD,
		EXIT_CODE_SCENE_FILE_ERROR,
    };
};

/**
 * A wrapper of ostream that inserts a prefix string (DEBUG, INFO, ERROR, etc.) 
 * before every output.
 */
class LogStream {
public:
    LogStream();
    LogStream(ostream &out);    
        
    template <typename T>    
    friend LogStream& operator<< (LogStream& stream, const T &v) {
        if (! stream.is_enable) return stream;
        
        stream.out << v;

        if (stream.has_file)
            (*stream.file) << v;

        return stream;
    }
        
    void set_enable(bool is_enable);

    void set_file_stream(ofstream *file) {
        this->file = file;
        this->has_file = true;
    }

protected:
    ostream &out;

    ofstream *file;
    bool has_file;

    bool is_enable;
};


/**
 * A singleton class to collect statistics, errors and warnings from the program.
 */
class Log {
protected:
	Log();
	static Log *log;
	
public: 
    ~Log();
	static Log *instance();
    static void open(string filename);
    
    static void set_debug(bool is_debug);
    void set_display_time(bool time);

    /**
     * Three levels of log:
     * 1. Info: information for coding and testing.
     * 2. Warn: notice needed. Error might occur, but program continues.
     * 3. Error: error encountered. Program still can handle and continue. 
     * 4. Critical: error encountered and program exits immediately.
     */
    static LogStream& info();    
    static LogStream& warn();
    static LogStream& error();    
    static void critical(string s, ExitCode::kExitCode code);
    static LogStream& debug();

protected:
	void _open(string filename);

    LogStream& get_info();    
    LogStream& get_warn();
    LogStream& get_error();    
    LogStream& get_critical();    

protected:
	ofstream file;
	
    LogStream *out_stream;
    LogStream *err_stream;

    LogStream *debug_stream;

protected:
    inline string get_time_string();
    bool display_time;
};

inline string Log::get_time_string() {
    return display_time ? "[" + Time::get_system_time_string() + "] " : "";
}

// TODO: implement a chain to support output to both file and console (require decorator pattern or intercept pattern)

} // end namespace 

#endif
