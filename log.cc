#include "log.h"
#include "time.h"

#include <sstream>
using namespace std;

namespace Renzoku {


LogStream::LogStream() : out(cout), is_enable(true), has_file(false) {
}

LogStream::LogStream(ostream &out) : out(out), is_enable(true), has_file(false) {
}

void LogStream::set_enable(bool is_enable) {
    this->is_enable = is_enable;
}

Log *Log::log = NULL;

Log::Log() : display_time(false) {
	out_stream = new LogStream(cout);
    err_stream = new LogStream(cerr);
    debug_stream = new LogStream();
}

void Log::open(string filename) {
    instance()->_open(filename);
}

void Log::_open(string filename) {
    if (file.is_open())
        file.close();

    file.open(filename.c_str());
    if (! file.is_open()) {
        
        ostringstream oss;
        oss << "Cannot open log file " << filename << ". All logs output to screen.";
        instance()->warn() << oss.str();
        return;
    }

    out_stream->set_file_stream(&file);
    err_stream->set_file_stream(&file);
}

Log::~Log() {
    if (file.is_open())
	    file.close();

    delete debug_stream;
    delete out_stream;
    delete err_stream;
}

Log *Log::instance() {
	if (log == NULL) log = new Log();
	return log;
}

void Log::set_debug(bool is_debug) {
    instance()->debug_stream->set_enable(is_debug);
}

void Log::set_display_time(bool time) {
    display_time = time;
}

LogStream &Log::get_error() {        
    return *err_stream;
}

LogStream &Log::get_warn() {
    return *err_stream;
}

LogStream &Log::get_info() {    
    return *out_stream;
}

LogStream &Log::get_critical() {
    return *err_stream;
}

LogStream &Log::warn() {    
    instance()->get_warn() << instance()->get_time_string() << "WARNING : ";
    return instance()->get_warn();
}

LogStream &Log::info() {
    instance()->get_info() << instance()->get_time_string() << "INFO    : ";
    return instance()->get_info();
}

LogStream &Log::error() {
    instance()->get_error() << instance()->get_time_string() << "ERROR   : ";
    return instance()->get_error();
}

void Log::critical(string s, ExitCode::kExitCode code) {
    instance()->get_critical() << instance()->get_time_string() << "CRITICAL: ";
    instance()->get_critical() << s;
    exit(code);
}

LogStream &Log::debug() {
    *(instance()->debug_stream) << instance()->get_time_string() << "DEBUG   : ";
    return *(instance()->debug_stream);
}

} // end namespace
