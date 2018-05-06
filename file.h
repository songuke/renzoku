#ifndef _FILE_H_
#define _FILE_H_

#include <string>
using namespace std;

namespace Renzoku {

class File {
public:
    static void split_path(string path, string &file_name, string &scene_path, string &scene_name) {
        string file = path;
        unsigned int found_dot = file.rfind(".json");
	    if (found_dot == string::npos) {
		
		    if (file[file.length() - 1] == '/') {
			    file = file.substr(0, file.length() - 1);
		    }
		    unsigned int found_slash = file.find_last_of('/');
		    scene_name = file.substr(found_slash + 1);
		    scene_path = file + "/";
	    } else {
		    unsigned int found_slash = file.find_last_of('/');
		    scene_name = file.substr(found_slash + 1, found_dot - found_slash - 1);
		    scene_path = file.substr(0, found_slash + 1);
	    }
        file_name = scene_path + scene_name + ".json";
    }
};

}   // end namespace

#endif