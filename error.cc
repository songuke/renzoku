#include <stdio.h>
#include <stdlib.h>

#include "error.h"

namespace Renzoku {
void throw_error(const char *error, const char *file, int line) {
	fprintf(stderr, "%s:%d: error: %s\n", file, line, error);
	exit(1);
}

void throw_warning(const char *warn, const char *file, int line) {
	fprintf(stderr, "%s:%d: warning: %s\n", file, line, warn);
}

void error_alloc(const char *file, int line) {
	throw_error("memory allocation failed.", file, line);
}

void error_file(const char *file, int line, const char *filename) {
	char err[256];
	sprintf(err, "failed to open/read/write/close file %s", filename);
	throw_error(err, file, line);
}
};


