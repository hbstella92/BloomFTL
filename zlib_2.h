#include "zlib.h"
#include "zconf.h"

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

#define CHUNK 16384

int def(char* source, size_t src_len, char* dest, size_t* dst_len, int level);
int inf(char* source, size_t src_len, char* dest, size_t* dst_len);
void zerr(int ret);
