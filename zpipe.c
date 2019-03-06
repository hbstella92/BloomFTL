/* zpipe.c: example of proper use of zlib's inflate() and deflate()
   Not copyrighted -- provided to the public domain
   Version 1.4  11 December 2005  Mark Adler */

/* Version history:
   1.0  30 Oct 2004  First version
   1.1   8 Nov 2004  Add void casting for unused return values
                     Use switch statement for inflate() return values
   1.2   9 Nov 2004  Add assertions to document zlib guarantees
   1.3   6 Apr 2005  Remove incorrect assertion in inf()
   1.4  11 Dec 2005  Add hack to avoid MSDOS end-of-line conversions
                     Avoid some compiler warnings for input and output buffers
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "zlib.h"

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

#define CHUNK 16384

/*
 * Modified zlib.c
 */
int def(char* source, size_t src_len, char* dest, size_t* dst_len, int level) {
    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK) {
        return ret;
    }

    size_t rem_len = src_len;
    size_t s_ofst=0, d_ofst=0;
    do {
        if(src_len < CHUNK) {
            strm.avail_in = src_len;
            flush = Z_FINISH;
        }
        else {
            if(rem_len < CHUNK) {
                strm.avail_in = CHUNK/2;
                rem_len = 0;
                flush = Z_FINISH;
            }
            else {
                strm.avail_in = CHUNK;
                rem_len -= CHUNK;
                flush = Z_NO_FLUSH;
            }
        }
        memset(in, 0, CHUNK);
        memcpy(in, source+s_ofst, strm.avail_in);
        s_ofst += strm.avail_in;
       
        strm.next_in = in;
        
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = deflate(&strm, flush);
            
            assert(ret != Z_STREAM_ERROR);
            
            have = CHUNK - strm.avail_out;
            memcpy(dest+d_ofst, out, have);
            d_ofst += have;
        } while (strm.avail_out == 0);
        
        assert(strm.avail_in == 0);

    } while (flush != Z_FINISH);
    
    assert(ret == Z_STREAM_END);

    *dst_len = d_ofst;
    (void)deflateEnd(&strm);
    return Z_OK;
}

int inf(char* source, size_t src_len, char* dest, size_t* dst_len) {
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;
    
    size_t rem_len = src_len;
    long int s_ofst=0, d_ofst=0;
    do {
        if(src_len < CHUNK) {
            strm.avail_in = src_len;
            rem_len = 0;
        }
        else {
            strm.avail_in = CHUNK;
            rem_len -= CHUNK;
        }
        memset(in, 0, CHUNK);
        memcpy(in, source+s_ofst, strm.avail_in);
        s_ofst += strm.avail_in;

        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
    
            assert(ret != Z_STREAM_ERROR);

            switch (ret) {
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    (void)inflateEnd(&strm);
                    return ret;
            }

            have = CHUNK - strm.avail_out;
            memcpy(dest+d_ofst, out, have);
            d_ofst += have;
        } while (strm.avail_out == 0);

    } while (ret != Z_STREAM_END);

    *dst_len = d_ofst;
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

void zerr(int ret)
{
    fputs("zpipe: ", stderr);
    switch (ret) {
    case Z_ERRNO:
        if (ferror(stdin))
            fputs("error reading stdin\n", stderr);
        if (ferror(stdout))
            fputs("error writing stdout\n", stderr);
        break;
    case Z_STREAM_ERROR:
        fputs("invalid compression level\n", stderr);
        break;
    case Z_DATA_ERROR:
        fputs("invalid or incomplete deflate data\n", stderr);
        break;
    case Z_MEM_ERROR:
        fputs("out of memory\n", stderr);
        break;
    case Z_VERSION_ERROR:
        fputs("zlib version mismatch!\n", stderr);
    }
}

