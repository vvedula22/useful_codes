
	#include <stdio.h>
	#include <string.h>
	#include <assert.h>
	#include "zlib.h"

	int inf ( unsigned char* in, int nin, unsigned char* out, int nout );
	
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
		
	extern int infzlibdata_( unsigned char* in, int* m, unsigned char* out, int* n, int* ierr ) {
		int ret;
		int nin = *m;
		int nout = *n;
		
		*ierr = 0;
		ret = inf( in, nin, out, nout);
		if ( ret!=Z_OK ) {
			zerr(ret);
			*ierr = -1;
		}
		return ret;
	}
	
	int inf ( unsigned char* in, int nin, unsigned char* out, int nout ) {
		int ret;
		unsigned have;
		z_stream strm;
		
	    /* allocate inflate state */
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;
		ret = inflateInit(&strm);
		if (ret != Z_OK)
			return ret;

		strm.avail_in = nin;
		strm.next_in = in;

		/* run inflate() on input until output buffer not full */
		do {
			strm.avail_out = nout;
			strm.next_out = out;
			ret = inflate(&strm, Z_NO_FLUSH);
			assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
			switch (ret) {
			case Z_NEED_DICT:
				fprintf(stdout,"zpipe: need dictionary..\n");
				ret = Z_DATA_ERROR;     /* and fall through */
			case Z_DATA_ERROR:
				fprintf(stdout,"zpipe: data error..\n");
			case Z_MEM_ERROR:
				(void)inflateEnd(&strm);
				return ret;
			}
			have = nout - strm.avail_out;
		} while (strm.avail_out == 0);

		/* clean up and return */
		(void)inflateEnd(&strm);
		return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
	}

