
	#include <stdio.h>
	#include <stdlib.h>
	
	#define USE_ZLIB 1
	
	#ifdef USE_ZLIB
		#include <zlib.h>
	#else
		#define gzopen fopen
		#define gzprintf fprintf
		#define gzFile FILE*
		#define gzclose fclose
		#define gzeof feof
		#define gzclearerr clearerr
		#define gzflush(p1,p2) fflush((p1))
		#define gzrewind rewind
		#define gzread(p1,p2,p3) fread((p2),(p3),1,(p1))
		#define gzwrite(p1,p2,p3) fwrite((p2),(p3),1,(p1))
		#define gzseek fseek
		#define gzgets(p1,p2,p3) fgets((p2),(p3),(p1))
		#define Z_NULL NULL
	#endif
	
	#define IOK 0
	#define IERR 1
	#define MAXFILE 1000

	extern "C" {
		
	int openfile_( const char* fName, const char* mode, int* fID );
	void closefile_( int* fID ); 
	void writestring_( int* fID, const char* str );
	void writedata_( int* fID, void* arr, int* narr );
	void readstring_( int* fID, char* str );
	void readdata_( int* fID, void* arr, int* narr );
	}

	class solverIO {

	public:
		
		solverIO();
		~solverIO();
		
		int openFile (const char* fName, const char *mode);
		int closeFile ();
		int writeString (const char* str);
		int writeData(void* data, int ndata);
		int readString (char* str);
		int readData(void* data, int ndata);
		char* stringStripper ( const char* str );

	private:

		gzFile fPtr_;
		char* mode_;
	};







