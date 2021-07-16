
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "myzlib.h"

	solverIO** fpSolverIO=NULL;
	
	solverIO::solverIO () {
		mode_ = NULL;
	}

	solverIO::~solverIO () {
		if ( mode_ != NULL ) delete[] mode_;
	}

	int solverIO::openFile (const char *fName, const char *mode) {

		fPtr_ = Z_NULL ;
		mode_ = stringStripper(mode);
		
		if ( !strcmp ( mode, "read" ) )
			fPtr_ = gzopen ( fName, "rb" );
		else if ( !strcmp ( mode_, "write" ) )
			fPtr_ = gzopen ( fName, "wb" );
		else if ( !strcmp ( mode, "append" ) )
			fPtr_ = gzopen ( fName, "ab" );
		else {
			fprintf ( stdout, "Error: invalid mode '%s' \n", mode);
			return IERR;
		}

		if ( fPtr_ == Z_NULL ) {
			fprintf(stdout,"Error opening file..\n");
			return IERR;
		}
		else
			return IOK;
	}

	int solverIO::closeFile() {
		if ( strcmp(mode_, "write") || strcmp(mode_,"append") ) {
		  gzflush(fPtr_, Z_FULL_FLUSH);
		}
		gzclose(fPtr_);
		return IOK;
	}

	int solverIO::writeString( const char* str ) {
		gzprintf(fPtr_,"%s",str);
		return IOK;
	}

	int solverIO::readString( char* str ) {
		while(gzgets(fPtr_,str,1024) != Z_NULL) {
			return IOK;
		}
	}

	int solverIO::writeData( void* data, int ndata ) {
		size_t typeSize = sizeof(double);
		fprintf(stdout,"typeSize: %d\n",typeSize);
        gzwrite(fPtr_,static_cast< char* >(data),typeSize*ndata);
		return IOK;
	}

	int solverIO::readData( void* data, int ndata ) {
		size_t typeSize = sizeof(double);
		fprintf(stdout,"typeSize: %d\n",typeSize);
        gzread(fPtr_,data,typeSize*ndata);
		return IOK;
	}

	int openfile_( const char* fName, const char* mode, int* fID ) {

		if (fpSolverIO == NULL) {
			fpSolverIO = new solverIO* [MAXFILE];
			for ( int i = 0; i < MAXFILE; i++ ) {
				fpSolverIO[i] = NULL;
			}
		}
		
		for ( int i = 1; i < MAXFILE; i++ ) {
			if ( fpSolverIO[i] == NULL ) {
				fpSolverIO[i] = new solverIO();
				if ( fpSolverIO[i]->openFile(fName, mode) == IERR) {
					delete fpSolverIO[i];
					fpSolverIO[i] = NULL;
					*fID = 0;
					return IERR;
				}
				*fID = i;
				fprintf(stdout,"File %s opened in %s mode with id %d\n", \
					fName, mode, *fID);
				return IOK;
			}
		}

		fprintf(stderr,"Error: Could not open file\n");
		fprintf(stderr,"       Exceeded file limit %d \n", MAXFILE);
		exit(-1);
		return IERR;

	}

	void closefile_( int* fID ) {
		fpSolverIO[*fID]->closeFile();
		delete fpSolverIO[*fID];
		fpSolverIO[*fID] = NULL; 
		fprintf(stdout,"File with id %d is closed\n",*fID);
		return;
	}

	void writestring_( int* fID, const char* str ) {
		fpSolverIO[*fID] -> writeString(str);
	}

	void writedata_( int* fID, void* arr, int* narr ) {
		int num = *narr;
		fpSolverIO[*fID] -> writeData(arr, num);
	}

	void readstring_( int* fID, char* str ) {
		fpSolverIO[*fID] -> readString(str);
	}

	void readdata_( int* fID, void* arr, int* narr ) {
		int num = *narr;
		fpSolverIO[*fID] -> readData(arr, num);
	}

	char* solverIO::stringStripper ( const char* str ) {
		char* fName;
		int slen = strcspn( str, " " );
		fName = new char [ slen+1 ];
		strncpy( fName, str , slen );
		fName [ slen ] = '\0';
		return fName;
	}

