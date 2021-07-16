!**************************************************
		
		module modGzip
		use, intrinsic :: ISO_C_BINDING
		implicit none
		integer, parameter :: Z_NULL=0
		interface
			integer(kind=c_int) function gzOpen(fName,mode) &
				bind(C, NAME='gzopen')
			use, intrinsic :: ISO_C_BINDING
			character(c_char) fName(*),mode(*)
			end function
		end interface

		interface
			type(c_ptr) function gzGets(fid,str,slen) &
				bind(C, NAME='gzgets')
			use, intrinsic :: ISO_C_BINDING
			type(c_ptr), value :: str
			integer(c_int),value :: fid,slen
			end function
		end interface

		interface
			integer(kind=c_int) function gzRead(fid,buf,blen) &
				bind(C, NAME='gzread')
			use, intrinsic :: ISO_C_BINDING
			integer(c_int), value :: fid,blen
			type(c_ptr), value :: buf
			end function gzRead
		end interface
		
		interface
			integer(kind=c_int) function gzRewind(fid) &
				bind(C, NAME='gzrewind')
			use, intrinsic :: ISO_C_BINDING
			integer(c_int), value :: fid
			end function
		end interface

		interface
			integer(kind=c_int) function gzClose(fid) &
				bind(C, NAME='gzclose')
			use, intrinsic :: ISO_C_BINDING
			integer(c_int), value :: fid
			end function
		end interface
		
		end module modGzip

!**************************************************
