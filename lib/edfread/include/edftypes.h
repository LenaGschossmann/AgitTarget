/*  EYELINK APPLICATION CODE      */
/*  (c) 1995 by SR Research Ltd.  */
/*  9 Dec '95 by Dave Stampe      */
/*  PLATFORM-PORTABLE TYPES       */
/*  MAY NEED TO BE MODIFIED       */
/*  FOR Mac, etc.                 */

/*#define FARTYPE _far */  /* for some mixed-model builds */
#define FARTYPE            /* make blank for most DOS, 32-bit, ANSI C */

#ifdef __cplusplus 	/* For C++ definitions */
extern "C" {
#endif

#ifndef _BASETSD_H_ /* windows header */
	#ifndef BYTEDEF
		#define BYTEDEF 1

		typedef unsigned char  byte;
		typedef short          INT16;
		typedef long           INT32;
		typedef unsigned short UINT16;
		typedef unsigned long  UINT32;
	#endif
#endif


#ifndef MICRODEF
	#define MICRODEF 1
	typedef struct {
		INT32  msec;	/* SIGNED for offset computations */
		INT16  usec;
	} MICRO ;
#endif

#ifdef __cplusplus	/* For C++ definitions */
}
#endif
