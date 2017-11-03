#
# mexopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Note: For the version of system compiler supported with this release,
#       refer to Technical Note 1601 at:
#       http://www.mathworks.com/support/tech-notes/1600/1601.html
#
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building MEX-files via the system ANSI compiler
#
# Copyright 1984-2006 The MathWorks, Inc.
# $Revision: 1.78.4.12 $  $Date: 2006/10/14 11:47:40 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lmwservices -lut"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat"
    fi
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        maci)
#----------------------------------------------------------------------------
            CC='gcc-4.0'
            CFLAGS='-fno-common -no-cpp-precomp -fexceptions'
            CLIBS="$MLIBS -lstdc++"
            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX=g++-4.0
            CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch i386'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            FC='g95'
            FFLAGS="-fexceptions"
            FC_LIBDIR=`$FC -print-file-name=libf95.a 2>&1 | sed -n '1s/\/*libf95\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lf95"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci'
            LDFLAGS="-bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
