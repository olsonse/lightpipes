#!/bin/bash
# fullpath.sh
# Given a relative path, create a full path.

fullpath() {
    # check for one arg
    if test "$#" -ne "1" ; then
        echo "Usage : fullpath {relative path}"
        exit -1
    fi

    # check to see if $1 exists
    if test -e "$1" ; then
        local B=`basename $1`
        local P=`dirname $1`
        #  echo BASE:$B  PATH:$P
        pushd $P > /dev/null
    
        if test `pwd` != "/" ; then
            FULLPATH_RETVAL=`pwd`/$B
        else
            FULLPATH_RETVAL=/$B
        fi
        popd > /dev/null
    
        return 0
    else
        return 1
    fi
}

fullpath "$@"
echo $FULLPATH_RETVAL
