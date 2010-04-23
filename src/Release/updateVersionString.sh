#!/bin/bash

export REVISION=`svn info ../flameSys.cpp | grep Revision | perl -p -e 's/^.*?(\d+)$/$1/'`
perl -pi -e "s/std::string REVISION = \".*\";/std::string REVISION = \"$REVISION\";/" ../strainedFlame.cpp
export DATE=`date +%F\ %T`
perl -pi -e "s/std::string BUILDDATE = \".*\";/std::string BUILDDATE = \"$DATE\";/" ../strainedFlame.cpp
