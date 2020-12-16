#! /bin/bash 
# Last edited on 2020-12-08 02:32:34 by jstolfi

svn status | egrep -e '^[AM]' > .libs-git
for d in `cat .libs-git` ; do 
  ???