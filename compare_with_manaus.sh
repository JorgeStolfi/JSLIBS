#! /bin/bash
# Last edited on 2024-01-10 02:38:01 by stolfi
 
find_all_files_cksum_size_date.sh ./ | sort > .baikal-files.csdf
# rsync -avzu ${SMAN}:programs/c/JSLIBS/.manaus-files.csdf ./

for mac in baikal manaus ; do 
    
  # List files in SAVE directories that are junk:
  cat .${mac}-files.csdf \
    | egrep -e 'SAVE[/]' \
    | egrep -e '(([.]svn.*|JUNK|out|OLD)[/])|([/ ]([.][a-fh-zA-CE-Z]|deps)$)|(([.](a|ho|o)|[~])$)' \
    | sort \
    > .${mac}-SAVE-junk.csdf

  # list files that are not junk:
  cat .${mac}-files.csdf \
    | sort \
    | bool 1-2 - .${mac}-SAVE-junk.csdf \
    > .${mac}-nojunk.csdf
    
  # List files in SAVE directories that are duplicated in other places:
  cat .${mac}-files.csdf | egrep -e '[/ ]SAVE/' | sort -b -k4 > .${mac}-SAVE.csdf
  cat .${mac}-SAVE.csdf | sed -e 's:^.*[/]:/:g' | sort | uniq > .${mac}-SAVE.fnames
  fgrep -F -f .${mac}-SAVE.fnames .${mac}-files.csdf | sort -b -k4 > .${mac}-SAVE-origs.csdf

  cat .${mac}-SAVE.csdf .${mac}-SAVE-origs.csdf \
    | sort | uniq \
    > .${mac}-SAVE-relevant.csdf
    
  cat .${mac}-SAVE-relevant.csdf \
    | gawk \
        ' BEGIN { cprev = ""; sprev = ""; linprev = ""; nd = 0 } 
          //{ if (($1 == cprev) && ($2 == sprev)) { 
                if (nd == 0) { print ""; print linprev; nd = 1; }
                print $0;
                nd++;
              } else { 
                nd = 0;
              } 
              cprev=$1; sprev = $2;
              linprev=$0;
            }
        ' \
    > .${mac}-SAVE-dups.csdf
done

join -j4 -a1 -a2 -e '???' -o 1.1,2.1,1.2,2.2,1.3,2.3,0 .{baikal,manaus}-nojunk.csdf \
  | gawk \
      ' (NF != 7) { printf "** bug NF\n[[%s]]", $0 > "/dev/stderr";  exit(1); } 
        // { 
          if ($1 == "???") { $1 = "??????????"; } 
          if ($2 == "???") { $2 = "??????????"; } 
          if ($3 == "???") { $3 = "-"; } 
          if ($4 == "???") { $4 = "-"; } 
          if ($5 == "???") { $5 = "????-??-??-??????"; } 
          if ($6 == "???") { $6 = "????-??-??-??????"; } 
          printf "%10s %10s  %14s %14s  %s %s  %s\n", $1,$2,$3,$4,$5,$6,$7; next;
        }
      ' \
  > .both-nojunk.ccssddf

