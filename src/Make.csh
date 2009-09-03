# Make.csh = update Makefile.lib and Makefile.list
# use current list of *.cpp and *.h files in src dir

if ($1 == "Makefile.lib") then

  set list1 = `ls [abcde]*.cpp`
  set list2 = `ls [fghij]*.cpp`
  set list3 = `ls [klmno]*.cpp | sed s/^main\.cpp//`
  set list4 = `ls [pqrst]*.cpp`
  set list5 = `ls [uvwxyz]*.cpp`
  sed -i -e "s/SRC =\t.*/SRC =\t$list1/" Makefile.lib
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list2/" Makefile.lib
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list3/" Makefile.lib
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list4/" Makefile.lib
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list5/" Makefile.lib

  set list1 = `ls [abcde]*.h`
  set list2 = `ls [fghij]*.h`
  set list3 = `ls [klmno]*.h`
  set list4 = `ls [pqrst]*.h`
  set list5 = `ls [uvwxyz]*.h`
  sed -i -e "s/INC =\t.*/INC =\t$list1/" Makefile.lib
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list2/" Makefile.lib
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list3/" Makefile.lib
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list4/" Makefile.lib
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list5/" Makefile.lib

else if ($1 == "Makefile.list") then

  set list1 = `ls [abcde]*.cpp`
  set list2 = `ls [fghij]*.cpp`
  set list3 = `ls [klmno]*.cpp`
  set list4 = `ls [pqrst]*.cpp`
  set list5 = `ls [uvwxyz]*.cpp`
  sed -i -e "s/SRC =\t.*/SRC =\t$list1/" Makefile.list
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list2/" Makefile.list
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list3/" Makefile.list
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list4/" Makefile.list
  sed -i -e "s/SRC =\t\(.*\)/SRC =\t\1 $list5/" Makefile.list

  set list1 = `ls [abcde]*.h`
  set list2 = `ls [fghij]*.h`
  set list3 = `ls [klmno]*.h`
  set list4 = `ls [pqrst]*.h`
  set list5 = `ls [uvwxyz]*.h`
  sed -i -e "s/INC =\t.*/INC =\t$list1/" Makefile.list
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list2/" Makefile.list
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list3/" Makefile.list
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list4/" Makefile.list
  sed -i -e "s/INC =\t\(.*\)/INC =\t\1 $list5/" Makefile.list

else if ($1 == "style") then

  set list = `grep -l SolveClass solve_*.h`
  if (-e style_solve.tmp) then
    rm style_solve.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_solve.tmp
  end
  if (! -e style_solve.h) then
     mv style_solve.tmp style_solve.h
  else if (`diff style_solve.h style_solve.tmp` != "") then
     mv style_solve.tmp style_solve.h
  else
     rm style_solve.tmp
  endif

endif
