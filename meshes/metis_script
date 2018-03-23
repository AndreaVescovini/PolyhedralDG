#!/bin/bash

# Author: Andrea Vescovini
#
# Two parameters are required
# $1 is the mesh file to be partitioned
# $2 is the number of partition I want metis to do

# 1296p is generated with 400
# 3072p is generated with 979

# extract E2V connectivity matrix needed by metis
sed -n '/Tetrahedra/,/Triangles/p' $1 | sed '1d;$d' | sed '$d' | sed 's/ 0//' > metis_input

# launch metis
mpmetis -ncommon=3 -minconn -contig metis_input $2

# check metis output
echo "Checking missing polyhedra..."
end=$(($2-1))
ii=0
while [ $ii -le $end ]; do
  grep -q -w $ii metis_input.epart.$2
  # if there isnt't the polyhedron ii I have to decrease the total number and to
  # substitute the previous last polyhedron with the number ii
  if [ $? -eq 1 ]; then
    awk -i inplace -v miss=$ii '{if($1>miss) $1--; print}' metis_input.epart.$2
    let end-=1
  fi
  let ii+=1
done
echo "Polyhedra after check: "$((end+1))
echo "******************************************************************************"

# create the file with the final mesh
final_name="$(echo $1 | cut -f 1 -d '.')_metis.mesh"
sed '/End/d' $1 > $final_name

# add metis output
echo "Polyhedra" >> $final_name
echo $((end+1)) >> $final_name
# I sum 1 to the metis output i order to have polyhedra numbered from 1 to N
awk '{print $1+1}' metis_input.epart.$2 >> $final_name
echo "" >> $final_name
echo End >> $final_name

# delete metis output files
rm metis_input.npart.$2
rm metis_input.epart.$2
rm metis_input
