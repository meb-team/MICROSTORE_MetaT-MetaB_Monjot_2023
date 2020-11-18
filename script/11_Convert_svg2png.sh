#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 20/10/2020
#
#!/bin/bash

cd ../result
for result in $(ls)
do
  for file in $(ls $result)
  do
    for image in $(ls $result/$file)
    do
      new=$(echo $image | sed 's/.svg/.png/g')
      if [ $(echo $image | grep ".svg$" | wc -l) == 1 ]
      then
        convert -density 300 $result/$file/$image $result/$file/$new
      fi
    done
  done
done
