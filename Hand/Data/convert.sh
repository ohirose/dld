for i in `seq -f %03g 1 40`; do
  cat hand$i.txt |tr ',' '\t' > tmp;
  mv tmp hand$i.txt
done;
