values='7 8 9 10 11 12 13 14 15'
for value in $values
do
    echo ${value}
    wget http://pasha.csail.mit.edu/sets_july24/k${value}/PASHA${value}_100.txt
    wget http://pasha.csail.mit.edu/sets_july24/k${value}/decyc${value}.txt
done
