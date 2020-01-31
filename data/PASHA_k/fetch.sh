values='20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200'
for value in $values
do
    echo ${value}
    wget http://pasha.csail.mit.edu/sets_july24/k13/PASHA13_${value}.txt
done
wget http://pasha.csail.mit.edu/sets_july24/k13/decyc13.txt
