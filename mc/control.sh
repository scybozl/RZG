for file in $(find ~/mc/ -maxdepth 1 -type d -name "00*") #or you can use the awk one here
do

echo $file
find $file -name "*.yoda" | wc -l

done
