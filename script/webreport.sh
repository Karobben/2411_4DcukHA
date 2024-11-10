mkdir WebReport

for i in $(ls */*/outs/web_summary.html); do 
    NAME=$(echo $i| awk -F"/" '{print $2}') 
    cp $i WebReport/$NAME.html
done 
