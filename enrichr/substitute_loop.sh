# prep gene list for enrichr input
cut -d, -f3 limmaresults.csv | grep -v NA | grep -v Symbol | head -n5000 > enrichr.txt
# remove quotation marks for a single word
for g in $(cat enrichr.txt); do
	sed 's/\"//g' <<< $g >> enrichr2.txt
done