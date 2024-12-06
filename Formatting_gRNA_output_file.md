(1) Reorder the columns
	awk '{print $1,'\t',$3,'\t',$2,'\t',$3,'\t',$4,'\t',$5,'\t',$6,'\t',$7,'\t',$8}' <filename>

(2) Sort a column by number of identical occurrences 
	awk  'NR==FNR{a[$1]++;next}{ print a[$1],$0}'<filename>  <filename> |sort -V|sed 's/[0-9]* //'

(3) Get out the repeats
	awk '!a[$1]++'  <filename>

(4) Add a 'g' in front of 19 bases gRNA designed
	awk '{sub("^","g",$3)}; 1' <filename>
	
(5) Only print $1, $3
	awk '{print $1,'\t',$3}' <filename>
	
Remember to replace the column titles to OligoName Sequence and
(6) Add barcode (Esp3I) to gRNAs with 73bp_script
	Cacc g 5' -> 3'
		C 3' -> 5' caaa
