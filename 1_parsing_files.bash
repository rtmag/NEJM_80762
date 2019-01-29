# Create GSM2136658 patient relationship
more acc.cgi\?acc\=GSE80762 |grep GSM -A1| \
perl -pe 's/.+geoaxema_recenter\)\"\>//g'| \
perl -pe 's/\<\/a\>\<\/td\>//g'| \
perl -pe 's/.+\"top\"\>//g'| \
perl -pe 's/\<\/td\>\n//g' | \
perl -pe 's/\n/\t/g' | \
perl -pe 's/\-\-/\n/g' | \
perl -pe 's/\tGSM/GSM/g' | \
perl -pe 's/Patient\s/Patient\_/g'| \
perl -pe 's/\,\sDay.+/\_Untreated/g'| \
perl -pe 's/\,\sCycle.+/\_Treated/g' > gsm_patient_table.txt

# Unzip IDATs
gunzip *idat.gz
ls -1 *idat|perl -pe 's/\_/\t/g'|cut -f1,2,3|sort|uniq > ../idat_table.txt

## R Code for matching GSM with proper labels ##
gsm = read.table("gsm_patient_table.txt")
idat = read.table("idat_table.txt")
