###filter the BAF_LRR  files by the PFB to make sure the Names in the signal files are in the PFB
BEGIN { FS="\t"; OFS="\t" }
NR==FNR { if ($1!="Name") n[$1]=1; next }
FNR==1 && $1=="Name" { print; next }
($1 in n)
