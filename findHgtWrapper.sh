#!/bin/bash


####################################################
d=`dirname $0`


####################################################
indir_str=''
blast_file=''
ref_file=''
is_alpha=''
is_fungi=''
cpu=8


####################################################
errorMsg(){
	echo "$1 has to be provided! Exiting ......" >&2
	exit 1
}


####################################################
while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			indir_str="$indir_str --indir $2"
			shift
			;;
		-b)
			blast_file=$2
			shift
			;;
		--ref)
			ref_file="$ref_file --ref $2"
			shift
			;;
		--res_outdir)
			res_outdir=$2
			shift
			;;
		--alpha)
			is_alpha="--alpha"
			;;
		--fungi)
			is_fungi="--fungi"
			;;
		--cpu)
			cpu=$2
			shift
			;;
	esac
	shift
done


####################################################
[ -z "$indir_str" ] && errorMsg "--indir"
[ -z $blast_file ] && errorMsg "-b"
[ -z "$ref_file" ] && errorMsg "--ref"
[ -z $res_outdir ] && errorMsg "--res_outdir"


####################################################
b=`basename $blast_file`
c=${b%%.*}

res=$res_outdir/$c.res

echo "getProt2Taxon ......"
ruby $d/getProt2Taxon.rb --cpu $cpu $indir_str -b $blast_file > $res;

echo "parseTaxonInBlast ......"
ruby $d/parseTaxonInBlast.rb $ref_file -i $res --cpu $cpu $is_alpha $is_fungi > $res_outdir/taxa.$c.res;


if [ $? == 0 ]; then
	echo "Done!"
else
	echo "Error!"
fi


