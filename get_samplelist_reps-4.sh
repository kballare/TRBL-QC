##########
# Alisa Vershinina
# 12 April 2017
# Goal: get a list of names by a pattern.
# Usage:
#	$1 - pattern
#	$2 - d (dirs) or -f (files)
# kim: change -f 1,2 to -f 1,2,34 to include the reps
##########

for x in $(find . -name "$1*" -type $2)
do
	echo "$(cut -d _ -f 1,2,3,4 <<< $(basename $x))"
done | sort -u | while read y; do echo -n "$y "; done
echo ""
