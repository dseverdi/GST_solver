
# experimenting and log generating on ANT instances

GST_SOLVER="./GST_solver"
method=""
type=""

# passing arguments from stdin

while getopts ":m:i:" opt; do
   case $opt in
	m) method=$OPTARG
	;;
	i) instance_path=$OPTARG
	;;
   esac
done

# check methods which are available
case $method in
      OPT)    type="GR"
      ;;
      APPROX) type="NR"
      ;;
esac

# detecting input files
for entry in "$instance_path"/*
do
     if [ -f "$entry" ]; then 
	$GST_SOLVER -"$method" -"$type" -f $entry
     fi	
done
	   
