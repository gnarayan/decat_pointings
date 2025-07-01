#export OBSPLAN_ROOTDIR="/Users/arest/observing/decat-obsplan-generator"
if [[ $HOSTNAME =~ arminmac* ]] || [[  $HOSTNAME =~ lswlan* ]]; then
	export OBSPLAN_ROOTDIR="/Users/arest/observing/decat_pointings"
	alias cdjson='cd $OBSPLAN_ROOTDIR/json_files'
	alias cdobsplan='cd $OBSPLAN_ROOTDIR'
elif [[ $HOSTNAME =~ plhstproc* ]]; then
	export OBSPLAN_ROOTDIR="/astro/armin/observing/decat_pointings"
	alias cdjson='cd $OBSPLAN_ROOTDIR/json_files'
	alias cdobsplan='cd $OBSPLAN_ROOTDIR'
else
   echo "Hostname $HOSTNAME is not defined yet in the sourceme file!"
   return 1;
fi


