qcinvfile=$1


echo "uploading now, please wait..."

curl -F file=@$qcinvfile -F "initial_comment=qcInv" -F channels=C069PPDNQ2H -H "Authorization: Bearer xoxb-1734759934679-6355554410112-0kiSsFIH7zqebY6hkSuPvATO" https://slack.com/api/files.upload > upload.log

