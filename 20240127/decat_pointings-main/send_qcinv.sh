qcinvfile=$1
token=`cat token.txt`

echo "uploading now, please wait..."

curl -F file=@$qcinvfile -F "initial_comment=qcInv" -F channels=C069PPDNQ2H -H "Authorization: Bearer $token" https://slack.com/api/files.upload > upload.log

