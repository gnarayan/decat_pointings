for i in [0-9]*json; do
    sed -i '.bak' s'/"filter": "z"/"filter": "r"/' $i
done
rm -f *json.bak
