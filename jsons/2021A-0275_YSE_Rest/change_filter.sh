for i in *json; do
    sed -i '.bak' s'/"filter": "g"/"filter": "z"/' $i
done
rm -f *json.bak
