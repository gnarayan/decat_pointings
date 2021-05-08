for i in *json; do
    sed -i '.bak' s'/"filter": "z"/"filter": "i"/' $i
done
rm -f *json.bak
