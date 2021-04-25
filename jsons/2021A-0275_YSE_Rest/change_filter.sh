for i in *json; do
    sed -i '.bak' s'/"filter": "r"/"filter": "z"/' $i
done
rm -f *json.bak
