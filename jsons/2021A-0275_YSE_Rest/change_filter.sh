for i in *json; do
    sed -i '.bak' s'/"filter": "z"/"filter": "r"/' $i
done
rm -f *json.bak
