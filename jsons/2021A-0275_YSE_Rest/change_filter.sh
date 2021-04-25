for i in *json; do
    sed -i '.bak' s'/"filter": "g"/"filter": "r"/' $i
done
rm -f *json.bak
