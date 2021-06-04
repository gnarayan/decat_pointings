for i in *json; do
    sed -i '.bak' s'/"filter": "r"/"filter": "i"/' $i
done
rm -f *json.bak
