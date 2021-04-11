for i in *json; do
    sed -i '.bak' s'/"filter": "i"/"filter": "r"/' $i
done
rm -f *json.bak
