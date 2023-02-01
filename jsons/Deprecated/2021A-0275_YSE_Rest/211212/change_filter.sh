for i in *json; do
    sed -i '.bak' s'/"filter": "g"/"filter": "i"/' $i
done
rm -f *json.bak
