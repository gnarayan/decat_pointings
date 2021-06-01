for i in *json; do
    sed -i '.bak' s'/"filter": "i"/"filter": "g"/' $i
done
rm -f *json.bak
