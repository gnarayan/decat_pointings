python json_from_ysepz_target_list.py
git pull
wget --output-file="log.log" "https://docs.google.com/spreadsheets/d/1JPIAXjcy-maVeNMkImRHFnhfoo2ulJQzHkCJOL0AbKs/export?format=csv&gid=0" -O "debass_sample.csv"
python prepare_debass.py $1
python total_time.py jsons/2020B-0053_DEBASS_Brout/EVERYTHING debass250205 > jsons/2020B-0053_DEBASS_Brout/EVERYTHING/EVERYTHING_RAS_DECS.txt
cat jsons/2020B-0053_DEBASS_Brout/EVERYTHING/EVERYTHING_RAS_DECS.txt
git add  jsons/2020B-0053_DEBASS_Brout/EVERYTHING/*
git add  jsons/2020B-0053_DEBASS_Brout/TEMPLATE/*

git commit -am 'update debass jsons'
git push
