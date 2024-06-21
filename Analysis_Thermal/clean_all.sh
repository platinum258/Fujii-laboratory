
mkdir Clean
mv * ./Clean
mv ./Clean/Program .
mv ./Clean/*.sh .
mv ./Clean/README .
mv ./Clean/FEM_Data* .
mv ./Clean/Data .

rm -rf ./Clean

cd Program 
sh clean.sh
cd ..

