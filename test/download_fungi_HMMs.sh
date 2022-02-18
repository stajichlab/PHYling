if [ ! -d HMM ]; then
wget https://github.com/1KFG/Phylogenomics_HMMs/releases/download/v1.1/fungal_HMMs-1.0.zip
unzip fungal_HMMs-1.0.zip
mv Phylogenomics_HMMs-1.0/HMM .
unlink fungal_HMMs-1.0.zip
rm -rf Phylogenomics_HMMs-1.0
fi

