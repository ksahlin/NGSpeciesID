for all files in *.fastq; do
  bn = `basename $file .fastq`
  
  NGSpeciesID --ont --consensus --sample_size 500 --m 400 --s 50 --medaka --fastq $file --outfolder ${bn}

done
