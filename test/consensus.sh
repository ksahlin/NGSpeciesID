for file in *.fastq; do
  bn=`basename $file .fastq`
  NGSpeciesID --ont --consensus --sample_size 500 --m 800 --s 100 --medaka --primer_file primers.txt --fastq $file --outfolder ${bn}
done
