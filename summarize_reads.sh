#!/bin/sh

/bin/ls *.fastq *.fq > fqlist.txt

echo "Computing Read lens"
rm -f fqlist.lens.txt
for i in `cat fqlist.txt`
do
  if [ ! -r $i.lens ]
  then
    echo $i
    echo $i >> fqlist.lens.txt
  fi
done

if [ -s fqlist.lens.txt ]
then
  cat fqlist.lens.txt | parallel -t 'fastq_to_fasta_fast {} | getlengths > {}.lens'
  rm -f allreads.csv allreads.stats
fi

echo "Checking Read len stats"
rm -f fqlist.stats.txt
for i in `cat fqlist.txt`
do
  if [ ! -r $i.lens.stats ]
  then
    echo $i
    echo $i >> fqlist.lens.txt
  fi
done

if [ -s fqlist.lens.txt ]
then
  cat fqlist.lens.txt | parallel -t 'stats -f 2 {}.lens -big 10000,20000,30000,40000,50000 > {}.lens.stats'
fi

if [ ! -r allreads.stats ]
then
  echo "Summarizing all reads"
  awk '{print $2}' *.lens | stats -big 10000,20000,30000,40000,50000 > allreads.stats
fi

if [ ! -r allreads.csv ]
then
  echo -e "run\t#reads\tmax\tmean\tn50\tMb_total\tMb>10kb\tMb>20kbp\tMb>30kbp\tMb>40kb\tMb>50kb\n" > allreads.csv
  head *.stats | awk '{printf("%s", $0); if(NR%3==0){print}}' | tr '=' ' ' | tr -d '[' | tr -d ']' | \
        awk '{print $2,$6,$8,$9,$15,int($13/1000000),int($23/1000000),int($27/1000000),int($31/1000000),int($35/1000000),int($39/1000000)}' | tr ' ' '\t' >> allreads.csv
fi

column -t allreads.csv
