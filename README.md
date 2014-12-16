# GeCo #

GeCo: a flexible genomic [ACGT] compressor 

=> conditional compressor;
=> conditional exclusive compressor (transitive information)
=> self-information compressor (most common compression)]

## INSTALLATION ##

Simply run the following instructions at a Linux terminal:

<pre>
wget https://github.com/pratas/geco/archive/master.zip
unzip master.zip
cd geco-master
make
</pre>

## EXECUTION

### Run GeCo

Run GeCo using:

<pre>
./GeCo -tm 4:1:0 -tm 8:1:0 -tm 14:1:0 -tm 1:18:1 file.seq
</pre>

## PARAMETERS

To see the possible options type
<pre>
./GeCo
</pre>
or
<pre>
./GeCo -h
</pre>

These will print the following options:
<pre>
<p>
Usage: GeCo &#60OPTIONS&#62 ... -r &#60FILE&#62  [FILE]:&#60...&#62

  -v                     verbose mode             
  -f                     force (be sure!)             

  -rm &#60ctx&#62:&#60den&#62:&#60ir&#62  reference context model       
  -rm &#60ctx&#62:&#60den&#62:&#60ir&#62  reference context model
  ...
  -tm &#60ctx&#62:&#60den&#62:&#60ir&#62  target context model  
  -tm &#60ctx&#62:&#60den&#62:&#60ir&#62  target context model
  ...

  -g  &#60gamma&#62            gamma factor

  -r  &#60rFile&#62            reference file

[tFile1]:&#60tFile2&#62:&#60...&#62  target file(s)</p>
</pre>
## CITATION ##

On using this software/method please cite:

Armando J. Pinho, Diogo Pratas, and Paulo J.S.G. Ferreira. "Bacteria DNA sequence compression using a mixture of finite-context models." Statistical Signal Processing Workshop (SSP), 2011 IEEE, pp.125,128, 28-30 June 2011.

DOI: 10.1109/SSP.2011.5967637

## ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/GeCo/issues).

## LICENSE ##

GPL v2.

For more information:
<pre>http://www.gnu.org/licenses/gpl-2.0.html</pre>

                                                    

