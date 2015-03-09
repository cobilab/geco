# GeCo #
<p align="center"><img src="/logo.png" 
alt="EAGLE" width="350" height="250" border="0" /></p>
Compress and analyze genomic sequences. GeCo can afford individual compression and referential compression. As a compression tool, GeCo is able to provide additional compression gains over several top specific tools, while as an analysis tool, GeCo is hable to determine absolute measures, namely for many distance computations, and local measures, such as the information content contained in each element, providing a way to quantify and locate specific genomic events.

## INSTALLATION ##

Simply run the following instructions at a Linux terminal:

<pre>
wget https://github.com/pratas/geco/archive/master.zip
unzip master.zip
cd geco-master
cmake . 
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
  -rm &#60ctx&#62:&#60den&#62:&#60ir&#62   reference context model       
  -rm &#60ctx&#62:&#60den&#62:&#60ir&#62   reference context model
  ...
  -tm &#60ctx&#62:&#60den&#62:&#60ir&#62   target context model  
  -tm &#60ctx&#62:&#60den&#62:&#60ir&#62   target context model
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

                                                    

