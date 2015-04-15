# GeCo #
<p align="center"><img src="/logo.png" 
alt="GeCo" width="350" height="260" border="0" /></p>
Compress and analyze genomic sequences. As a compression tool, GeCo is able to provide additional compression gains over several top specific tools, while as an analysis tool, GeCo is able to determine absolute measures, namely for many distance computations, and local measures, such as the information content contained in each element, providing a way to quantify and locate specific genomic events. GeCo can afford individual compression and referential compression.

## INSTALLATION ##

Cmake is needed for installation (http://www.cmake.org/). You can download it directly from http://www.cmake.org/cmake/resources/software.html or use an appropriate packet manager. In the following instructions we show the procedure to install, compile and run GeCo:

### STEP 1

Download, install and resolve conflicts.

#### Linux 
<pre>
sudo apt-get install cmake
wget https://github.com/pratas/geco/archive/master.zip
unzip master.zip
cd geco-master
cmake .
make
</pre>

Alternatively, you can install (without cmake and only for linux) using

<pre>
wget https://github.com/pratas/geco/archive/master.zip
unzip master.zip
cd geco-master
mv Makefile.linux Makefile
make
</pre>

#### OS X
Install brew:
<pre>
ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"
</pre>
only if you do not have it. After type:
<pre>
brew install cmake
brew install wget
brew install gcc48
wget https://github.com/pratas/geco/archive/master.zip
unzip master.zip
cd geco-master
cmake .
make
</pre>
With some versions you might need to create a link to cc or gcc (after the *brew install gcc48* command), namely
<pre>
sudo mv /usr/bin/gcc /usr/bin/gcc-old   # gcc backup
sudo mv /usr/bin/cc /usr/bin/cc-old     # cc backup
sudo ln -s /usr/bin/gcc-4.8 /usr/bin/gcc
sudo ln -s /usr/bin/gcc-4.8 /usr/bin/cc
</pre>
In some versions, the gcc48 is installed over /usr/local/bin, therefore you might need to substitute the last two commands by the following two:
<pre>
sudo ln -s /usr/local/bin/gcc-4.8 /usr/bin/gcc
sudo ln -s /usr/local/bin/gcc-4.8 /usr/bin/cc
</pre>

#### Windows

In windows use cygwin (https://www.cygwin.com/) and make sure that it is included in the installation: cmake, make, zcat, unzip, wget, tr, grep (and any dependencies). If you install the complete cygwin packet then all these will be installed. After, all steps will be the same as in Linux.

## EXECUTION

### Run GeCo

Run GeCo using (lazy) level 5:

<pre>
./GeCo -l 5 File.seq
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
Usage: GeCo [OPTION]... -r [FILE]  [FILE]:[...]                        
Compress and analyze a genomic sequence (by default, compress).        
                                                                       
Non-mandatory arguments:                                               
                                                                       
  -h                     give this help,                                
  -x                     show several running examples,                 
  -s                     show GeCo compression levels,                  
  -v                     verbose mode (more information),               
  -V                     display version number,                        
  -f                     force overwrite of output,                     
  -l &#60level&#62             level of compression [1;9] (lazy -tm setup),   
  -g &#60gamma&#62             mixture decayment forgetting factor. It is     
                         a real value in the interval [0;1),            
  -c &#60cache&#62             maximum collisions for hash cache. Memory      
                         values are higly dependent of the parameter    
                         specification,                                 
  -e                     it creates a file with the extension ".iae"  
                         with the respective information content. If    
                         the file is FASTA or FASTQ it will only use    
                         the "ACGT" (genomic) data,                   

  -r &#60FILE&#62              reference file ("-rm" are loaded here),      

  -rm &#60c&#62:&#60d&#62:&#60i&#62:&#60m/e&#62  reference context model (ex:-rm 13:100:0:0/0), 
  -rm &#60c&#62:&#60d&#62:&#60i&#62:&#60m/e&#62  reference context model (ex:-rm 18:1000:0:1/1000),
  ...                                                                  
  -tm &#60c&#62:&#60d&#62:&#60i&#62:&#60m/e&#62  target context model (ex:-tm 4:1:0:0/0),      
  -tm &#60c&#62:&#60d&#62:&#60i&#62:&#60m/e&#62  target context model (ex:-tm 18:20:1:2/10),   
  ...                                                                  
                         target and reference templates use &#60c&#62 for    
                         context-order size, &#60d&#62 for alpha (1/&#60d&#62),    
                         &#60i&#62 (0 or 1) to set the usage of inverted     
                         repeats (1 to use) and &#60m&#62 to the maximum     
                         allowed mutation on the context without       
                         being discarded (usefull in deep contexts),   
                         under the estimator &#60e&#62,                      

Mandatory arguments:                                                   
                                                                       
  &#60FILE&#62                 file to compress (last argument). For more     
                         files use splitting ":" characters.          
                                                                       
Report bugs to &#60{pratas,ap,pjf}@ua.pt&#62. 
</pre>

To get several running examples type:
<pre>
./GeCo -x
</pre>
This will print the following:
<pre>
                                                                       
GeCo running examples:                                                 
                                                                       
Considerations: the decompression is symmetric, therefore the same     
resources, namely time and memory will be used as in the compression.  
The memory used, after creating the models, will be constant, even in  
deeper context models (cache-hash context model).                      
                                                                       
[A]=> Compressing sequences C(X) or C(X,Y):                            
                                                                       
1) Compression of a human genome (using 5.8 GB RAM memory):            
   ./GeCo -tm 6:1:0:0/0 -tm 13:20:1:0/0 -tm 19:50:1:2/10 -c 35 -g 0.8 HS
                                                                       
2) Compression of a human genome (using 3.8 GB of RAM memory):         
   ./GeCo -tm 6:1:0:0/0 -tm 13:10:1:0/0 -tm 18:50:1:3/10 -c 20 -g 0.8 HS
                                                                       
3) Compression of a human genome (using 0.5 GB of RAM memory):         
   ./GeCo -tm 6:1:0:0/0 -tm 13:10:1:0/0 -g 0.8 HS                      
                                                                       
   Decompression for A1, A2 and A3: ./GeDe HS.co                       
   The decompressed file will be HS.de                                 
                                                                       
4) Compression of a human chromosome Y (repetitive nature):            
   ./GeCo -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:1:0/0 -tm 8:1:0:0/0      
   -tm 11:10:1:0/0 -tm 14:10:0:1/10 -tm 14:50:1:0/0 -tm 18:30:1:6/10   
   -c 10 -g 0.88 CY.fasta                                              
   Decompression for A4: ./GeDe CY.fasta.co                            
   The decompressed file will be CY.fasta.de                           
                                                                       
5) Highly-redundant genomic sequence (full ACGT from fastq)            
   ./GeCo -tm 4:1:0:0/0 -tm 11:1:0:0/0 -tm 14:20:0:0/0 -tm 20:100:0:1/10
   -c 40 -g 0.8 SRR957627.fastq                                        
   Decompression for A5: ./GeDe SRR957627.fastq.co                     
   The decompressed file will be SRR957627.fastq.de                    
                                                                       
                                                                       
[B]=> Conditional (referential) exclusive compression C(X||Y):         
                                                                       
1) Compression of the gorilla (GG8) chromosome 8 given exclusively     
   information from chimpanzee (PT8):                                  
   ./GeCo -rm 4:1:0:0/0 -rm 20:1000:1:1/100 -c 20 -r PT8 GG8           
   Decompression for B1: ./GeDe -r PT8 GG8.co                          
   The decompressed file will be GG8.de                                
                                                                       
2) Compression of the same file (for identity studies):                
   ./GeCo -rm 20:1000:0:0/0 -c 30 -r File1.txt File1.txt               
   Decompression for B2: ./GeDe -r File1.txt File1.txt.co              
   The decompressed file will be File1.txt.de                          
                                                                       
3) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  
   chromsomes 5 given exclusively the gorilla (GG17) chromosome 17 as  
   reference:                                                          
   ./GeCo -rm 20:1000:1:1/100 -c 20 -r GG17 HS5:PT5:PA5                
   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           
   The decompressed files will be HS5.de, PT5.de and PA5.de            
                                                                       
                                                                       
[C]=> Conditional compression C(X|Y) [use reference and target]:       
                                                                       
1) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  
   chromsomes 5 given the gorilla (GG17) chromosome 17 as reference:   
   -rm 12:100:1:0/0 -rm 20:1000:1:1/100 -tm 4:1:0:0/0 -tm 14:20:1:1/10 
   -c 20 -g 0.85 -r GG17 HS5:PT5:PA5                                   
   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           
   The decompressed files will be HS5.de, PT5.de and PA5.de 
</pre>

If you are not interested in setting the template for each model, then use the levels mode. To see the possible levels type:
<pre>
./GeCo -s
</pre>
This will print the following:
<pre>
Level 1: -tm 13:1:0:0/0 -g 0.85 
Level 2: -tm 4:1:0:0/0 -tm 11:10:0:0/0 -tm 13:20:1:0/0 -tm 19:20:1:2/10 -c 30 -g 0.85 
Level 3: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 8:1:0:0/0 -tm 12:10:0:0/0 -tm 13:20:1:1/10 -tm 18:100:1:6/10 -c 9 -g 0.90 
Level 4: -tm 4:1:0:0/0 -tm 6:1:1:0/0 -tm 13:20:1:0/0 -tm 18:20:1:3/10 -c 20 -g 0.9 
Level 5: -tm 4:1:0:0/0 -tm 11:1:0:0/0 -tm 13:20:1:0/0 -tm 19:20:1:3/10 -c 30 -g 0.85 
Level 6: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 9:1:0:0/0 -tm 11:1:0:0/0 -tm 13:10:1:0/0 -tm 14:20:1:0/0 -tm 19:20:1:3/10 -c 30 -g 0.8 
Level 7: -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 9:1:0:0/0 -tm 11:1:0:0/0 -tm 13:10:0:0/0 -tm 14:20:1:0/0 -tm 19:20:1:3/10 -c 20 -g 0.8 
Level 8: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 9:1:0:0/0 -tm 11:1:0:0/0 -tm 13:10:0:0/0 -tm 14:20:1:0/0 -tm 19:20:1:3/10 -c 30 -g 0.8 
Level 9: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 9:1:0:0/0 -tm 11:1:0:0/0 -tm 13:10:0:0/0 -tm 14:20:1:0/0 -tm 19:20:1:3/10 -c 40 -g 0.8 
Level 10: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 9:1:0:0/0 -tm 11:1:0:0/0 -tm 13:10:0:0/0 -tm 14:20:1:0/0 -tm 19:20:1:3/10 -c 50 -g 0.8
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

                                                    

