# GeCo #
<p align="center"><img src="/logo.png" 
alt="EAGLE" width="350" height="260" border="0" /></p>
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
                                                                       
  -h                    give this help,                                
  -x                    show several running examples,                 
  -s                    show GeCo compression levels,                  
  -v                    verbose mode (more information),               
  -V                    display version number,                        
  -f                    force overwrite of output,                     
  -l &#60level&#62            level of compression [1;9] (lazy -tm setup),   
  -g &#60gamma&#62            mixture decayment forgetting factor. It is     
                        a real value in the interval [0;1),            
  -c &#60cache&#62            maximum collisions for hash cache. Memory      
                        values are higly dependent of the parameter    
                        specification,                                 
  -e                    it creates a file with the extension ".iae"  
                        with the respective information content. If    
                        the file is FASTA or FASTQ it will only use    
                        the "ACGT" (genomic) data,                   
  -r &#60FILE&#62             reference file ("-rm" are loaded here),      
                                                                       
Mandatory arguments:                                                   
                                                                       
  -rm &#60c&#62:&#60d&#62:&#60i&#62:&#60m&#62   reference context model (ex:-rm 13:100:0:0),   
  -rm &#60c&#62:&#60d>:&#60i&#62:&#60m&#62   reference context model (ex:-rm 18:1000:0:1),  
  ...                                                                  
  -tm &#60c&#62:&#60d&#62:&#60i&#62:&#60m&#62   target context model (ex:-tm 4:1:0:0),         
  -tm &#60c&#62:&#60d&#62:&#60i&#62:&#60m&#62   target context model (ex:-tm 18:20:1:1),       
  ...                                                                  
                        target and reference templates use &#60c&#62 for     
                        context-order size, &#60d&#62 for alpha (1/&#60d&#62),     
                        &#60i&#62 (0 or 1) to set the usage of inverted      
                        repeats (1 to use) and &#60m&#62 to the maximum      
                        allowed mutation on the context without        
                        being discarded (usefull in deep contexts),    
                                                                       
  &#60FILE&#62                file to compress (last argument). For more     
                        files use splitting ":" characters.          
                                                                       
Report bugs to &#60{pratas,ap,pjf}@ua.pt&#62. 
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

                                                    

