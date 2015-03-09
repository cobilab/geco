# GeCo #
<p align="center"><img src="/logo.png" 
alt="EAGLE" width="350" height="260" border="0" /></p>
Compress and analyze genomic sequences. As a compression tool, GeCo is able to provide additional compression gains over several top specific tools, while as an analysis tool, GeCo is hable to determine absolute measures, namely for many distance computations, and local measures, such as the information content contained in each element, providing a way to quantify and locate specific genomic events. GeCo can afford individual compression and referential compression.

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

                                                    

