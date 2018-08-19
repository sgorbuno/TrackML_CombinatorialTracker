Hello!  

Below you can find a outline of how to reproduce my solution for the TrackML competition.  
If you run into any trouble with the setup/code or have any questions please contact me at sergey.gorbuov.32@gmail.com  

#ARCHIVE CONTENTS  
Makefile            : file which steers the compilation  
*.cxx *.h           : reconstruction code  
analysis/*.cxx/*.h  : additional code which is used to analyse the data  
geoLayerField.txt,  
geoLayerSizes.txt   : configuration files with detector geometry and magnetic field approximation  


#HARDWARE: (The following specs were used to create the original solution)  
Mac OsX 10.13.6  
2,6 GHz Intel Core i5  
8 GB 1600 MHz DDR3  


#SOFTWARE:  

ROOT 6.14/00 or any other version of ROOT   
"root-config" script should be accesible in the command line   

Web-site of the package: https://root.cern.ch/   
Installation instructions can be found here: https://root.cern.ch/downloading-root  

The package is only used to analyse the data. The dependency can be removed by commenting out   
all the places in the code where TNtuple and TFile classes are used.  


#DATA SETUP (assumes the [Kaggle API](https://github.com/Kaggle/kaggle-api) is installed)  
# below are the shell commands used in each step, as run from the top level directory  

mkdir -p data/  
cd data/  
kaggle competitions download -c trackml-particle-identification -f train_sample.zip  
kaggle competitions download -c trackml-particle-identification -f test.zip  
unzip train_sample.zip  
unzip test.zip  

#DATA PROCESSING  
# The train/predict code will also call this script if it has not already been run on the relevant data.  

1. specify the path to the data and N of events to process in reconstruction.cxx file.   
By default, it is "data/train_100_events"  
One can also set there a flag for analysing the reconstruction efficiency with the truth data:  
bool analyseTruth = true;  

2. compile the code via:  
make  

3. run the reconstruction by typing:  
./reco  

the reconstruction output will be written to mysubmission.csv file  


#MODEL BUILD:  

The algorithm is combinatorial, there is no usual model build step.   
But it has some parameters, which are extracted from the test samples and stored in    
geoLayerField.txt and  geoLayerSizes.txt files.   

In order to rebuid these parameters, one should uncomment   
  //tracker.analyzeGeometry(0);     
  //continue;  
and  
  //tracker.analyzeGeometry(1);  
lines in the reconstruction.cxx  