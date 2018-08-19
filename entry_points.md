
1. specify the path to the data and N of events to process in reconstruction.cxx file. 
By default, it is "data/train_100_events"

One can also set there a flag for analysing the reconstruction efficiency with the truth data:
bool analyseTruth = true;

2. compile the code via:
make

3. run the reconstruction by typing:
./reco

the reconstruction output will be written to mysubmission.csv file
