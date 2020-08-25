# README #

### What is this repository for? ###

This repository contains the processing functions used for the processing of the UK-OSNAP mooring data. All these processing functions come originally from the AMOC-RAPID processing developped at NOC (they were located under the rapid/exec/"cruise"/ directory). At this time only the functions used during the last UK-OSNAP cruise (for ADCP, Nortek and microcat processing) were moved in this toolbox. Part of them were adapted/modified by the author. 

The idea behind that is to separate the core functions from the processing scripts to be able to better track the changes in the code using version control software as GIT by monitoring the core functions only, and for example being able to share the evolution of the processing code efficiently between SAMS and NOC. 

In the future this toolbox aims to gather other functions present like for the BPR, the Rapid widget, etc.... This will allow to keep track and record the changes in the code, and make easier sharing code between researcher inside SAMS but also with other institutes. 

Contributions are more than welcome, if one wants to make comments, create a new branch, adding processing functions for other instruments, widgets, etc...

### How do I get set up? ###

1- Pull the "master" branch and add the .m files to your MATLAB search path by using, for example, the command:


```
#!matlab
addpath(genpath(*path of the toolbox*))

```

2- [Download example processing scripts from here](https://bitbucket.org/Lhoupert/mooring_processing_toolbox/downloads/example_processing_scripts.zip). This will give you an idea how to call the functions.

3- Have a look to the notes about how to go through the different stages of the microcat processing by reading this [summary document](https://bitbucket.org/Lhoupert/mooring_processing_toolbox/downloads/summary_OSNAP_microcat_postprocessing.pdf). The processing of the current meter is detailled in this [processing report](https://bitbucket.org/Lhoupert/mooring_processing_toolbox/downloads/proc_report_curren_meter_OSNAP1_v1.pdf)

4- [An example of the directory structure (with control metadatafiles) can be downloaded here](https://bitbucket.org/Lhoupert/mooring_processing_toolbox/downloads/osnap_directory_structure.zip)  

### Who do I talk to? ###

For questions or comments, you can open a ticket ("Issue" tab in the left panel), or contact Loïc Houpert (loic.houpert at sams.ac.uk)