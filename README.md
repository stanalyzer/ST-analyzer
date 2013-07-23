ST-analyzer
===========

ST-analyzer is a standalone GUI toolset to perform various analyses of molecular dynamics simulation trajectories and provides a variety of analysis routines especially focused on membrane systems (e.g., lipid chain order parameter, lipid area, etc). Since trajectory files are generally too large to upload to a remote server, ST-analyzer has been developed in cross-platform by installing it into a server where trajectories are located. Once ST-analyzer is installed, user’s local machines governed by any types of existing OS can access the ST-analyzer through HTTP. ST-analyzer is also hosted at http://im.bioinformatics.ku.edu/st-analyzer/

![screenshot](http://people.eecs.ku.edu/~jjeong/images/STanalyzer/system_diagram.png)

![screenshot](http://people.eecs.ku.edu/~jjeong/images/STanalyzer/ST_Analyzer.png)

Installation
============
Although ST-analyzer is developed in Python codes, to cover the wide range of analysis demands and to maintain the cross-platform characteristics, some external python modules and programs are required. For the users who only need to run particular modules, this documentation clearly distinct ‘required modules’ and ‘optional modules’ by notice them as ‘required’ and ‘optional’ in the rest of context. 

Required modules & programs
----------------------------

### Python (*required)
> Python v2.7 or above is required.  
> Python is available through http://www.python.org/download/


### Django (*required)
> ST-analyzer is optimized Django v1.4.1.    
> Django is available through https://www.djangoproject.com/download/

### MDAnalysis (*required)
> ST-analyzer is optimized MDAnalysis v0.7.6 and above.   
> MDAnalysis requires other modules; therefore, to make the installation simple, we encourage installing all-in-one package. Following list of packages has their own copyright, so please visit their websites and check the eligibility prior to the installation. 

> #### All-in-one package
> Anaconda (http://continuum.io/downloads.html) for Linux, Windows and Mac.   
> Enthought Canopy (https://www.enthought.com/products/canopy/) for Linux, Windows and Mac.

> #### Install MDAnalysis
> Install one of the all-in-one packages listed above     
> Download and install MDAnalysis through https://code.google.com/p/mdanalysis/    
>	Details of installing MDAnalysis can be found in https://code.google.com/p/mdanalysis/wiki/Install    
> We have been reported about problems of the installation. Most problems are caused by outdated version of Python and GNU C compiler (http://gcc.gnu.org/). If you have problems with install, please check the version of your GNU C compiler and Python and discuss with your system administrator
> For more questions about installation, please use discussion board at https://code.google.com/p/mdanalysis/wiki/Install

### Pyhull (*optional)
> Pyhull is Python wrapper to qhull (http://www.qhull.org) used for calculating ‘area-per-lipid’ in ST-analyzer. If calculating area-per-lipid is not necessary, users are not required installing this module.  
> Details of instruction for installing Pyhull can be found at http://pythonhosted.org/pyhull/


