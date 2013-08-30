ST-analyzer
===========

ST-analyzer is a standalone GUI toolset to perform various analyses of molecular dynamics simulation trajectories and provides a variety of analysis routines especially focused on membrane systems (e.g., lipid chain order parameter, lipid area, etc). Since trajectory files are generally too large to upload to a remote server, ST-analyzer has been developed in cross-platform by installing it into a server where trajectories are located. Once ST-analyzer is installed, user’s local machines governed by any types of existing OS can access the ST-analyzer through HTTP. ST-analyzer is also hosted at http://im.bioinformatics.ku.edu/st-analyzer/

![screenshot](http://people.eecs.ku.edu/~jjeong/images/STanalyzer/system_diagram.png)

![screenshot](http://people.eecs.ku.edu/~jjeong/images/STanalyzer/ST_Analyzer.png)

Installation
============
Although ST-analyzer is developed in Python codes, to cover the wide range of analysis demands and to maintain the cross-platform characteristics, some external python modules and programs are required. For the users who only need to run particular modules, this documentation clearly distinct ‘required modules’ and ‘optional modules’ by notice them as ‘required’ and ‘optional’ in the rest of context. 
> For more about ST-analzyer, please refer ST-analyzer tutorial through http://im.bioinformatics.ku.edu/st-analyzer/tutorials.html
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


### ST-analyzer (*required)
> * Download ST-analyzer from Git-hub: choose one of methods shown below    
>     * Manual download: https://github.com/stanalyzer/ST-analyzer    
>     * Git clone (using commandline): git clone git@github.com:stanalyzer/ST-analyzer.git

Configuration
=============
Let's assume ST-analyzer is stored into /home/your_account/ST-analyzer/stanalyzer    
At the inside installed (unzipped) directory, you can see files and directories: 'manage.py' 'stanalyzer.db' and directories including 'gui', 'stanalyzer', 'templates', and 'trajectory'.     
* manage.py: required to run Django server 
* stanalyzer.db: database file used for ST-analyzer. ID: admin, Password: 12345 
* gui: diretory containing 'models' and 'views' 
* media: default directory storing the results 
* static: directory storing APIs and background modules
* stanalyzer: directory containing system setup files 
* templates: containing template files for ST-analyzer GUI
* trajectory: containing sample trajectory files

### Checking DB consistency
> At your system command line prompt, use followings:
> * user@stanalyzer> cd /home/your_account/ST-analyzer/stanalyzer 
> * user@stanalyzer> python manage.py syncdb


### Run Django to launch ST-analyzer
> * user@stanalyzer> cd /home/your_account/ST-analyzer/stanalyzer 
> * user@stanalyzer> python manage.py runserver 8000    
>     the number ‘8000’ are used as a port number communicating with ST-analyzer. Thus the port number can be changed


### Forwarding port
> Use ssh configuration to forward port
> * user@stanalyzer> cd ~/.ssh 
> * user@stanalyzer> vi config     
>     Edit 'config' file as following:     
>     Host any_name      
>     HostName your.server.com     
>     LocalForward 8000 127.0.0.1:8000     


### Connecting to ST-analyzer through your web-browser 
> * Open your terminal and connect server where ST-analyzer is installed by using ‘Forwarding port’ described above
> * Connect ST-analyzer through web-browser with http://127.0.0.1:8000    
> * You will see the ST-analyzer login.     
> * Initial account and password are 'admin' and '12345'   
