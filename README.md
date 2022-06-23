# pgTree repository
Welcome to the **pyGimli Tree**  repository. 
Scope of the project: 

* Forward and inverse modeling of realistic tree data for both 2D and 3D, using prior information on tree geometry



Table of Contents
--------
* [Introduction](README.md#introduction)
* [Installation](README.md#installation)
    * [Conflicting packages](README.md#conflicting-packages)
* [Download sample data](README.md#download-sample-data)
* [Project Development](README.md#project-development)
    * [Maintainers](README.md#maintainers)
* [License, use and attributions](README.md#license-use-and-attribution)
    

Introduction
-----------

More info about project here


Installation 
-----

First, you will need a new Python 3 environment. The easiest way is to use 
[Anaconda](https://www.anaconda.com/distribution/). 

In the following we provide detailed installation instructions.
1. Create a new anaconda environment
```
conda create -n pg
```
2. Activate the environment and install pygimli. 
3. Remember to activate the environment every time you want to work on the project.
```
conda activate pg

conda install -c gimli pygimli 
```

3. Clone the repository:
```
git clone https://github.com/shakasaki/pgTree.git

cd pgTree
```
5. Install the package locally by using the file `setup.py` by doing:
```
python -m pip install -e .
```

Download data from polybox using pooch
-------

In the terminal type:
```bash
python core/download_data.py
```

Running this code a second time will not trigger a download since the file already exists.


Models and Tutorials
-------------------



Project Development
-------------------

Write about development


### Maintainers

 
License, use and attribution
----------------------------

Feel free to download and use the Bedretto general model in your work! We do not provide any
warranty and any guarantee for the use, but we aim to help and answer questions posted as Issues on the repository 
page as quickly as possible.

The model is published under an **GNU Lesser General Public License v3.0**, which
means that you are free to use it. If you do any modifications, especially for scientific and educational use,
then please _provide them back to the main project_ in the form of a pull request.

If you have questions on the procedure, feel free to contact us about it.
