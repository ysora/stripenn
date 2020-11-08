========
STRIPENN
========
--------------------------------------------------------------------
Architectural stripe detection from 3D genome conformation data
--------------------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/pip.svg
   :target: https://pypi.org/project/stripenn/
   
.. image:: https://github.com/ysora/stripenn/blob/main/image/example_call.png
   :height: 20px
   :width: 20px

Contents
########
* Introduction
* Installation
* Quick start
* Usage
* Stripenn algorithm

Introduction
############
**stripenn** is a command line interface python package developed for detection of atchitectural stripes from 3D genome conformation data. Then what are 3D genome conformation data and stripes?

* **3D genome conformation data**


* **Architectural stripe**

Installation
############
**Linux / Mac**

Open terminal and type following:
::

    pip install stripenn
    
If pip is not installed, install pip first: `pip installation <https://pip.pypa.io/en/stable/installing/>`_

**Windows**

Windows users will fail to install stripenn using the above command. This is because the pypairix package installed with the cooler package does not support Windows. But we can still install stripenn on Windows! The strategy is (1) to install all dependency packages except pypairix and then (2) to install stripenn without dependency packages. To do this, type followings on command windows:
::

   pip install simplejson multiprocess pandas numpy matplotlib opencv-python scikit-image scipy joblib tqdm Bottleneck
   pip install --no-deps cooler
   pip install --no-deps stripenn

If you have any trouble with installation, please leave it on issue board.

Quick start (Example run)
#########################
Let's check if stripenn is working or not with a simple example.
::

   cd <Test_Directory> # Move to your test directory
   wget https://data.4dnucleome.org/files-processed/4DNFISA93XFU/@@download/4DNFISA93XFU.mcool -O Vian_aB_30hrs.mcool
   stripenn Vian_aB_30hrs.mcool::resolutions/5000 output_dir/ -k chr19 -m 0.99
   
Here, the example mcool file contains Hi-C data of mouse activated B cell (`Vian et al., Cell, 2018 <https://www.sciencedirect.com/science/article/pii/S0092867418304045>`_).

Stripes were searched from chromosome 19 of 5kb-resolution data for short running time.

**Output**
::

   cat output_dir/result_filtered.txt
   
will show a result table including 12 columns like below.

.. csv-table:: result_filtered.txt
   :header: "chr_1", "pos1","pos2","chr_2","pos3","pos4","length","width","Mean","maxpixel","pvalue","Stripiness"
   
   "chr19", "-", "-", "chr19", "-", "-", "-", "-", "-", "99.0%", "-", "-"
   "chr19", "-", "-", "chr19", "-", "-", "-", "-", "-", "99.0%", "-", "-"
   "chr19", "-", "-", "chr19", "-", "-", "-", "-", "-", "99.0%", "-", "-"

Each line represents the coordinates and other information of a vertical stripe.

.. image:: https://github.com/ysora/stripenn/blob/main/image/readme1.png

* chr_1: chromosome
* pos1: x1 position of stripe
* pos2: x2 position of stripe
* chr_2: chromosome (same as chr_1)
* pos3: y1 position of stripe
* pos4: y2 position of stripe
* length: Length of vertical stripe (y2-y1+1)
* width: Width of vertical stripe (x2-x1+1)
* mean: Average pixel value within stripe
* maxpixel: 뭐라고 해야되지? 
* pvalue: P-value of the stripe
* Stripiness: Score of the stripe
