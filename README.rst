========
STRIPENN (ver 1.1.65)
========
--------------------------------------------------------------------
Architectural stripe detection from 3D genome conformation data
--------------------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/pip.svg
   :target: https://pypi.org/project/stripenn/

.. image:: https://github.com/ysora/stripenn/blob/main/image/example_call.png
   :height: 500px
   :width: 500px

Contents
########
* Introduction
* Installation
* Quick start
* Usage

Introduction
############
**Stripenn** is a command line interface python package developed for detection of atchitectural stripes from chromatin conformation capture (3C) data. Then what are 3C data and stripes?

* **Chromatin conformation capture technique**
    Basically, the chromatin conformation capture technique measures the frequency of DNA interactions. There are several sequencing techniques to see the genome-wide DNA interactions such as `Hi-C <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2858594/>`_ , `HiChIP <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501173/>`_,  `PLAC-seq <https://www.nature.com/articles/cr2016137>`_ and `micro-C <https://www.cell.com/fulltext/S0092-8674(15)00638-8>`_. Currently, imaging-based methods are also developed to construct the spatial distance map of DNA interaction ( `Su et al., 2020 <https://www.sciencedirect.com/science/article/pii/S0092867420309405>`_ ). These data have revealed that our genome is highly organized including special features such as compartment, TAD, loops and stripes.

* **Architectural stripe**
    Architectural stripes are shown in the figure above (highlighted in green border). Stripe is one of the distinctive features found from 3C data. It was first predicted by `Fudenberg et al., Cell Reports, 2016 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4889513/>`_. As shown in the figure above, it is a linear pattern where an anchor forms interaction with contiguous genomic region. It is known that stripes are enriched with super-enhancers, and also it is a hopspot for DNA damage. To address more characteristics of stripes, we need accurate stripe caller and that's why we've developed Stripenn.

Requirement
############
**Python 3.8 or higher version**

Installation
############
**Linux / Mac**

Open terminal and type following:
::

    pip install stripenn

If pip is not installed, install pip first: `pip installation <https://pip.pypa.io/en/stable/installing/>`_

When pip is not working, replace it with 'pip3'. Installation may take several minutes.

**Windows**

Windows users will fail to install stripenn using the above command. This is because the pypairix package installed with the cooler package does not support Windows. But we can still install stripenn on Windows! The strategy is (1) to install all dependency packages except pypairix and then (2) to install stripenn without dependency packages. To do this, type followings on command windows:
::

   pip install simplejson multiprocess pandas numpy matplotlib opencv-python scikit-image scipy joblib tqdm Bottleneck typer pathlib
   pip install --no-deps cooler
   pip install --no-deps stripenn

If you have any trouble with installation, please leave it on issue board.

Quick start (Example run)
#########################
Let's check if stripenn is working or not with a simple example which might take less than 5 minutes. This example .cool file is Smc1-HiChIP of only chr16 of mouse T lymphocyte (`Fasolino et al., Immunity, 2020 <https://www.sciencedirect.com/science/article/pii/S1074761320300303>`_).
::

   cd <Test_Directory> # Move to your test directory
   wget https://www.dropbox.com/s/1bb2npvrzp3by5y/BL6.DPT.chr16.mcool?dl=0 -O test.mcool --no-check-certificate
   stripenn compute --cool test.mcool::resolutions/5000 --out output_dir/ -k 16 -m 0.95,0.96,0.97,0.98,0.99

*Tip*: For those whose computer has not enough memory (e.g., < 24GB), we provide slow version of Stripenn as follow:
::

   stripenn compute --cool test.mcool::resolutions/5000 --out output_dir/ -k 16 -m 0.95,0.96,0.97,0.98,0.99 -s

In this example, stripes are searched from chromosome 16 of 5kb-resolution data for short running time.

**Output**
::

   cat output_dir/result_filtered.txt

will show a result table including 12 columns like below.

.. csv-table:: result_filtered.txt
   :header: "chr", "pos1","pos2","chr2","pos3","pos4","length","width","Mean","maxpixel","pvalue","Stripiness"

    "16", "23970001", "24015000", "16", "23970001", "24870000", "900000", "45000", "2.913", "98.0%", "0.025", "10.809"
    "16", "10550001", "10595000", "16", "10550001", "10870000", "320000", "45000", "6.638", "99.0%", "0.042", "10.020"
    "16", "32490001", "32525000", "16", "32490001", "32700000", "210000", "35000", "8.117", "98.0%", "0.012", "5.254"

Each line represents the coordinates and other information of a vertical stripe.

.. image:: https://github.com/ysora/stripenn/blob/main/image/readme1.png

* chr: chromosome
* pos1: x1 position of stripe
* pos2: x2 position of stripe
* chr2: chromosome (same as chr_1)
* pos3: y1 position of stripe
* pos4: y2 position of stripe
* length: Length of vertical stripe (y2-y1+1)
* width: Width of vertical stripe (x2-x1+1)
* mean: Average pixel value within stripe
* maxpixel: Pixel saturation parameter. See below.
* pvalue: P-value of the stripe
* Stripiness: Score of the stripe

Usage
#####

Stripenn has three functions.

* compute
* score
* seeimage

**compute**
:It is main function of stripenn that detects stripes using image-processing method. There are several options in it.

Options:
  --cool TEXT             Path to cool file  [required]
  -o, --out TEXT          Path to output directory  [required]
  --norm TEXT             Normalization method. It should be one of the column
                          name of Cooler.bin(). Check it with
                          Cooler.bins().columns (e.g., KR, VC, VC_SQRT)
                          [default: KR]

  -k, --chrom TEXT        Set of chromosomes. e.g., 'chr1,chr2,chr3', 'all'
                          will generate stripes from all chromosomes
                          [default: all]

  -c, --canny FLOAT       Canny edge detection parameter.  [default: 2.0]
  -l, --minL INTEGER      Minimum length of stripe.  [default: 10]
  -w, --maxW INTEGER      Maximum width of stripe.  [default: 8]
  -m, --maxpixel TEXT     Percentiles of the contact frequency data to
                          saturate the image. Separated by comma  [default:
                          0.95,0.96,0.97,0.98,0.99]

  -n, --numcores INTEGER  The number of cores will be used.  [default: 40]
  -p, --pvalue FLOAT      P-value cutoff for stripe.  [default: 0.1]
  --mask TEXT             Column coordinates to be masked. e.g.,
                          chr9:12345678-12345789  [default: 0]

  -s BOOLEAN_FLAG         Use this if system memory is low.  [default: False]
  -b --bfilter INTEGER    Kernel size of Mean filter.  [default: 3]


**score**
:It calculates p-value and stripiness of given stripes on given 3D genome conformation data. It is useful to compare stripiness of given stripes in two datasets.

Options:
  --cool TEXT             Path to cool file  [required]
  -c, --coord TEXT        Path to stripe coordinate table  [required]
  --norm TEXT             Normalization method. It should be one of the column
                          name of Cooler.bin(). Check it with
                          Cooler.bins().columns (e.g., KR, VC, VC_SQRT)
                          [default: KR]

  -h, --header            Does the stripe coordinate table have header?
                          [default: False]

  -n, --numcores INTEGER  The number of cores will be used.  [default: 40]
  -o, --out TEXT          Path to output file  [default: scores.out]
  
   In the output, four columns (O_Mean_added, O_Sum_added, O/E_Mean_added and O/E_Total_added) are added in addition to the stripiness and P-value, and each represents (1) Mean of observed contact frequency, (2) Total sum of observed contact frequency, (3) Mean of observed/expected contact frequency and (4) total sum of observed/expected contact frequency within stripe. 

**seeimage**
:This function was included to help users choose proper maximum-pixel-value. It draws heatmap image of given position for given maximum pixel paramter.

Options:
  --cool TEXT          Path to cool file  [required]
  -p, --position TEXT  Genomic position (e.g., chr1:135010000-136000000)
                       [required]

  -m, --maxpixel TEXT  Quantile for the pixel saturation. (e.g., 0.95)
                       [default: 0.95,0.96,0.97,0.98,0.99]

  -o, --out TEXT       Output prefix  [default: ./my_heatmap]
  --norm TEXT          Normalization method. It should be one of the column
                       name of Cooler.bin(). Check it with
                       Cooler.bins().columns (e.g., KR, VC, VC_SQRT)
                       [default: KR]

  -s                   Use if system memory is low.  [default: False]


