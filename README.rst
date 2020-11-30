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

Introduction
############
**stripenn** is a command line interface python package developed for detection of atchitectural stripes from 3D genome conformation data. Then what are 3D genome conformation data and stripes?

* **3D genome conformation data**
    The human DNA is about 2m in length and it is folded within small cells. However, the folding pattern is not random and this DNA 3D conformation is important in gene regulation. There are several sequencing techniques to see the genome-wide DNA interactions such as `Hi-C <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2858594/>`_ , `HiChIP <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501173/>`_,  `PLAC-seq <https://www.nature.com/articles/cr2016137>`_ and `micro-C <https://www.cell.com/fulltext/S0092-8674(15)00638-8>`_. Currently, imaging-based methods are also developed to construct the spatial distance map of DNA ( `Su et al., 2020 <https://www.sciencedirect.com/science/article/pii/S0092867420309405>`_ ). From these data, special structures such as compartment, TAD and loops have been revealed.

* **Architectural stripe**
    Architectural stripes are one of the distinctive features found from 3D genome conformation data. It was first reported from `Vian et al., Cell, 2018 <https://www.sciencedirect.com/science/article/pii/S0092867418304045>`_. As shown in the figure above, it has anchor domain and contiguous regions are subsequently interacting with the anchor. It is known that stripes contains a large number of super-enhancers, and also it is a hopspot for DNA damage. However, more systematic analysis of this feature is needed, and that's why we developed Stripenn. :)

Installation
############
**Linux / Mac**

Open terminal and type following:
::

    pip install stripenn

If pip is not installed, install pip first: `pip installation <https://pip.pypa.io/en/stable/installing/>`_

When pip is not working, replace it with 'pip3'

**Windows**

Windows users will fail to install stripenn using the above command. This is because the pypairix package installed with the cooler package does not support Windows. But we can still install stripenn on Windows! The strategy is (1) to install all dependency packages except pypairix and then (2) to install stripenn without dependency packages. To do this, type followings on command windows:
::

   pip install simplejson multiprocess pandas numpy matplotlib opencv-python scikit-image scipy joblib tqdm Bottleneck typer pathlib
   pip install --no-deps cooler
   pip install --no-deps stripenn

If you have any trouble with installation, please leave it on issue board.

Quick start (Example run)
#########################
Let's check if stripenn is working or not with a simple example.
::

   cd <Test_Directory> # Move to your test directory
   wget https://data.4dnucleome.org/files-processed/4DNFISA93XFU/@@download/4DNFISA93XFU.mcool -O Vian_aB_30hrs.mcool
   stripenn compute --cool ../hic/Vian_aB_30hrs.mcool::resolutions/5000 --out output_dir/ -k chr19 -m 0.99

Here, the example mcool file contains Hi-C data of mouse activated B cell (`Vian et al., Cell, 2018 <https://www.sciencedirect.com/science/article/pii/S0092867418304045>`_).

Stripes are searched from chromosome 19 of 5kb-resolution data for short running time.

**Output**
::

   cat output_dir/result_filtered.txt

will show a result table including 12 columns like below.

.. csv-table:: result_filtered.txt
   :header: "chr", "pos1","pos2","chr2","pos3","pos4","length","width","Mean","maxpixel","pvalue","Stripiness"

    "chr19", "24285001", "24320000", "chr19", "24285001", "24640000", "355000", "35000", "10.3924", "98.0%", "0.04653", "3.6686"
    "chr19", "6305001", "6345000", "chr19", "6125001", "6345000", "220000", "40000", "17.5088", "98.0%", "0.059462", "2.0324"
    "chr19", "53905001", "53940000", "chr19", "53750001", "53940000", "190000", "35000", "15.7701", "98.0%", "0.03981", "0.5934"

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

1) compute: It is main function of stripenn that detects stripes using image-processing method. There are several options in it.

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

  -c, --canny FLOAT       Canny edge detection parameter.  [default: 2.5]
  -l, --minL INTEGER      Minimum length of stripe.  [default: 10]
  -w, --maxW INTEGER      Maximum width of stripe.  [default: 8] --> we recommend to adjust it to 16 using 5kb-resolution data
  -m, --maxpixel TEXT     Percentiles of the contact frequency data to
                          saturate the image. Separated by comma  [default:
                          0.95,0.96,0.97,0.98,0.99]

  -n, --numcores INTEGER  The number of cores will be used.  [default: 40]
  -p, --pvalue FLOAT      P-value cutoff for stripe.  [default: 0.1]
  --help                  Show this message and exit.

2) score: It calculates p-value and stripiness of given stripes on given 3D genome conformation data. It is useful to compare stripiness of given stripes in two datasets.

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
  --help                  Show this message and exit.

3) seeimage: This function was included to help users choose proper maximum-pixel-value.

  --cool TEXT           Path to cool file  [required]
  -p, --position TEXT   Genomic position (e.g., chr1:135010000-136000000)
                        [required]

  -m, --maxpixel FLOAT  Quantile for the pixel saturation. (e.g., 0.95)
                        [default: 0.95]

  -o, --out TEXT        Path to output directory  [default: ./heatmap.png]
  --norm TEXT           Normalization method. It should be one of the column
                        name of Cooler.bin(). Check it with
                        Cooler.bins().columns (e.g., KR, VC, VC_SQRT)
                        [default: KR]

  --help                Show this message and exit.

