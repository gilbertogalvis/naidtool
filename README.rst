#################################################
naidtool - Nonrigid Assited Image Distortion Tool
#################################################

:Authors: John Kay, Gilberto Galvis
:Email: john.kays2020@gmail.com, galvisgilberto@gmail.com
:Version: $revision: 0.1.1 $

This project provides an interface tool for making non-rigid assisted distortions on images. Assistance is possible through a grid of points that can be dragged. The dragging of the points generates the non-rigid distortion on the image. The interface is based on MatLab GUI

Requirements
------------

- Clone this repository on your machine

- MatLab: The interface tool is based on MatLab, so it is required to have MatLab installed

Installation
------------

There are two ways to use the tool: 1) optimized mode and 2) non-optimized mode. The optimized mode is much faster.

Using optimized mode
====================

Requires the compilation of the .mex files found in the functions_nonrigid folder. Then, open MatLab and run the following in the command window

.. code:: matlab

    compile_c_files

Once executed, the optimized functions will be ready to use the optimized mode.

Note: Please note that, depending on your operating system, compiling the .mex files may require some additional steps such as compiler settings, for example.

Using non-optimized mode
========================

You only need to change the path to the files that the interface will use. Please open the nonrigid_assited_image_distortion.m file with the text editor of your choice. Then change line 68 to ``addpath ([functiondir '/ functions_nonrigid_matlab'])``. Save the changes and close the edition.

With this you can already use the tool in non-optimized mode

Usage
-----

To run the tool, you just need to open MatLab and run the following line in the MatLab Command Window

.. code:: matlab

    nonrigid_assited_image_distortion