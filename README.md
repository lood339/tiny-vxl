# tiny-xvl
C++ headfile-only tiny VXL for Computer Vision and Image Understanding  
The main purpose is for easy use of VXL.       

## Introduction: What is (original) VXL?
VXL (the Vision-something-Libraries) is a collection of C++ libraries designed for computer vision research and implementation. It was created from TargetJr and the IUE with the aim of making a light, fast and consistent system. VXL is written in ANSI/ISO C++ and is designed to be portable over many platforms. The core libraries in VXL are:

A comprehensive description of the VXL project can be views at https://vxl.github.io/

# What is the difference between tiny-vxl and the vxl  
1. Only keep vgl, vnl and vpgl from the original VXL project  
2. Only rely on Eigen https://github.com/eigenteam/eigen-git-mirror 
3. The project is on-going, please contact me (jhchen14@cs.ubc.ca) if you would like to contribute.  

# What this project do
1. Move code from .cxx, .cpp to the .h file
2. The original souce code has a good record of author and history, keep them in the file as much as possible.  
4. Rewrite the testing code using the Gtest Framework https://github.com/google/googletest    

# How to use the code
Add the vxl folder to your head file directory 

# Test  
1. Mac OS Mojave (10.14.12)   
    cd build   
    cmake ..  
    make -j4  


 






