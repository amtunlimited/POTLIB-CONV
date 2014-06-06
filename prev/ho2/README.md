#Potlib Conversion Tutorial

This tutorial is a description of a process for converting potential energy subroutines into a potlib-compliant program. It is a modified version of the process developed by Dr. Al Wagner. 

The process, in general, takes these stages:

1. Make a test case and compile native code
3. Transfer code into potlib template
4. Check if potlib compiles
5. Run tester and see if you get same results.

I'll try to outline all of the steps. This document is meant to be both a reference and tutorial. It is assumed you have reasonable knowledge of chemistry (mostly for generating the test case) and have a good handle on FORTRAN, both writing and  compiling, the make system, and that you have all of the template files that come with this tutorial.

###The Makefile

The Make sstem is a was to help automate the compiling of code.The Makefile provided with this template, while not strictly necessary, will make things a little easier. I will be referring to it throughout this tutorial.

In the makefile there are specific variables that need to be set. The `FC` variable needs to be set to the FORTRAN compiler of the computer you're using (`gfortran` by default), the `SYSTEM` variable need to be set to the system you're working on (i.e `FH2` or `H2O`), and the `FULLNAME` variable needs to be set to the full name of the subroutine you're working with (i.e. `FH25sec`).

You'll also need to follow naming conventions. The original code should be named `$FULLNAME_original.f`, the test code should be named `origintest.f`, the potlib-compliant code should be named `$FULLNAME.f` and your input file for the potlib tester code should be named `$SYSTEM.dat`

##Making a Test Case

This step is more chemistry than computer science, so I won't be able to go in quite as much detail, but I'll try to lay down a good starting point. In general, you want to find the stable geometry for the molecule you're working with and subtly change it. You don't want the energy to be zero, but you want it to reasonably close. This way, if the answer is too high or zero (or somehow negative), you'll know something is wrong. There are also many reasonable test cases in the POTLIB repository for systems that are already in use.

Make sure you save the input and output for this case, as this is what we will be testing against when we have the completed potlib code. If you use the provided makefile, running "make origin" will compile `$FULLNAME_original.f` and `nativetest.f` code to `origintest`, then run the program and save the output to `origintest.out`

##Making the code POTLIB-compliant

Once you've confirmed that the native code runs correctly, you can start making the code POTLIB compliant. It's highly suggested that you use the template.f file provided, and it's assumed you'll be using it for the purposes of this tutorial.

There are five sections of the template.f file:

1. The Header
2. Prepot
3. Pot
4. PTPACM
5. Additional Subroutines

I'll now go over them in detail.

###The Header

The header is a series of comments at the beginning of the file. It contains the reference to the original POTLIB specification and fields for the name of the system, the full name of the method, the internal coordinates it uses, and any other uses it may have. It also has a place for the original subroutine's reference as well as a place for the transcriber (that is, *you*) to put your name and contact info.

While this is a small section, it is still important to document these things, so this section **should not be skipped**.

###Prepot

This subroutine is called once before any calls to `pot()` are made. It performs various initialization functions and any calculation that only need to be done once.

The first section to look at are the common blocks. Because of how inheritance and scope work in FORTRAN, it is important that you make a common block (I usually use `POTCM`, but it's arbitrary) for *any and all* variables that are initialized in `prepot()` to be used in `pot()`. Also, if any common blocks are used in the original code, add them here too.

Next are a series of variable definitions. These variables are both used as meta-data and for later conversion methods, so it's important that this data be correct. The first definition, `REF`, is a text array meant to store the reference for the original subroutine.

The next two variables, `IRCTNT` and `INDEXES`, are important for the internal coordinate system. The input and output of POTLIB is only done with Cartesian coordinates, so these variables, along with a few function is `utility.f`, are used to convert the Cartesian inputs to the internal coordinates (which will be defined later). For now is it important to know what to store in these variables. `INDEXES` stores the atomic number of each atom in the same order the inputs are made. For example, BrHCl would be `{35,1,17}`. `IRCTNT` is used to define where the split would be if these atoms were in a bimolecular reaction. Atoms 1 to IRCTNT-1 are in reactant one and IRCTNT to NATOMS (a variable containing the number of atoms, which will be discussed later) are in reactant 2.

Next, after the calls to `POTINFO` and `ANCNVRT`, is where any code that needs to only be run once goes. This includes calculation and initializations. After you've added that, you can proceed to the next section.

###Pot

This subroutine is called everytime a calculation needs to be made. This part is reletively staright forward. First, you copy all of the nessesary common blocks into place. This includes any common blocks you created to get data from prepot as well as any common blocks from the original code. Then you copy all of the code used in the calculation (minus any extra subroutines, there's a seperate spot for those) into the marked section.

The most important (and difficult part) is making sure all of the inputs are received from `R`, the energy is stored in `ENGYGS`, and the derivatives (if any are calculated) are stored in `DEGSDR`. From there it's just making sure that none of the variable names aren't the same as variables used by POTLIB (it comes up surprisingly often).

###PTPACM

This section is where all of the DATA sections of the program go. This section also sets a few options for the the conversion functions. The options that need to be set are `NATOMS`, `ICARTR`, `MDER`, and `MSURF`.

`NATOMS` is the number of atoms in the system. `ICARTR` tells the the conversion function what coordinate system you're using:

* ICARTR = 0, no action is taken, the array R not filled
* ICARTR = 1, R = CARTNU, i.e., POT routine uses Cartesian coordinates
* ICARTR = 2, R = r<sub>ij</sub> in canonical order: R(1)=r<sub>12</sub>, R(2)=r<sub>13</sub>,…,R(NATOMS)=r<sub>1,NATOMS</sub>
* ICARTR = 3 (only for 3-atom systems), R(1) = r<sub>12</sub>, R(2) = r<sub>23</sub>, R(3) = r<sub>13</sub> 
* ICARTR = 4 (only for 4-atom systems), R(1) = r<sub>com</sub> (center of mass distance between diatom 1-2 and 3-4), R(2) = r<sub>12</sub>, R(3) = r<sub>34</sub>,  R(4) = r<sub>12</sub>-r<sub>com</sub> angle, R(5) = r<sub>34</sub>-r<sub>com</sub> angle, R(6) = dihedral angle

`MDER` is how many derivatives the code can handle, and `MSURF` is how many excited state energies the code can handle.

The place where all of the data statements should go is marked in the comments. Be sure to also include any common blocks that are necessary.

###Additional Subroutines

There is a place at the end of the template for additional subroutines. It's very important that any subroutines in the original code be copied over. Please also take care that they are not the same name as any of the functions in `template.f` or `utility.f`

##Testing the Converted Code

Once you've converted all the original code, you must test the code against the original code. **YOUR TASK IS NOT COMPLETE UNTIL YOU'VE DONE THIS STEP!** 

First you'll need to make a test data file. The format is:

* 3 x (number of atoms)
* each cartesian coordinate of each atom, one line at a time.
* (this is all in one line) # of atoms, # of derivatives (if unknown put 0), the excited states of the atom (if unknown, put `1 35*0`, then the nflags (if unknown, put `1 1 0 0 0`)

So an example of a `.dat` file would be the following:

    9     
    -0.91666666667d0
    1.77756075064d0
    0.0d0
    0.0d0
    0.0d0
    0.0d0
    1.5d0
    0.0d0
    0.0d0
    3 0 1 35*0 1 1 0 0 0

If you're using the provided makefile, executing `make potlibtest` will try to compile the code along with `utility.f` and `tester.f`. If you're like me, it will inevitably fail the first time. This is where you fix any syntax or glaring errors. The three problems I usually run into are not putting variables that were initialized in the `potlib` or `ptpacm` in a common block so it can transfer, something was named the same as something in the potlib code, or I didn't name one of the potlib variables correctly.

Once you can compile your converted code and get the same answers as the the original code, congratulate yourself because you're done!

If you're not already involved with the project, please feel free to fork this repo, convert some of this code, and make a pull request. If you have any questions you can email me at <aaron.tagliaboschi@gmail.com>