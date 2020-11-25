////////// README LBMS LIBRARY //////////

---------- Introduction ----------

1- LBMS library has been developed for research. None of the authors of the LBMS library claim that it is well suited for
your personal application. Moreover, the respect of PALABOS programming guidelines and C++ good practices are not guaranteed.

---------- Get started with PALABOS files ----------

1- Go to modified PALABOS files folder and replace all files with the same name into PALABOS library.
They contain new dynamics and related external forces implementations. We can notice that none of initial PALABOS
dynamics, classes or functions have been modified. If you want to keep your version of PALABOS (v2.0r0 mandatory),
please copy the classes and implementations in related PALABOS files. This procedure is under your own responsability.
None of the authors will be liable for changes that can occur in your PALABOS native files.

---------- Used version of PALABOS ----------

1- PALABOS v2.0r0.

---------- Compiler options and others ----------

1- We recommend to copy all the LBMS library in the same folder where the folder /src of PALABOS is.
2- You must add new sources files to your compiler search directories.

---------- Dependencies ----------

1- PALABOS v2.0r0.

---------- MPI ----------

1- Only sequential jobs are developed. Using MPI can lead to errors.

---------- IDE ----------

1- We recommend to use CODEBLOCKS IDE (because the project file is provided with PALABOS).
