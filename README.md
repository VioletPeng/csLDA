# csLDA
The implementation of the paper:

**Learning Polylingual Topic Models from Code-Switched Social Media Documents**  
Nanyun Peng, Yiming Wang and Mark Dredze  
*The Annual Meeting of the Association for Computational Linguistics (ACL)*, 2014  

If you use the code, please kindly cite the following bibtex:

@inproceedings{peng2014learning,  
  title={Learning Polylingual Topic Models from Code-Switched Social Media Documents.},  
  author={Peng, Nanyun and Wang, Yiming and Dredze, Mark},  
  booktitle={ACL (2)},  
  pages={674--679},  
  year={2014}  
}

Implementation adapted from GibbsLDA++, thus use the same liscence. The README of GibbsLDA++ are attached for your reference.

## Compile 
the same as GibbsLDA++.

cd src

make

## Input data format 
almost the same except we require an additional indicator on the language.
    
    [M]
    
    [document_1]
    
    [document_2]
    
    ...
    
    [document_M]

  in which the first line is the total number for documents [M]. Each line 
  after that is one document. [document_i] is the i^th document of the dataset 
  that consists of a list of Ni words/terms.

    [document_i] = Language [word_i1] [word_i2] ... [word_iNi]

  in which Language represents the language of this document. ``codeS" represents code-switched document. All [word_ij] (i=1..M, j=1..Ni) are text strings and they are separated by the space character.

## Run the code
The running option and parameters are mostly inherent from GibbsLDA++, we also added some specific parameters. Please see the code in lda.cpp for more details of the options. run.sh provides an example.

