# Risk Factors Combinations

This program runs several combinations of factors and compare them, based on the R squared adjusted, to return the best ones. The user can customize the number of factors to include and can request specific models as well.

## Content
- The code is composed by 5 sections:
    Initialization - we upload the starting dataframe from an excel file with all the useful data
    1. In this section the main vectors of the code are defined. The user has to fill up some inputs (such as the number of days he wants in the vectors and the
       R squared target). The vector ENPACL represents our y vector of the linear regression.
    2. This section is a control section, it looks for NaN values inside the vectors and assures that there's none
    3. The third section requires the user to choose which kind of analysis he wants. He can choose either all the possible combinations for a certain
       number of factors or a customized model.
    4. Here are defined all the functions that are used in section 3. These functions are three: one for all the possible combinations of a certain number of
       risk factors (AllCombiations), one for the customized model (ModelloPersonalizzato) and one for the output that has to be shown in the terminal (Operazioni)   
- This code helps the user to run hundreds of combinatons of risk factors to find the best models

## How to use it
The fudamental part of this code is to have an Excel file with all the necessary data. In this Excel file we need a sheet called "Dataset" which contains our time series
of our y vector on the second column and a matrix of our risk factors. The user will be asked to insert tha starting column, meaning he has to indicate which column of the excel file
contains the first risk factor of the matrix. He will also be asked to write how many factors are there in the matrix.
Once the user has done this, he just has to run the matllab code, setting an R squared target to get only the models that have a better R squared.
