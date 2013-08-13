rkt-matrix
==========

An interpretation of a matrix operations calculator, based on elementary matrices.

##Introduction
Created and distributed for study purposes for a course in linear algebra.

##Definitions
<p>In this module, a number of definitions are used. </p>
An `L-Vector` is defined as any `(listof Num)`. <br/>
An `L-Matrix` is defined as any `(listof L-Vector)`. <br/>
A `Matrix` is defined as any `(matrix L-Matrix L-Matrix)`. <br/>
<br/>
We use two `L-Matrices` in order to store a n by m matrix and an additional m by n matrix, which we will define as the reverse matrix, in which each entry in the main diagonal is 1. When n = m, the reverse matrix is an identity matrix.
##Basics
###new-matrix
To create a new `Matrix` from an equivalent `L-Matrix`, we apply the `new-matrix` function on it. `NM` is an equivalent short form.<br/>

    Syntax: (new-matrix <L-Matrix>)
            (NM <L-Matrix>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix l-matrix))
     (define A-2 (NM l-matrix))

###display-matrix

To display the contents of the matrix, we use `display-matrix`.

    Syntax: (display-matrix <Matrix>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix l-matrix))
     (display-matrix A)
    (1 2 3)
    (4 5 6)
    (7 8 9)
    
###display-reverse

To display the contents of the reverse matrix, we use `display-reverse`.

    Syntax: (display-reverse <Matrix>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix l-matrix))
     (display-reverse A)
    (1 0 0)
    (0 1 0)
    (0 0 1)
##Elementary Matrices

The following functions produce an elementary `L-Matrix` of the same dimension of the given `L-Matrix`.
###EM-swap

`EM-swap` produces an elementary matrix of row exchange between `row1` and `row2`.

    Syntax: (EM-swap <Matrix> <row1> <row2>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix (EM-swap l-matrix 1 2)))
     (display-matrix A)
    (0 1 0)
    (1 0 0)
    (0 0 1)

###EM-mult

`EM-swap` produces an elementary matrix of scalar row multiplication on `row` by the given `scalar`.

    Syntax: (EM-mult <Matrix> <row> <scalar>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix (EM-mult l-matrix 2 5)))
     (display-matrix A)
    (1 0 0)
    (0 5 0)
    (0 0 1)

###EM-plus

`EM-plus` produces an elementary matrix of scalar row addition of `scalar`*`row2` onto `row1`.

    Syntax: (EM-plus <Matrix> <row1> <row2> <scalar>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix (EM-plus l-matrix 1 2 3)))
     (display-matrix A)
    (1 3 0)
    (0 1 0)
    (0 0 1)

##Basic Operations
The following functions perform simple operations on `L-Matrices` and `Matrices`.
###transpose/transpose!

`transpose` performs matrix transposition on the given `L-Matrix`. 

    Syntax: (transpose <L-Matrix>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (transpose l-matrix)
    '((1 4 7) (2 5 8) (3 6 9))

`transpose!` is the equivalent for `Matrices`.

    Syntax: (transpose! <Matrix>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix (EM-plus l-matrix 1 2 3)))
     (transpose! A)
     (display-matrix A)
    (1 4 7)
    (2 5 8)
    (3 6 9)

###matrix-multiply/matrix-multiply!

`matrix-multiply` performs matrix-matrix multiplication AxB.

    Syntax: (matrix-multiply <L-Matrix> <L-Matrix>)
    >(define l-matrix '((1 5 9)
                        (0 1 0)
                        (0 3 2)))
     (define l-matrix2 '((1 0 0)
                         (0 1 0)
                         (0 3 2)))
     (matrix-multiply l-matrix l-matrix2)
    '((1 32 18) (0 1 0)(0 9 4))

`matrix-multiply!` is the equivalent for `Matrices`.

    Syntax: (matrix-multiply! <Matrix> <Matrix>)
    >(define l-matrix '((1 5 9)
                        (0 1 0)
                        (0 3 2)))
     (define l-matrix2 '((1 0 0)
                         (0 1 0)
                         (0 3 2)))
     (define A (new-matrix l-matrix))
     (define B (new-matrix l-matrix2))
     (matrix-multiply! A B)
     (display-matrix A)
    (1 32 18)
    (0 1 0)
    (0 9 4)

##Matrix Property Operations

The following functions perform more complex operations based on the properties of a `Matrix`.

###ERO-swap!

`ERO-swap!` performs the elementary row operation of row exchange on a `Matrix`, between `row1` and `row2`. The operation detail and resulting matrix is shown. `Quiet-ERO-swap!` is the non-verbose equivalent. `swp` and `qswp` are their short forms.

    Syntax: (ERO-swap! <Matrix> <row1> <row2>)
            (swp <Matrix> <row1> <row2>)
            (Quiet-ERO-swap! <Matrix> <row1> <row2>)
            (qswp <Matrix> <row1> <row2>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix l-matrix))
     (ERO-swap! A 1 2)
     (display-matrix A)
    (4 5 6)
    (1 2 3)
    (7 8 9)

###ERO-mult!

`ERO-mult!` performs the elementary row operation of scalar row multiplication on `row` by the given `scalar`. The operation detail and resulting matrix is shown. `Quiet-ERO-mult!` is the non-verbose equivalent. `mlt` and `qmlt` are their short forms.

    Syntax: (ERO-mult! <Matrix> <row> <scalar>)
            (mlt <Matrix> <row> <scalar>)
            (Quiet-ERO-mult! <Matrix> <row> <scalar>)
            (qmlt <Matrix> <row> <scalar>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix l-matrix))
     (ERO-mult! A 1 2)
     (display-matrix A)
    (2 4 6)
    (4 5 6)
    (7 8 9)

###ERO-plus!

`EM-plus` produces an elementary matrix of scalar row addition of `scalar`*`row2` onto `row1`. The operation detail and resulting matrix is shown. `Quiet-ERO-plus!` is the non-verbose equivalent. `qpls` and `qpls` are their short forms.

    Syntax: (ERO-plus! <Matrix> <row1> <row2> <scalar>)
            (pls <Matrix> <row1> <row2> <scalar>)
            (Quiet-ERO-plus! <Matrix> <row1> <row2> <scalar>)
            (qpls <Matrix> <row1> <row2> <scalar>)
    >(define l-matrix '((1 2 3)
                        (4 5 6)
                        (7 8 9)))
     (define A (new-matrix l-matrix))
     (ERO-plus! A 2 1 -1)
     (display-matrix A)
    (1 2 3)
    (3 3 3)
    (7 8 9)

###row-reduce!/RREF!
`row-reduce` performs a full Gaussian Elimination algorithm on the given `Matrix`. `RREF!` is an equivalent short form.

    Syntax: (row-reduce! <Matrix>)
            (RREF! <Matrix>)
    > (RREF! A)
    Reducing Column 1
    Operation: Multiplying R1 by 1
    (1 2 3)
    (4 5 6)
    (7 8 9)
    Search result: 0
    Operation: Adding -4xR2 to R1
    (1 2 3)
    (0 -3 -6)
    (7 8 9)
    Operation: Adding -7xR3 to R1
    (1 2 3)
    (0 -3 -6)
    (0 -6 -12)
    End of Column 1 Reduction
    Reducing Column 2
    Operation: Multiplying R2 by -1/3
    (1 2 3)
    (0 1 2)
    (0 -6 -12)
    Search result: 1
    Operation: Adding -2xR1 to R2
    (1 0 -1)
    (0 1 2)
    (0 -6 -12)
    Operation: Adding 6xR3 to R2
    (1 0 -1)
    (0 1 2)
    (0 0 0)
    End of Column 2 Reduction
    Reducing Column 3
    Search result: #f
    No suitable rows found for Column 3
    End of Column 3 Reduction
    Row Reduction Complete
    Reduced Row Echelon Form:
    (1 0 -1)
    (0 1 2)
    (0 0 0)

###REF?
`REF?` determines whether or not the given `Matrix` is in Row Echelon Form.

    Syntax: (REF? <Matrix>)
    >(define l-matrix '((1 0 0)
                        (0 2 2)
                        (0 3 7)))
     (define l-matrix2 '((1 0 2)
                         (0 1 1)
                         (0 0 0)))
     (define A (new-matrix l-matrix))
     (define B (new-matrix l-matrix2))
    >(REF? A)
    #f
    >(REF? B)
    #t

###RREF?
`RREF?` determines whether or not the given `Matrix` is in Reduced Row Echelon Form.

    Syntax: (RREF? <Matrix>)
    >(define l-matrix '((1 5 9 0)
                        (0 0 0 0)
                        (0 1 2 3)
                        (0 0 1 0)))
     (define l-matrix2 '((1 0 0 5)
                         (0 1 0 0)
                         (0 0 1 3)
                         (0 0 0 0)))
     (define A (new-matrix l-matrix))
     (define B (new-matrix l-matrix2))
    >(RREF? A)
    #f
    >(RREF? B)
    #t

###rank
`rank` determines the rank of the given `Matrix`.

    Syntax: (rank <Matrix>)
    >(define l-matrix '((1 0 0)
                        (3 1 0)
                        (0 0 0)))
     (define A (new-matrix l-matrix))
     (rank A)
    2
