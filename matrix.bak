#lang racket

;; **********************************************
;; Matrix Calculation
;; Version 1.0
;; By Alex Truong
;; **********************************************

; Introduction:
; Created and distributed for study purposes for a course in linear algebra.
;
;      Operations were based off of the interactions of elementary matrices.
; Thus, the functions in this module are likely not the most efficient, but
; will function very quickly for most typical matrixes students would likely
; type out by hand for the purpose of verification.
;
;      Matrices as defined in this module consist of a mutable structure 
; containing an L-Matrix (a matrix written in the form of a row vector list),
; and a pseudo-identity L-Matrix of equal dimensions. 
;
;     This L-Matrix is used for finding inverses, as any row operations
; performed on one L-Matrix is also performed on the other. Hence, the 
; pseudo-identity L-matrix is merely piggy-back riding. Inverses are invalid
; for any L-Matrix of dimension n by m where n != m.
;
;     A number of functions have been defined as shortcuts for interaction
; purposes. Non-mutative functions are intended for L-Matrices, while the
; mutative functions are intended for Matrices.

(provide 
 matrix
 M
 new-matrix
 NM
 transpose
 transpose!
 matrix-multiply
 matrix-multiply!
 EM-swap
 EM-mult
 EM-plus
 ERO-swap!
 ERO-mult!
 ERO-plus!
 row-reduce!
 RREF!
 display-matrix
 display-reverse)

;; An L-Vector is a (listof Num)
;; An L-Matrix is a (listof Vector)
;; A Matrix is a (matrix L-Matrix L-Matrix)
(define-struct matrix (entries reverse) #:mutable)
;; Shorthand for Interactions
(define M matrix)

;; non-zero? : L-Vector -> Boolean
;; PRE:   true
;; POST:  Returns true if L-Vector a is non-zero.
(define (non-zero? a)
  (cond [(empty? a) true]
        [(zero? (first a)) false]
        [else (non-zero? (rest a))]))

;; build-identity : Nat Nat Nat -> L-Matrix
;; PRE:   true
;; POST:  Produces an L-Matrix of size n by m.
;; Notes: i is used as incrementation.
(define (build-identity i n m)
  (cond [(= m i) empty]
        [else
         (cons (build-list n (lambda (x)
                               (cond [(= x i) 1]
                                     [else 0])))
               (build-identity (add1 i) n m))]))

;; new-matrix : L-Matrix -> Matrix
;; PRE:   true
;; POST:  Produces a new Matrix.
(define (new-matrix A)
  (define n (length A))
  (define m (length (first A)))
  (M A (build-identity 0 n m)))
;; Shorthand for Interactions
(define NM new-matrix)

;; vector-multiply : L-Vector L-Vector -> Num
;; PRE:   true
;; POST:  Returns the vector product of a and b.
(define (vector-multiply a b)
  (cond [(empty? b) 0]
        [else
         (+ (* (first a) (first b))
            (vector-multiply (rest a) (rest b)))]))

;; row-multiply : L-Vector L-Matrix -> Num
;; PRE:   true
;; POST:  Returns the vector-matrix product of a and B.
(define (row-multiply a B)
  (cond [(empty? B) empty]
        [else
         (cons (vector-multiply a (first B))
               (row-multiply a (rest B)))]))

;; transpose : L-Matrix -> L-Matrix
;; PRE:   true
;; POST:  Returns the transpose matrix of A.
(define (transpose A)
  (apply map list A))

;; transpose : Matrix -> Void
;; PRE:   true
;; POST:  Sets A to be the transpose matrix of A.
;;        Also transposes the reverse component of A.
(define (transpose! A)
  (set-matrix-entries! A (apply map list (matrix-entries A)))
  (set-matrix-reverse! A (apply map list (matrix-reverse A))))

;; matrix-multiply : L-Matrix L-Matrix -> L-Matrix
;; PRE:   true
;; POST:  Returns the matrix-matrix product of A and B.
(define (matrix-multiply A B)
  (define (multiply A B)
    (cond [(empty? A) empty]
          [else
           (cons (row-multiply (first A) B)
                 (multiply (rest A) B))]))
  (multiply A (transpose B)))

;; matrix-multiply! : Matrix Matrix -> Void
;; PRE:   true
;; POST:  Sets A to be the matrix-matrix product of A and B.
;; Notes: B remains unchanged.
(define (matrix-multiply! A B)
  (set-matrix-entries! A 
                       (matrix-multiply (matrix-entries A)
                                        (matrix-entries B))))

;; EM-swap : L-Matrix Nat Nat -> L-Matrix
;; PRE:   row1, row2 are within range.
;; POST:  Returns an Elementary L-Matrix of row exchange of the same dimensions
;;        as L-Matrix A, for row1 and row2.
(define (EM-swap A row1 row2)
  (define m (length A))
  (define n (length (first A)))
  (define (entry-ones-fn n)
    (lambda (x)
      (cond [(= x (sub1 n)) 1]
            [else 0]))) 
  (define (row-create i)
    (cond [(= i (sub1 row1))
           (build-list n (entry-ones-fn row2))]
          [(= i (sub1 row2))
           (build-list n (entry-ones-fn row1))]
          [else
           (build-list n (entry-ones-fn (add1 i)))]))
  (define (matrix-create i)
    (cond [(= i m) empty]
          [else
           (cons (row-create i)
                 (matrix-create (add1 i)))]))
  (matrix-create 0))

;; EM-mult : L-Matrix Nat Num -> L-Matrix
;; PRE:   row is within range.
;;        scalar is non-zero.
;; POST:  Returns an Elementary L-Matrix of scalar multiplication of the
;;        same dimensions as L-Matrix A, for scalar*row.
(define (EM-mult A row scalar)
  (define m (length A))
  (define n (length (first A)))
  (define (entry-ones-fn n one)
    (lambda (x)
      (cond [(= x (sub1 n)) one]
            [else 0])))
  (define (row-create i)
    (cond [(= i (sub1 row))
           (build-list n (entry-ones-fn row scalar))]
          [else
           (build-list n (entry-ones-fn (add1 i) 1))]))
  (define (matrix-create i)
    (cond [(= i m) empty]
          [else
           (cons (row-create i)
                 (matrix-create (add1 i)))]))
  (matrix-create 0))

;; EM-plus : L-Matrix Nat Nat Num -> L-Matrix
;; PRE:   row1, row2 are within range.
;;        scalar is non-zero.
;; POST:  Returns an Elementary L-Matrix of row addition of the same dimensions
;         as L-Matrix A, for row1+scalar*row2.
(define (EM-plus A row1 row2 scalar)
  (define m (length A))
  (define n (length (first A)))
  (define (entry-ones-fn n scalar)
    (lambda (x)
      (cond [(and (= n row1)
                  (= x (sub1 row2)))
             scalar]
            [(= x (sub1 n)) 1]
            [else 0])))
  (define (row-create i)
    (cond [(= i (sub1 row1))
           (build-list n (entry-ones-fn row1 scalar))]
          [else
           (build-list n (entry-ones-fn (add1 i) 0))]))
  (define (matrix-create i)
    (cond [(= i m) empty]
          [else
           (cons (row-create i)
                 (matrix-create (add1 i)))]))
  (matrix-create 0))

;; ERO-swap! : Matrix Nat Nat -> Void
;; PRE:   row1, row2 are within range.
;; POST:  Performs the Elementary Row Operation of row exchange on Matrix A,
;         for row1 <-> row2.
(define (ERO-swap! A row1 row2)
  (define EM (EM-swap (matrix-entries A) row1 row2))
  (set-matrix-entries! A (matrix-multiply EM (matrix-entries A)))
  (set-matrix-reverse! A (matrix-multiply EM (matrix-reverse A)))
  (printf "Operation: Swapping R~a <-> R~a\n" row1 row2)
  (display-matrix A))

;; ERO-mult! : Matrix Nat Num -> Void
;; PRE:   row is within range.
;;        scalar is non-zero.
;; POST:  Performs the Elementary Row Operation of scalar multiplication on
;;        Matrix A, for scalar*row.
(define (ERO-mult! A row scalar)
  (define EM (EM-mult (matrix-entries A) row scalar))
  (set-matrix-entries! A (matrix-multiply EM (matrix-entries A)))
  (set-matrix-reverse! A (matrix-multiply EM (matrix-reverse A)))
  (printf "Operation: Multiplying R~a by ~a\n" row scalar)
  (display-matrix A))

;; ERO-plus! : Matrix Nat Nat Num -> Void
;; PRE:   row1, row2 are within range.
;;        scalar is non-zero.
;; POST:  Performs the Elementary Row Operation of row addition on Matrix A,
;         for row1+scalar*row2.
(define (ERO-plus! A row1 row2 scalar)
  (define EM (EM-plus (matrix-entries A) row1 row2 scalar))
  (set-matrix-entries! A (matrix-multiply EM (matrix-entries A)))
  (set-matrix-reverse! A (matrix-multiply EM (matrix-reverse A)))
  (printf "Operation: Adding ~axR~a to R~a\n" scalar row1 row2)
  (display-matrix A))

;; Shorthand for Interactions
(define swp ERO-swap!)
(define mlt ERO-mult!)
(define pls ERO-plus!)

;; display-L-matrix : L-Matrix -> Void
;; PRE:   true
;; POST:  Prints the L-Vectors of L-Matrix A row by row.
;; Notes: Normally displaying a list by itself results in a horizontal
;;        display, which is more difficult to read.
(define (display-L-matrix A)
  (cond [(empty? A) (display "")]
        [else
         (printf "~a\n" (first A))
         (display-L-matrix (rest A))]))

;; display-matrix : Matrix -> Void
;; PRE:   true
;; POST:  Prints the L-Vectors of Matrix A row by row.
;; Notes: Normally displaying a list by itself results in a horizontal
;;        display, which is more difficult to read.
(define (display-matrix A)
  (display-L-matrix (matrix-entries A)))

;; display-reverse : Matrix -> Void
;; PRE:   true
;; POST:  Prints the L-Vectors of Matrix A reverse row by row.
;; Notes: Normally displaying a list by itself results in a horizontal
;;        display, which is more difficult to read.
(define (display-reverse A)
  (display-L-matrix (matrix-reverse A)))

;; row-reduce! : Matrix -> Void
;; PRE:   true
;; POST:  Row-reduces Matrix A.
;; Notes: For the algorithm below, n is used to record the number of
;;        rows.
(define (row-reduce! A)
  (define rows (length (matrix-entries A)))
  ;; row-search : Nat Nat L-Matrix -> (Union false Nat)
  ;; PRE:   true
  ;; POST:  Returns the index of the first suitable row for which the n-th
  ;;        entry is non-zero, or false if such a row does not exist.
  ;;        Performs scalar multiplication on that row, if it exists.
  (define (row-search n i entries)
    (cond [(empty? entries) false]
          [(<= i n)
           (row-search n (add1 i) (rest entries))]
          [else
           (define nth-entry (list-ref (first entries) n))
           (cond [(not (= 0 nth-entry))
                  (mlt A (add1 i) (/ nth-entry))
                  i]
                 [else
                  (row-search n (add1 i) (rest entries))])]))
  ;; row-zero : Nat Nat Nat L-Matrix -> Void
  ;; PRE:   true
  ;; POST:  Performs row reduction by row addition using the row found by
  ;;        row-search, for each row in the L-Matrix.
  (define (row-zero n i origin entries)
    (cond [(empty? entries) void]
          [(= i origin)
           (row-zero n (add1 i) origin (rest entries))]
          [else
           (define elim-co-eff
             (- (list-ref (first entries) n)))
           (pls A (add1 i) (add1 origin) elim-co-eff)
           (row-zero n (add1 i) origin (rest entries))]))
  ;; row-swap : Nat Nat -> Void
  ;; PRE:   true
  ;; POST:  Performs row exchange if necessary.
  (define (row-swap n i)
    (cond [(= n i) void]
          [else
           (swp A (add1 n) (add1 i))]))
  ;; reduce : Nat -> Void
  ;; PRE:   true
  ;; POST:  Enacts the following algorithm:
  ;;        Let n be the number of completed columns via row reductions.
  ;;        Search n-th entry in each column for a suitable row for
  ;;        row addition, via row-search.
  ;;        Performs row reduction if such a row exists for n, via row-zero.
  ;;        Performs row exchange if required, via row-swap.
  ;;        Increment n, and repeat until n = rows of Matrix A.
  (define (reduce n)
    (define entries (matrix-entries A))
    (printf "Reducing Column ~a\n" (add1 n))
    (cond [(= n rows) void]
          [else
           (define search-result (row-search n 0 entries))
           (cond [search-result
                  (row-zero n 0 search-result entries)
                  (row-swap n search-result)
                  (display-matrix A)]
                 [else
                  (printf "No leading numbers found for Column ~a\n" n)])
           (printf "End of Column ~a Reduction\n" n)
           (reduce (add1 n))]))
  (reduce 0)
  (printf "Row-Reduced Eschelon Form:\n")
  (display-matrix A))
;; Shorthand for Interactions
(define RREF! row-reduce!)


(define a '((1 4 7) 
            (2 5 8) 
            (3 6 9)))
(define b '((1 0 0) 
            (0 0 1) 
            (0 1 0)))

(define A (NM a))