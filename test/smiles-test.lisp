
(in-package :chemicl-test)

(5am:def-suite :chemicl-smiles)
(5am:in-suite :chemicl-smiles)

(defun check-elements (molecule element-count-list)
  (equal (sort (count-elements molecule)
               #'string<
               :key (lambda (x) (id (car x))))
         (sort (loop for element-count in element-count-list
                  collect
                  (cons (get-element (car element-count))
                        (cadr element-count)))
               #'string<
               :key (lambda (x) (id (car x))))))

(defun check-smiles-strings (molspecs)
    (map nil
     (lambda (molspec)
       (destructuring-bind (smiles mass elements &optional charge)
           molspec
         (let ((mol (parse-smiles-string smiles)))
           (5am:is (< (abs (- (mass mol) mass .1))))
           (5am:is (check-elements mol elements))
           (5am:is (= (or charge 0) (charge mol))))))
     molspecs))


;;; Simple one non-carbon atom molecules containing an atom from the
;;; "organic" subset and the appropriate number of hydrogens.

(5am:test
    (chemicl-smiles.1 :compile-at :definition-time) 
  "Simple one non-carbon atom molecules containing an atom from the
organic subset and the appropriate number of hydrogens."
  (check-smiles-strings
   `(("C" 16 ((:c 1) (:h 4)))
     ("N" 17 ((:h 3) (:n 1)))
     ("O" 18 ((:h 2) (:o 1)))
     ("P" 33 ((:h 3) (:p 1)))
     ("S" 34 ((:h 2) (:s 1)))
     ("Cl" 36 ((:cl 1) (:h 1))))))

;;; Bracketed atoms, including those with hydrogens and charges.
(5am:test
 (chemicl-smiles.2 :compile-at :definition-time)
 (check-smiles-strings
  `(("[Au]" 197 ((:au 1)))
    ("[H+]" 1 ((:h 1)) 1)
    ("[OH-]" 17 ((:h 1) (:o 1)) -1)
    ("[OH3+]" 19 ((:h 3) (:o 1)) 1)
    ("[Fe+2]" 56 ((:fe 1)) 2)
    ("[Fe+++]" 56 ((:fe 1)) 3)
    ("[NH4+]" 18 ((:h 4) (:n 1)) 1))))

;;; Simple explicitly (and implicitly) specified bonds.
(5am:test
 (chemicl-smiles.3 :compile-at :definition-time)
 (check-smiles-strings
  `(("CC" 30 ((:c 2) (:h 6)))
    ("C=C" 28 ((:c 2) (:h 4)))
    ("COC" 46 ((:c 2) (:h 6) (:o 1)))
    ("CCO" 46 ((:c 2) (:h 6) (:o 1)))
    ("C=O" 32 ((:c 1) (:h 2) (:o 1)))
    ("O=C=O" 44 ((:c 1) (:o 2)))
    ("C#N" 27 ((:c 1) (:h 1) (:n 1)))
    ("[H][H]" 2 ((:h 2))))))

;;; three ways to say 6-hydroxy-1,4-hexadiene
(5am:test
 (chemicl-smiles.4 :compile-at :definition-time)
 (check-smiles-strings
  `(("C=CCC=CCO" 98 ((:c 6) (:h 10) (:o 1)))
    ("C=C-C-C=CCO" 98 ((:c 6) (:h 10) (:o 1)))
    ("OCC=CCC=C" 98 ((:c 6) (:h 10) (:o 1))))))

;;; branches
(5am:test
 (chemicl-smiles.5 :compile-at :definition-time)
 (check-smiles-strings
  `(("CCN(CC)CC" 101 ((:c 6) (:h 15) (:n 1)))
    ("CC(C)C(=O)O" 88 ((:c 4) (:h 8) (:o 2)))
    ("C=CC(CCC)C(C(C)C)C" 154 ((:c 11) (:h 22))))))

;;; cycles

(5am:test
 (chemicl-smiles.6 :compile-at :definition-time)
 (check-smiles-strings
  `(("C1CCCCC1" 84 ((:c 6) (:h 12))))))

(5am:test
 (chemicl-smiles.7 :compile-at :definition-time)
 (check-smiles-strings
  `(("CC1=CC(Br)CCC1" 175 ((:c 7) (:h 11) (:br 1)))
    ("CC1=CC(CCC1)Br" 175 ((:c 7) (:h 11) (:br 1))))))

(5am:test
 (chemicl-smiles.8 :compile-at :definition-time)
 (check-smiles-strings
  `(("c1(Cl)c(O)cc(Cl)c(Cl)c1" 197 ((:c 6) (:cl 3) (:h 3) (:o 1)))
    ("Clc1cc(O)c(Cl)cc1Cl" 197 ((:c 6) (:cl 3) (:h 3) (:o 1))))))

;;; anthracene
(5am:test
 (chemicl-smiles.9 :compile-at :definition-time)
 (check-smiles-strings
  `(("c1cccc2cc3ccccc3cc12" 178 ((:c 14)  (:h 10))))))

;;; cubane
(5am:test
 (chemicl-smiles.10 :compile-at :definition-time)
 (check-smiles-strings
  `(("C12C3C4C1C5C4C3C25" 104 ((:c 8) (:h 8))))))

(5am:test
 (chemicl-smiles.11 :compile-at :definition-time)
 (check-smiles-strings
  `(("O1CCCCC1N1CCCCC1" 169 ((:c 10) (:h 19) (:n 1) (:o 1))))))

;;; sodium phenoxide salt
(5am:test
 (chemicl-smiles.12 :compile-at :definition-time)
 (let ((salt (parse-smiles-string "[Na+].[O-]c1ccccc1")))
   (5am:is
    (= (length (graph:find-connected-components salt))
       2))
   (5am:is
    (= (reduce #'+
               (mapcar #'charge (graph:find-connected-components salt)))
       0))
   (5am:is
    (= (reduce #'+
               (mapcar
                #'abs
                (mapcar #'charge (graph:find-connected-components salt))))
       2))))

;;; sodium phenoxide salt with ring closure on both sides of the
;;; separated ions
(5am:test
 (chemicl-smiles.13 :compile-at :definition-time)
 (let ((salt (parse-smiles-string "c1cc([O-].[Na+])ccc1")))
   (5am:is
    (= (length (graph:find-connected-components salt))
       2))
   (5am:is
    (= (reduce #'+
               (mapcar #'charge (graph:find-connected-components salt)))
       0))
   (5am:is
    (= (reduce #'+
               (mapcar
                #'abs
                (mapcar #'charge (graph:find-connected-components salt))))
       2))))

;;; sulfuric acid
(5am:test
 (chemicl-smiles.14 :compile-at :definition-time)
 (check-smiles-strings
  `(("OS(=O)(=O)O" 98 ((:h 2) (:o 4) (:s 1))))))

(defparameter *tamoxifen-smiles-string*
  "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3")

(defparameter *smiles-test-strings*
  `(("CC(C)=O" . "CC(C)=O")
    ("CC(=O)C" . "CC(C)=O")
    (,*tamoxifen-smiles-string*
     . ,*tamoxifen-smiles-string*)))

(5am:test
 (chemicl-smiles.15 :compile-at :definition-time)
 (map nil (lambda (x)
            (destructuring-bind (smiles . result)
                x
              (5am:is
               (equal
                (write-smiles-string
                 (parse-smiles-string
                  smiles))
                result))))
      *smiles-test-strings*))
