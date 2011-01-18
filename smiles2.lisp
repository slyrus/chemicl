
(defpackage #:chemicl-smiles
  (:use #:cl #:chemicl #:parser-combinators))

(in-package #:chemicl-smiles)

;;;
;;; parser-combinator utilities
(defun as-digit? ()
  (hook? #'digit-char-p (digit?)))


;;;
;;; Aliphatic Organic Subset Atoms (Cl Br B C N O S P F I)
(defun <aliphatic-organic-atom> ()
  (hook? #'make-atom
         (apply #'choices1
                (map 'list
                     #'string?
                     '("Cl" "Br" "B" "C" "N" "O" "S" "P" "F" "I")))))

;;;
;;; Aromatic Organic Subset Atoms (b c n o s p)
(defun <aromatic-atom-matcher> (str)
  (hook? #'string-upcase (string? str)))

(defun <aromatic-organic-atom> ()
  (hook? #'make-atom
         (apply #'choices1
                (map 'list
                     #'<aromatic-atom-matcher>
                     '("b" "c" "n" "o" "s" "p")))))

;;;
;;; Bracketed atoms e.g. [Na]
(defun <isotope> ()
  (nat?))

(defun <hydrogen-count> ()
  (opt?
   (named-seq*
    #\H
    (<- num (nat?))
    (or num 1))))

(defun <charge> ()
  (choice1 (named-seq*
            (char? #\+)
            (<- charge (opt? (choices1 (named-seq* (char? #\+) 2)
                                       (nat*))))
            (or charge 1))
           (named-seq*
            (char? #\-)
            (<- charge (opt? (choices1 (named-seq* (char? #\-) 2)
                                       (nat*))))
            (- (or charge 1)))))

(defun <bracket-aliphatic-atom-symbol> ()
  (named-seq? (<- pre (upper?))
              (<- suff (atmost? (lower?) 2))
              (format nil "~A~{~A~}" pre suff)))

(defun <bracket-aromatic-atom-symbol> ()
  (apply #'choices1
         (named-seq? (string? "se") "Se")
         (named-seq? (string? "as") "As")
         (map 'list
              #'<aromatic-atom-matcher>
              '("c" "n" "o" "p" "s"))))

(defun <bracket-atom> ()
  (named-seq? #\[ 
              (<- atm (choice1
                       (<bracket-aliphatic-atom-symbol>)
                       (<bracket-aromatic-atom-symbol>)))
              (<- hydrogen-count (<hydrogen-count>))
              (<- charge (opt? (<charge>)))
              #\]
              (let ((atom (apply #'make-atom atm
                                 (append
                                  (when charge `(:charge ,charge))))))
                (when hydrogen-count
                  (print (list 'hydrogens hydrogen-count)))
                atom)))

;;;
;;; Atoms
(defun <atom> ()
  (choices1 (<bracket-atom>)
            (<aliphatic-organic-atom>)
            (<aromatic-organic-atom>)))
;;;
;;; Bonds
(defun <single-bond> () (char? #\-))
(defun <double-bond> () (char? #\=))
(defun <triple-bond> () (char? #\#))
(defun <quadruple-bond> () (char? #\$))
(defun <aromatic-bond> () (char? #\:))
(defun <up-bond> () (char? #\/))
(defun <down-bond> () (char? #\\))

(defun <bond> ()
  (choices
   (named-seq? (<single-bond>) 1)
   (named-seq? (<double-bond>) 2)
   (named-seq? (<triple-bond>) 3)
   (named-seq? (<quadruple-bond>) 4)
   (named-seq? (<aromatic-bond>) :aromatic)
   (named-seq? (<up-bond>) :up)
   (named-seq? (<down-bond>) :down)))

(defun <bond-or-dot> ()
  (choice1 (<bond>)
           (<dot>)))

;;;
;;; Core atom rule
;;;
;;; This rule encompasses the atom itself, all connected branches,
;;; ring-bond symbols, and a bond or dot symbol indicating arity (or
;;; disconnectedness) of the next (non-chain) bond from this atom
(defun <atom-expr> (mol)
  (mdo (<- atom (<atom>))
       (<- branches1 (many? (<branch> mol atom)))
       (<- rings (many? (<ring-bond>)))
       (<- branches2 (many? (<branch> mol atom)))
       (opt? (<bond-or-dot>))
       (result (list atom branches1 rings branches2))))


;;;
;;; Branches
(defun <ring-bond> ()
  (choice
   (named-seq? (opt? (<bond>)) 
               #\%
               (<- digit1 (as-digit?))
               (<- digit2 (as-digit?))
               (+ (* 10 digit1) digit2))
   (named-seq? (opt? (<bond>))
               (<- digit1 (as-digit?))
               digit1)))

(defun <branched-atom> (mol)
  (mdo 
    (<- bond (opt? (<bond>)))
    (<- atom (<atom>))
    (<- ring-bonds (many? (<ring-bond>)))
    (<- branches (many? (<branch> mol atom)))
       (result (list bond atom ring-bonds branches))))
  
(defun <branch> (mol atom)
  (bracket? #\(
            (choices
             (<chain> mol)
             (named-seq?
              (<- bond (<bond>))
              (<- chain (curtail? (<chain> mol)))
              (list bond chain))
             (named-seq?
              (<- dot (<dot>))
              (<- chain (curtail? (<chain> mol)))
              (list dot chain)))
            #\)))

;;; use order 0 for disconnected atoms!
(defun <dot> () (char? #\.) )

(defun unzip (l)
  (if (and (car l) (listp (car l)))
      (cons (cdr l) (unzip (car l)))
      (cons (cdr l) (list (car l)))))
;;;
;;; Core chain rule -- a (possibly branched) chain of atoms
(defun <chain> (mol)
  (hook? #'unzip (chainl1?
                  (<branched-atom> mol)
                  (result #'cons))))

;;;
;;; Main (only?) (useful) entry point to all of this code
;;; 
;;; Parsess a SMILES string and returns a molecule (should this really
;;; return a set of molecules?)
(defun parse-smiles-string (str)
  (let ((mol (make-molecule)))
    (parse-string* (<chain> mol) str :complete t)))
