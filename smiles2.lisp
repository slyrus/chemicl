
(defpackage #:chemicl-smiles
  (:use #:cl #:chemicl #:parser-combinators))

(in-package #:chemicl-smiles)

;;;
;;; Aliphatic Organic Subset Atoms (Cl Br B C N O S P F I)
(defun <aliphatic-chlorine> () (string? "Cl"))
(defun <aliphatic-bromine> () (string? "Br"))
(defun <aliphatic-boron> () (string? "B"))
(defun <aliphatic-carbon> () (string? "C"))
(defun <aliphatic-nitrogen> () (string? "N"))
(defun <aliphatic-oxygen> () (string? "O"))
(defun <aliphatic-sulfur> () (string? "S"))
(defun <aliphatic-phosphorus> () (string? "P"))
(defun <aliphatic-fluorine> () (string? "F"))
(defun <aliphatic-iodine> () (string? "I"))

(defun <aliphatic-organic-atom> ()
  (named-seq* (<- element (choices (<aliphatic-chlorine>)
                                    (<aliphatic-boron>)
                                    (<aliphatic-carbon>)
                                    (<aliphatic-nitrogen>)
                                    (<aliphatic-oxygen>)
                                    (<aliphatic-sulfur>)
                                    (<aliphatic-phosphorus>)
                                    (<aliphatic-fluorine>)
                                    (<aliphatic-iodine>)))
              (make-atom element)))

;;;
;;; Aromatic Organic Subset Atoms (b c n o s p)
(defun <aromatic-boron> () (hook? #'string-upcase (string? "b")))
(defun <aromatic-carbon> () (hook? #'string-upcase (string? "c")))
(defun <aromatic-nitrogen> () (hook? #'string-upcase (string? "n")))
(defun <aromatic-oxygen> () (hook? #'string-upcase (string? "o")))
(defun <aromatic-sulfur> () (hook? #'string-upcase (string? "s")))
(defun <aromatic-phosphorus> () (hook? #'string-upcase (string? "p")))

(defun <aromatic-organic-atom> ()
  (named-seq* (<- element (choices (<aromatic-boron>)
                                    (<aromatic-carbon>)
                                    (<aromatic-nitrogen>)
                                    (<aromatic-oxygen>)
                                    (<aromatic-sulfur>)
                                    (<aromatic-phosphorus>)))
              (make-atom element)))

;;;
;;; Bracketed atoms e.g. [Na]

(defun <charge> ()
  (choice (named-seq*
           (char? #\+)
           (<- charge (atmost?
                       (choices (named-seq* (char? #\+) 2)
                                (nat*))
                       1))
           (or (and charge (first charge)) 1))
          (named-seq*
           (char? #\-)
           (<- charge (atmost?
                       (choices (named-seq* (char? #\-) 2)
                                (nat*))
                       1))
           (- (or (and charge (first charge)) 1)))))

(defun <bracket-aliphatic-atom-symbol> ()
  (named-seq* (<- pre (upper?))
              (<- suff (atmost? (lower?) 2))
              (format nil "~A~{~A~}" pre suff)))

(defun <bracket-aromatic-atom-symbol> ()
  (choices (named-seq* (char? #\s) (char? #\e) "Se")
           (named-seq* (char? #\a) (char? #\s) "As")
           (named-seq* (char? #\c) "C")
           (named-seq* (char? #\n) "N")
           (named-seq* (char? #\o) "O")
           (named-seq* (char? #\p) "P")
           (named-seq* (char? #\s) "S")))

(defun <bracket-atom> ()
  (named-seq* #\[ 
              (<- atm (choice
                       (<bracket-aliphatic-atom-symbol>)
                       (<bracket-aromatic-atom-symbol>)))
              (<- charge (atmost? (<charge>) 1))
              #\]
              (apply #'make-atom atm
                     (when charge `(:charge ,(first charge))))))

;;;
;;; Atoms
(defun <atom> ()
  (choices (<bracket-atom>)
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
   (named-seq* (<single-bond>) 1)
   (named-seq* (<double-bond>) 2)
   (named-seq* (<triple-bond>) 3)
   (named-seq* (<quadruple-bond>) 4)
   (named-seq* (<aromatic-bond>) :aromatic)
   (named-seq* (<up-bond>) :up)
   (named-seq* (<down-bond>) :down)))

(defun <atom-seq> ()
  (many? (<atom>)))

