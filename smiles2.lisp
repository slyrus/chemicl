
(defpackage #:chemicl-smiles
  (:use #:cl #:chemicl #:parser-combinators))

(in-package #:chemicl-smiles)

(defun <aliphatic-chlorine> ()
  (named-seq* (char? #\C) (char? #\l)
              "Cl"))

(defun <aliphatic-bromine> ()
  (named-seq* (char? #\B) (char? #\r)
              "Br"))

(defun <aliphatic-boron> () (named-seq* (char? #\B) "B"))
(defun <aliphatic-carbon> () (named-seq* (char? #\C) "C"))
(defun <aliphatic-nitrogen> () (named-seq* (char? #\N) "N"))
(defun <aliphatic-oxygen> () (named-seq* (char? #\O) "O"))
(defun <aliphatic-sulfur> () (named-seq* (char? #\S) "S"))
(defun <aliphatic-phosphorus> () (named-seq* (char? #\P) "P"))
(defun <aliphatic-fluorine> () (named-seq* (char? #\F) "F"))
(defun <aliphatic-iodine> () (named-seq* (char? #\I) "I"))

(defun <aliphatic-organic-atom> ()
  (named-seq* (<- element (choices1 (<aliphatic-chlorine>)
                                    (<aliphatic-boron>)
                                    (<aliphatic-carbon>)
                                    (<aliphatic-nitrogen>)
                                    (<aliphatic-oxygen>)
                                    (<aliphatic-sulfur>)
                                    (<aliphatic-phosphorus>)
                                    (<aliphatic-fluorine>)
                                    (<aliphatic-iodine>)))
              (make-atom element)))

(defun <aromatic-boron> () (named-seq* (char? #\b) "B"))
(defun <aromatic-carbon> () (named-seq* (char? #\c) "C"))
(defun <aromatic-nitrogen> () (named-seq* (char? #\n) "N"))
(defun <aromatic-oxygen> () (named-seq* (char? #\o) "O"))
(defun <aromatic-sulfur> () (named-seq* (char? #\s) "S"))
(defun <aromatic-phosphorus> () (named-seq* (char? #\p) "P"))

(defun <aromatic-organic-atom> ()
  (named-seq* (<- element (choices1 (<aromatic-boron>)
                                    (<aromatic-carbon>)
                                    (<aromatic-nitrogen>)
                                    (<aromatic-oxygen>)
                                    (<aromatic-sulfur>)
                                    (<aromatic-phosphorus>)))
              (make-atom element)))

(defun <atom> ()
  (choice (<aliphatic-organic-atom>)
          (<aromatic-organic-atom>)))

(defun <atom-seq> ()
  (between* (<atom>) nil nil))

