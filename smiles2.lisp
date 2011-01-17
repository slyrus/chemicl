
(defpackage #:chemicl-smiles
  (:use #:cl #:chemicl #:parser-combinators))

(in-package #:chemicl-smiles)

;;;
;;; special variables to hold state during the parsing and building of molecules/atoms
(defvar *current-molecule*)
(defvar *last-atom*)
(defvar *atom-counts*)
(defvar *pending-rings*)
(defvar *open-rings*)
(defvar *aromatic*)
(defvar *aromatic-atoms*)
(defvar *direction*)
(defvar *configurations*)

;;;
;;; Aliphatic Organic Subset Atoms (Cl Br B C N O S P F I)
(defun <aliphatic-organic-atom> ()
  (hook? #'make-atom
         (apply #'choices
                (map 'list
                     #'string?
                     '("Cl" "Br" "B" "C" "N" "O" "S" "P" "F" "I")))))

;;;
;;; Aromatic Organic Subset Atoms (b c n o s p)
(defun <aromatic-atom-matcher> (str)
  (hook? #'string-upcase (string? str)))

(defun <aromatic-organic-atom> ()
  (hook? #'make-atom
         (apply #'choices 
                (map 'list
                     #'<aromatic-atom-matcher>
                     '("b" "c" "n" "o" "s" "p")))))

;;;
;;; Bracketed atoms e.g. [Na]
(defun <charge> ()
  (choice (named-seq*
           (char? #\+)
           (<- charge (opt? (choices (named-seq* (char? #\+) 2)
                                     (nat*))))
           (or charge 1))
          (named-seq*
           (char? #\-)
           (<- charge (opt? (choices (named-seq* (char? #\-) 2)
                                     (nat*))))
           (- (or charge 1)))))

(defun <bracket-aliphatic-atom-symbol> ()
  (named-seq* (<- pre (upper?))
              (<- suff (atmost? (lower?) 2))
              (format nil "~A~{~A~}" pre suff)))

(defun <bracket-aromatic-atom-symbol> ()
  (apply #'choices
         (named-seq* (string? "se") "Se")
         (named-seq* (string? "as") "As")
         (map 'list
              #'<aromatic-atom-matcher>
              '("c" "n" "o" "p" "s"))))

(defun <bracket-atom> ()
  (named-seq* #\[ 
              (<- atm (choice
                       (<bracket-aliphatic-atom-symbol>)
                       (<bracket-aromatic-atom-symbol>)))
              (<- charge (opt? (<charge>)))
              #\]
              (apply #'make-atom atm
                     (when charge `(:charge ,charge)))))

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

(defun <branch> ()
  (named-seq* #\(
              #\)
              ))

(defun <atom-seq> ()
  (many? (named-seq*
          (<- atom (<atom>))
          (opt? (<branch>))
          (setf *last-atom* atom))))

(defun parse-smiles-string (str)
  (let (*current-molecule*
        *last-atom*
        *atom-counts*
        *pending-rings*
        *open-rings*
        *aromatic*
        *aromatic-atoms*
        *direction*
        *configurations*)
    (parse-string* (<atom-seq>) str :complete t)))
