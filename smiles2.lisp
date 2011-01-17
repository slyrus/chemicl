
(defpackage #:chemicl-smiles
  (:use #:cl #:chemicl #:parser-combinators))

(in-package #:chemicl-smiles)

;;;
;;; special variables to hold state during the parsing and building of molecules/atoms
(defvar *current-molecule*)
(defvar *last-atom*)
(defvar *bond-order*)
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
  (named-seq* (<- pre (upper?))
              (<- suff (atmost? (lower?) 2))
              (format nil "~A~{~A~}" pre suff)))

(defun <bracket-aromatic-atom-symbol> ()
  (apply #'choices1
         (named-seq* (string? "se") "Se")
         (named-seq* (string? "as") "As")
         (map 'list
              #'<aromatic-atom-matcher>
              '("c" "n" "o" "p" "s"))))

(defun <bracket-atom> ()
  (named-seq* #\[ 
              (<- atm (choice1
                       (<bracket-aliphatic-atom-symbol>)
                       (<bracket-aromatic-atom-symbol>)))
              (<- charge (opt? (<charge>)))
              #\]
              (apply #'make-atom atm
                     (when charge `(:charge ,charge)))))

;;;
;;; Atoms
(defun <atom> ()
  (hook?
   (lambda (atom) (setf *last-atom* atom))
   (choices1 (<bracket-atom>)
             (<aliphatic-organic-atom>)
             (<aromatic-organic-atom>))))
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
  (named-seq*
   (<- order (choices1
              (named-seq* (<single-bond>) 1)
              (named-seq* (<double-bond>) 2)
              (named-seq* (<triple-bond>) 3)
              (named-seq* (<quadruple-bond>) 4)
              (named-seq* (<aromatic-bond>) :aromatic)
              (named-seq* (<up-bond>) :up)
              (named-seq* (<down-bond>) :down)))
   (setf *bond-order* order)))

(defun <branch> ()
  (bracket? #\( (delayed? (<chain>)) #\) ))

(defun as-digit? ()
  (hook? #'digit-char-p (digit?)))

(defun <ring-bond> ()
  (choice1
   (named-seq* (opt? (<bond>)) 
               #\%
               (<- digit1 (as-digit?))
               (<- digit2 (as-digit?))
               (+ (* 10 digit1) digit2))
   (named-seq* (opt? (<bond>))
               (<- digit1 (as-digit?))
               digit1)))

;;; use order 0 for disconnected atoms!
(defun <dot> ()
  (named-seq* #\.
              (setf *bond-order* 0)))

(defun <bond-or-dot> ()
  (choice1 (<bond>)
           (<dot>)))

(defmacro always (p)
  `(named-seq* (result nil)
               ,p))

(defun <atom-expr> ()
  (named-seq* (always (print *last-atom*))
              (<- atom (<atom>))
              (<- branches1 (many? (<branch>)))
              (always (print *bond-order*))
              (<- rings (many? (<ring-bond>)))
              (always (setf *bond-order* nil))
              (<- branches2 (many? (<branch>)))
              (opt? (<bond-or-dot>))
              (list atom branches1 rings branches2)))

(defun <chain> ()
  (many? (<atom-expr>)))

(defun parse-smiles-string (str)
  (let (*current-molecule*
        *last-atom*
        *bond-order*
        *atom-counts*
        *pending-rings*
        *open-rings*
        *aromatic*
        *aromatic-atoms*
        *direction*
        *configurations*)
    (parse-string* (<chain>) str :complete t)))
