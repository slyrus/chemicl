
(cl:eval-when (:compile-toplevel :load-toplevel :execute)
 (ql:quickload 'named-readtables))

(cl:defpackage #:smiles3
  (:use #:cl #:parser-combinators #:chem)
  (:shadowing-import-from #:cl #:atom)
  (:shadow #:parse-smiles-string)
  (:import-from #:chemicl
                #:molecule
                #:get-element
                #:make-atom
                #:element
                #:add-atom
                #:id
                #:atoms
                #:atom1
                #:atom2
                #:bond
                #:bond-order))

(in-package :smiles3)

(defclass smiles-atom (chemicl:atom)
  ((aromatic :initarg :aromatic :accessor aromatic :initform nil)
   (explicit-hydrogen-count :initarg :explicit-hydrogen-count
                            :accessor explicit-hydrogen-count
                            :initform nil)
   (orientation :initarg :orientation :accessor orientation :initform nil)))

(defun make-smiles-atom (element-identifier &rest args)
  (apply #'make-instance
         'smiles-atom
         :element (get-element element-identifier)
         args))

;;;
;;; parser-combinator utilities
(defun as-digit? ()
  (hook? #'digit-char-p (digit?)))

;;;
;;; Aliphatic Organic Subset Atoms (Cl Br B C N O S P F I)
(defun <aliphatic-organic-atom> ()
  (hook? #'make-smiles-atom
         (apply #'choices1
                (map 'list
                     #'string?
                     '("Cl" "Br" "B" "C" "N" "O" "S" "P" "F" "I")))))

;;;
;;; Aromatic Organic Subset Atoms (b c n o s p)
(defun <aromatic-atom-matcher> (str)
  (hook? #'string-upcase (string? str)))

(defun <aromatic-organic-atom> ()
  (hook? (lambda (element)
           (make-smiles-atom element :aromatic t))
         (apply #'choices1
                (map 'list
                     #'<aromatic-atom-matcher>
                     '("b" "c" "n" "o" "s" "p")))))

;;;
;;; Bracketed atoms 
;;; e.g. [H] [Na] [13C] [O-]
;;; (not yet working: [C@] [C@@]
(defun <isotope> ()
  (nat?))

(defun <hydrogen-count> ()
  (named-seq?
   #\H
   (<- num (opt? (nat?)))
   (or num 1)))

(defun <charge> ()
  (choice1 (named-seq?
            (char? #\+)
            (<- charge (opt? (choices1 (named-seq? (char? #\+) 2)
                                       (nat?))))
            (or charge 1))
           (named-seq?
            (char? #\-)
            (<- charge (opt? (choices1 (named-seq? (char? #\-) 2)
                                       (nat?))))
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

(defun <orientation> ()
  (choices1 (chook? :clockwise "@@")
            (chook? :counterclockwise "@")))

(defun <bracket-modifier> ()
  (choices1 (named-seq? (<- hydrogens (<hydrogen-count>))
                        (list :hydrogen-count hydrogens))
            (named-seq? (<- charge (<charge>))
                        (list :charge charge))
            (named-seq? (<- orientation (<orientation>))
                        (list :orientation orientation))))

(defun <bracket-atom> ()
  (named-seq? #\[ 
              (<- isotope-number (opt? (nat?)))
              (<- atom-aromaticity (choice1
                                    (seq-list?
                                     (<bracket-aliphatic-atom-symbol>)
                                     (result nil))
                                    (seq-list?
                                     (<bracket-aromatic-atom-symbol>)
                                     (result t))))
              (<- mods (validate?
                        (many? (<bracket-modifier>))
                        (lambda (x)
                          (if (<= (apply #'max
                                         (map 'list
                                              #'(lambda (f) (funcall f x :key #'car))
                                              (map 'list 
                                                   #'(lambda (y) (alexandria:curry #'count y))
                                                   '(:hydrogen-count :charge :orientation))))
                                  1)
                              t
                              (error "Bracket modifier appeared more than once!")))))
              #\]
              (destructuring-bind (element aromaticity)
                  atom-aromaticity
                (destructuring-bind (&key charge hydrogen-count orientation)
                    (apply #'append mods)
                  (let* ((isotope
                          (when isotope-number
                            (chem::get-isotope element isotope-number)))
                         (atom (apply #'make-smiles-atom element
                                      (append
                                       (when isotope `(:isotope ,isotope))
                                       (when charge `(:charge ,charge))
                                       (when aromaticity `(:aromatic ,aromaticity))
                                       (when hydrogen-count
                                         `(:explicit-hydrogen-count ,hydrogen-count))
                                       (when orientation `(:orientation ,orientation))))))
                    atom)))))

;;;
;;; Atoms
(defun <atom> ()
  (choices1 (<bracket-atom>)
            (<aliphatic-organic-atom>)
            (<aromatic-organic-atom>)))
;;;
;;; Bonds
(defun <single-bond> () #\-)
(defun <double-bond> () #\=)
(defun <triple-bond> () #\#)
(defun <quadruple-bond> () #\$)
(defun <aromatic-bond> () #\:)
(defun <up-bond> () #\/)
(defun <down-bond> () #\\)

(defun <bond> ()
  (choices
   (chook? 1 (<single-bond>))
   (chook? 2 (<double-bond>))
   (chook? 3 (<triple-bond>))
   (chook? 4 (<quadruple-bond>))
   (chook? :aromatic (<aromatic-bond>))
   (chook? :up (<up-bond>))
   (chook? :down (<down-bond>))
   (result nil)))

(defun collect-ring-tags (atom ring-tags map)
  (if (null ring-tags)
      map
      (collect-ring-tags atom
                         (cdr ring-tags)
                         (fset:with map (car ring-tags) atom))))

(defun make-root-atom-structure (root-atom ring-tags)
  (fset:map (:root-atom root-atom)
            (:tip-atom root-atom)
            (:edge-set (fset:empty-set))
            (:atom-set (fset:set root-atom))
            (:ring-tags (collect-ring-tags root-atom ring-tags (fset:empty-map)))))

(defun branch-merge (root-atom branch)
  (multiple-value-bind (new-edges new-ring-tags)
      (handle-ring-tags root-atom branch)
    (fset:map (:root-atom (fset:@ root-atom :root-atom))
              (:tip-atom (fset:@ root-atom :tip-atom))
              (:edge-set (fset:set (fset:$ (fset:@ root-atom :edge-set))
                                   (fset:$ (fset:@ branch :edge-set))
                                   (fset:$ new-edges)
                                   (make-instance 'bond
                                                  :atom1 (fset:@ root-atom :tip-atom)
                                                  :atom2 (fset:@ branch :root-atom)
                                                  :order (fset:@ branch :root-bond))))
              (:atom-set (fset:union (fset:@ root-atom :atom-set)
                                     (fset:@ branch :atom-set)))
              (:ring-tags new-ring-tags))))

(defun <branch> (subchain-parser)
  (bracket? #\(
            (named-seq? (<- bond (<bond>))
                        (<- subchain subchain-parser)
                        (fset:map (fset:$ subchain) (:root-bond bond)))
            #\)))

(defun <ring-tag> ()
  (choice
   (named-seq? (opt? (<bond>)) 
               #\%
               (<- digit1 (as-digit?))
               (<- digit2 (as-digit?))
               (+ (* 10 digit1) digit2))
   (named-seq? (opt? (<bond>))
               (<- digit1 (as-digit?))
               digit1)))

(defun <atom-with-branches> (subchain-parser)
  (named-seq? (<- atom (<atom>))
              (<- ring-tags1 (many? (<ring-tag>)))
              (<- branches (many? (<branch> subchain-parser)))
              (<- ring-tags2 (many? (<ring-tag>)))
              (reduce #'branch-merge branches
                      :initial-value
                      (make-root-atom-structure atom (append ring-tags1 ring-tags2)))))

(defun handle-ring-tags (atom1 atom2)
  (let ((tags1 (fset:@ atom1 :ring-tags))
        (tags2 (fset:@ atom2 :ring-tags)))
    (let* ((edges (fset:map-intersection tags1 tags2
                                         #'(lambda (source dest)
                                             (make-instance 'bond
                                                            :atom1 source
                                                            :atom2 dest
                                                            :order 1))))
           (tags (fset:set-difference (fset:union (fset:domain tags1)
                                                  (fset:domain tags2))
                                      (fset:domain edges))))
      (values (fset:range edges)
              (gmap:gmap :map #'(lambda (tag)
                                  (values tag (or (fset:@ tags1 tag)
                                                  (fset:@ tags2 tag))))
                         (:set tags))))))

(defun <bond-function> ()
  (named-seq? (<- bond (<bond>))
              #'(lambda (atom1 atom2)
                  (multiple-value-bind (new-edges new-ring-tags)
                      (handle-ring-tags atom1 atom2)
                    (fset:map (:root-atom (fset:@ atom1 :root-atom))
                              (:tip-atom (fset:@ atom2 :tip-atom))
                              (:edge-set
                               (fset:set (fset:$ (fset:@ atom1 :edge-set))
                                         (fset:$ (fset:@ atom2 :edge-set))
                                         (fset:$ new-edges)
                                         (make-instance 'bond
                                                        :atom1 (fset:@ atom1 :tip-atom)
                                                        :atom2 (fset:@ atom2 :tip-atom)
                                                        :order bond)))
                              (:atom-set (fset:union (fset:@ atom1 :atom-set)
                                                     (fset:@ atom2 :atom-set)))
                              (:ring-tags new-ring-tags))))))

(defun <chain> ()
  (named? chain
    (chainl1? (<atom-with-branches> chain)
              (<bond-function>))))


(define-condition smiles-parsing-error (error)
  ((smiles-string
    :initarg :smiles-string
    :reader smiles-string))
  (:report (lambda (condition stream)
             (format stream "Error parsing SMILES string ~S"
                     (smiles-string condition)))))

(defun smiles-reader-error (stream control &rest args)
  (error 'reader-error
         :stream stream
         :format-control control
         :format-arguments args))

(defun parse-smiles-string (str)
  (let ((parsed (parse-string* (<chain>) str :complete t)))
    (if parsed
      (let ((element-counts (make-hash-table)))
        (flet ((element-number (atom)
                 (setf (gethash (element atom) element-counts)
                       (let ((cur (gethash (element atom) element-counts)))
                         (if cur (1+ cur) 1)))))
          (let ((mol (make-instance 'molecule)))
            (fset:image (lambda (atom)
                          (add-atom mol atom
                                    (format nil "~A~A"
                                            (id (element atom))
                                            (element-number atom))))
                        (fset:@ parsed :atom-set))
            (fset:image (lambda (edge) (epigraph:add-edge mol edge))
                        (fset:@ parsed :edge-set))
            ;; further postprocessing needed here:
            ;; 1. fixup aromatic bonds
            ;; 2. add implicit hydrogens
            mol)))
      (error 'smiles-parsing-error :smiles-string str))))

(defun emit-dot (molecule)
  (let ((atom-number-map (fset:empty-map))
        (counter 0))
    (loop for atom in (atoms molecule)
       do (setf atom-number-map
            (fset:with atom-number-map atom (incf counter))))
    (with-output-to-string (str)
      (format str "graph smiles {~&")
      (loop for edge in (chemicl::bonds molecule)
         do (format str "\"~a ~a\" -- \"~a ~a\" [label = \"~a\"]~&"
                    (id (element (chemicl::atom1 edge)))
                    (fset:@ atom-number-map (chemicl::atom1 edge))
                    (id (element (chemicl::atom2 edge)))
                    (fset:@ atom-number-map (chemicl::atom2 edge))
                   (bond-order edge)))
      (format str "}"))))


