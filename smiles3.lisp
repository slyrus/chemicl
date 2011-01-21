
(cl:defpackage smiles-minimal
  (:use #:cl #:parser-combinators #:chemicl))

(in-package :smiles-minimal)

(defclass the-atom ()
  ((name :accessor name-of :initarg :name))
  (:documentation "Immutable atom instance"))

(defclass the-edge ()
  ((source :accessor source-of :initarg :source)
   (dest   :accessor dest-of   :initarg :dest)
   (kind   :accessor kind-of   :initarg :kind))
  (:documentation "Immutable edge instance"))

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
  (fset:map (:root-atom (fset:@ root-atom :root-atom))
            (:tip-atom (fset:@ root-atom :tip-atom))
            (:edge-set (fset:set (fset:$ (fset:@ root-atom :edge-set))
                                 (fset:$ (fset:@ branch :edge-set))
                                 (make-instance 'the-edge
                                                :source (fset:@ branch :root-atom)
                                                :dest (fset:@ root-atom :root-atom)
                                                :kind (fset:@ branch :root-bond))))
            (:atom-set (fset:union (fset:@ root-atom :atom-set)
                                   (fset:@ branch :atom-set)))
            (:ring-tags (fset:@ root-atom :ring-tags))))

(defun <branch> (subchain-parser)
  (bracket? #\(
            (named-seq? (<- bond (choices #\- #\: #\= #\# #\\ #\/ #\$  (result nil)))
                        (<- subchain subchain-parser)
                        (fset:map (fset:$ subchain) (:root-bond bond)))
            #\)))

(defun <atom-with-branches> (subchain-parser)
  (named-seq? (<- atom (<atom>))
              (<- ring-tags (many? (digit?)))
              (<- branches (many? (<branch> subchain-parser)))
              (reduce #'branch-merge branches
                      :initial-value (make-root-atom-structure atom ring-tags))))

(defun handle-ring-tags (atom1 atom2)
  (let ((tags1 (fset:@ atom1 :ring-tags))
        (tags2 (fset:@ atom2 :ring-tags)))
    (let* ((edges (fset:map-intersection tags1 tags2
                                         #'(lambda (source dest)
                                             (make-instance 'the-edge
                                                            :source source
                                                            :dest dest
                                                            :kind nil))))
           (tags (fset:set-difference (fset:union (fset:domain tags1)
                                                  (fset:domain tags2))
                                      (fset:domain edges))))
      (values (fset:range edges)
              (gmap:gmap :map #'(lambda (tag)
                                  (values tag (or (fset:@ tags1 tag)
                                                  (fset:@ tags2 tag))))
                         (:set tags))))))

(defun <bond-function> ()
  (named-seq? (<- bond (choices #\- #\: #\= #\# (result nil)))
              #'(lambda (atom1 atom2)
                  (multiple-value-bind (new-edges new-ring-tags)
                      (handle-ring-tags atom1 atom2)
                    (fset:map (:root-atom (fset:@ atom1 :root-atom))
                              (:tip-atom (fset:@ atom2 :tip-atom))
                              (:edge-set
                               (fset:set (fset:$ (fset:@ atom1 :edge-set))
                                         (fset:$ (fset:@ atom2 :edge-set))
                                         (fset:$ new-edges)
                                         (make-instance 'the-edge
                                                        :source (fset:@ atom2 :root-atom)
                                                        :dest (fset:@ atom1 :tip-atom)
                                                        :kind bond)))
                              (:atom-set (fset:union (fset:@ atom1 :atom-set)
                                                     (fset:@ atom2 :atom-set)))
                              (:ring-tags new-ring-tags))))))

(defun <chain> ()
  (named? chain
    (chainl1? (<atom-with-branches> chain)
              (<bond-function>))))

(defun emit-dot (molecule)
  (let ((atom-set (fset:@ molecule :atom-set))
        (edge-set (fset:@ molecule :edge-set))
        (atom-number-map (fset:empty-map))
        (counter 0))
    (fset:do-set (atom atom-set)
      (setf atom-number-map
            (fset:with atom-number-map atom (incf counter))))
    (with-output-to-string (str)
      (format str "graph smiles {~&")
      (fset:do-set (edge edge-set)
        (format str "\"~a ~a\" -- \"~a ~a\" [label = \"~a\"]~&"
                (id (element (source-of edge)))
                (fset:@ atom-number-map (source-of edge))
                (id (element (dest-of edge)))
                (fset:@ atom-number-map (dest-of edge))
                (kind-of edge)))
      (format str "}"))))

(defun parse-smiles-string (str)
  (parse-string* (<chain>) str))
