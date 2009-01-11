;;; file: chemicl.lisp
;;;
;;; Copyright (c) 2008-2009 Cyrus Harmon (ch-lisp@bobobeach.com)
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms prohibited.
;;;

(in-package :chemicl)

(defclass atom ()
  ((name :initarg :name :accessor atom-name :initform nil)
   (element :initarg :element :accessor element :initform nil)
   (charge :initarg :charge :accessor charge :initform 0)
   (isotope :initarg :isotope :accessor isotope :initform nil)
   (hybridization :initarg :hybridization
                  :accessor hybridization
                  :initform nil))
  (:documentation "A class for representing individual atoms. For
  example, a molecule of hydrogen class would contain two atom
  instances, each of whose element slots would contain the (same)
  element instance for the hydrogen. An atom can be associated with at
  most one molecule at a time, as sepcified by its molecule slot,
  which can be NIL, indicating that the atom is not associated with
  any molecule."))

(defparameter *atom-print-verbosity* 0)

(defmethod print-object ((object atom) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (format stream " ~S" 
            (if (slot-boundp object 'element)
                (when (element object)
                  (if (> *atom-print-verbosity* 0)
                      (element object)
                      (id (element object))))
                "Element: unbound"))
    (when (slot-boundp object 'name)
      (format stream " ~S" (atom-name object)))
    (cond ((plusp (charge object))
           (format stream " +~S" (charge object)))
          ((minusp (charge object))
           (format stream " ~S" (charge object))))))

(defmethod atomic-number ((atom atom))
  (atomic-number (element atom)))

(defmethod id ((atom atom))
  (id (element atom)))

(defmethod mass ((atom atom))
  (mass (element atom)))

(defmethod exact-mass ((atom atom))
  (if (isotope atom)
      (isotope-exact-mass (isotope atom))
      (isotope-exact-mass (car (isotopes (element atom))))))

(defmethod get-normal-valence ((atom atom))
  (get-normal-valence (element atom)))

(defun make-atom (element-identifier &rest args)
  "Returns a new instance of the atom class as specified by an
element-identifer, which can be a number representing the atomic
number of the element (such as 6 for Carbon), or a string or a lisp
symbol containing an element symbol (such as Fe or :fe for Iron)."
  (apply #'make-instance
         'atom
         :element (get-element element-identifier)
         args))

;;; We need a way of explicitly storing configurations around double
;;; bonds and configurations at chiral centers. Hrm...
(defclass molecule (graph:simple-edge-list-graph)
  ((name :initarg :name :accessor name :initform nil)
   (atom-name-hash :accessor atom-name-hash
                   :initform (make-hash-table :test 'equal)))
  (:documentation "A class for representing molecules."))

(defun make-molecule (&rest args)
  (apply #'make-instance 'molecule args))

(defmethod graph:copy-graph ((graph molecule) &key copy-edges)
  (declare (ignore copy-edges))
  (let ((new-graph (call-next-method)))
    (when (slot-boundp graph 'name)
      (setf (name new-graph) (name graph)))
    (when (slot-boundp graph 'atom-name-hash)
      (setf (atom-name-hash new-graph)
            (alexandria:copy-hash-table (atom-name-hash graph))))
    new-graph))

(defgeneric copy-molecule (molecule)
  (:method ((molecule molecule))
    (graph:copy-graph molecule)))

(defgeneric get-atom (molecule atom-identifier)
  (:method ((molecule molecule) (atom atom))
    atom)
  (:method ((molecule molecule) atom-identifier)
    (gethash atom-identifier (atom-name-hash molecule)))
  (:documentation "If atom-identifier is an atom, returns
  atom-identifier as the atom-object. Otherwise, the atom named by
  atom-identifier is returned."))

(defgeneric add-atom (molecule element-identifier name)
  (:method ((molecule molecule) element-identifier name)
    (let ((atom (make-atom element-identifier :name name)))
      (graph:add-node molecule atom)
      (setf (gethash name (atom-name-hash molecule)) atom)
      atom))
  (:documentation "Adds a new instance of ATOM of the element
  specified by the element-identifier, which can be an atomic number,
  or a string or symbol whose value is the two letter elemental symbol
  of the desired element."))

(defgeneric remove-atom (molecule atom)
  (:method ((molecule molecule) (atom atom))
    (map nil
         (lambda (edge) (graph:remove-edge molecule edge))
         (graph:find-edges-containing molecule atom))
    (graph:remove-node molecule atom)
    (remhash (atom-name atom) (atom-name-hash molecule))
    atom)
  (:documentation "Removes the atom (and any edges containing it) from
  molecule."))

(defgeneric atom-count (molecule)
  (:method ((molecule molecule))
    (graph:node-count molecule))
  (:documentation "Returns the number of atoms in molecule."))

(defun map-atoms (fn molecule)
  (graph:map-nodes fn molecule))

(defun map-atoms->list (fn molecule)
  (graph:map-nodes->list fn molecule))

(defparameter *bond-orders* 
  '((:single . 1)
    (:aromatic . 1.5)
    (:double . 2)
    (:triple . 3)))

(defun get-bond-order-number (keyword-or-number)
  (etypecase keyword-or-number
    (keyword (cdr (assoc keyword-or-number *bond-orders*)))
    (number (cdr (rassoc keyword-or-number *bond-orders*)))))

(defun get-bond-order-keyword (keyword-or-number)
  (etypecase keyword-or-number
    (keyword (car (assoc keyword-or-number *bond-orders*)))
    (number (car (rassoc keyword-or-number *bond-orders*)))))

(defclass bond (graph:edge)
  ((type :accessor bond-type :initarg :type :initform :single)
   (order :accessor bond-order :initarg :order :initform 1)
   (direction :accessor bond-direction :initarg :direction :initform nil)))

(defmethod graph:copy-edge ((object bond))
  (let ((new-edge (call-next-method)))
    (setf (bond-type new-edge) (bond-type object) 
          (bond-order new-edge) (bond-order object) 
          (bond-direction new-edge) (bond-direction object))
    new-edge))

(defmethod graph:print-edge-data :after ((object bond) stream)
  (format stream " ~S ~S" 
          (if (slot-boundp object 'type)
              (bond-type object)
              "Type: unbound")
          (if (slot-boundp object 'order)
              (bond-order object)
              "Order: unbound"))
  (when (and (slot-boundp object 'direction)
           (bond-direction object))
    (format stream " ~S" (bond-direction object))))

(defgeneric atom-bond-order (molecule atom))

(defmethod atom-bond-order ((molecule molecule) (atom atom))
  (let ((bonds (graph:find-edges-containing molecule atom)))
    (reduce #'+ bonds :key 'bond-order)))

(defmethod atom-bond-order ((molecule molecule) atom-identifier)
  (let ((atom (get-atom molecule atom-identifier)))
    (when atom
      (atom-bond-order molecule atom))))

(defgeneric add-bond (molecule atom-identifier-1 atom-identifier-2
                               &key type order direction))

(defgeneric remove-bond (molecule atom-or-bond &rest args))

(defgeneric find-bonds-containing (molecule atom))

(defmethod add-bond ((molecule molecule) atom-identifier-1 atom-identifier-2
                     &key (type :single) (order 1) direction)
  (let ((atom-1 (get-atom molecule atom-identifier-1))
        (atom-2 (get-atom molecule atom-identifier-2)))
    (let ((bond (apply #'make-instance 'bond
                       :node1 atom-1
                       :node2 atom-2
                       (append
                        (when type `(:type ,type))
                        (when order `(:order ,order))
                        (when direction `(:direction ,direction))))))
      (graph:add-edge molecule bond)
      bond)))

(defmethod remove-bond ((molecule molecule) (atom-identifier-1 atom) &rest args)
  (destructuring-bind (atom-identifier-2)
      args
    (graph:remove-edge-between-nodes molecule atom-identifier-1 atom-identifier-2)))

(defmethod remove-bond ((molecule molecule) (bond bond) &rest args)
  (declare (ignore args))
  (graph:remove-edge molecule bond))

(defmethod find-bonds-containing ((molecule molecule) (atom atom))
  (graph:find-edges-containing molecule atom))


(defun double-bond-configuration (molecule bond)
  (let ((atom1 (graph:node1 bond))
        (atom2 (graph:node2 bond)))
    (let ((atom1-bonds
           (remove bond (graph:find-edges-containing molecule atom1)))
          (atom2-bonds
           (remove bond (graph:find-edges-containing molecule atom2))))
      (cons atom1-bonds atom2-bonds))))

(defun remove-keyword-args (keywords list)
  (if (listp keywords)
      (loop for (x y) on list by #'cddr
         append (unless (member x keywords)
                  (list x y)))
      (loop for (x y) on list by #'cddr
         append (unless (eq x keywords)
                  (list x y)))))

(defmethod mass ((molecule molecule))
  (let ((mass 0.0d0))
    (graph:map-nodes (lambda (atom)
                       (incf mass (mass atom)))
                     molecule)
    mass))

(defmethod exact-mass ((molecule molecule))
  (let ((exact-mass 0.0d0))
    (graph:map-nodes (lambda (atom)
                       (incf exact-mass (exact-mass atom)))
                     molecule)
    exact-mass))

(defmethod charge ((molecule molecule))
  (let ((charge 0))
    (graph:map-nodes (lambda (atom)
                       (incf charge (charge atom)))
                     molecule)
    charge))

(defun count-element (molecule element-id)
  (let ((element-count 0)
        (element (get-element element-id)))
    (graph:dfs-map molecule (graph:first-node molecule)
             (lambda (atom)
               (when (eq (element atom) element)
                 (incf element-count))))
    element-count))

(defun count-elements (molecule)
  (let ((element-count-hash (make-hash-table)))
    (graph:dfs-map molecule (graph:first-node molecule)
                   (lambda (atom)
                     (setf (gethash (element atom) element-count-hash)
                           (1+ (or (gethash (element atom)
                                            element-count-hash)
                                   0)))))
    (let (l)
      (maphash (lambda (k v)
                 (push (cons k v) l))
               element-count-hash)
      (sort l #'string< :key (lambda (x) (id (car x)))))))

(defun molecular-formula (molecule)
  (format nil "~{~{~A~A~}~}"
          (mapcar (lambda (x)
                    (list (id (car x))
                          (cdr x)))
                  (count-elements molecule))))

(defun excess-electrons (molecule atom)
  (cond ((eql (element atom) (get-element "C"))
         (if (< (atom-bond-order molecule atom) 4) 1 0))
        ((eql (element atom) (get-element "N"))
         (if (< (atom-bond-order molecule atom) (+ (charge atom) 3)) 1 0))
        ((eql (element atom) (get-element "O"))
         2)
        (t 0)))

(defun aromaticp (molecule ring)
  "Uses Hunckle's 4n+2 rule to decide if a ring in a molecule is
aromatic or not."
  (let ((excess-electrons
         (reduce #'+ (map 'list (lambda (atom) (excess-electrons molecule atom))
                          ring))))
    (and (> excess-electrons 5)
         (zerop (rem (- excess-electrons 2) 4)))))

(defclass spatial-arrangement () ())

(defclass tetrahedral-center
    (spatial-arrangement) ())

;;;; Hanser Ring Perception Algorithm

;;; For details see:
;;; Th. Hanser, Ph. Jauffret, and G. Kaufmann
;;; A New Algorithm for Exhaustive Ring Perception in a Molecular Graph 
;;; J. Chem. Inf. Comput. Sci. 1996, 36, 1146-1152 

(defun pairs (list)
  (loop for x on list
     append (loop with y = (car x) for z in (cdr x)
               collect (cons y z))))

(defun hanser-rings (graph)
  (let ((edge-hash (make-hash-table))
        rings)
    (labels ((append-paths (graph p1 p2 x)
               (when (> (length p2) (length p1)) (rotatef p1 p2))
               (if (eql x (car p1))
                   (append (graph:node-remove graph x p2)
                           (cons x (graph:node-remove graph x p1)))
                   (append (graph:node-remove graph x p1)
                           (cons x (graph:node-remove graph x p2)))))
             (hanser-convert (graph)
               (let ((graph (graph:copy-graph graph)))
                 (graph:map-edges (lambda (edge)
                                    (setf (gethash edge edge-hash)
                                          (graph:edge-nodes edge)))
                                  graph)
                 graph))
             (hanser-remove (graph x)
               (loop for (first . second)
                  in (pairs (remove-if (lambda (edge) (graph:self-edge-p graph edge))
                                       (graph:find-edges-containing graph x)))
                  do (let ((new-edge (graph:add-edge-between-nodes
                                      graph
                                      (graph:other-edge-node first x)
                                      (graph:other-edge-node second x))))
                       (setf (gethash new-edge edge-hash) 
                             (append-paths graph
                                           (gethash first edge-hash)
                                           (gethash second edge-hash)
                                           x))))
               (loop for path in (graph:find-edges-containing graph x)
                  do (when (graph:self-edge-p graph path) 
                       (push (gethash path edge-hash) rings))
                    (graph:remove-edge graph path))
               (graph:remove-node graph x)))
      (do ((v (hanser-convert graph)))
          ((zerop (graph:node-count v)))
        (hanser-remove v (car (sort (graph:map-nodes->list #'identity v)
                                    #'<
                                    :key (lambda (node)
                                           (length (graph:neighbors v node))))))))
    rings))
