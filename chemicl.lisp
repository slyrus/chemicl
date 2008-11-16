;;; file: chemicl.lisp
;;;
;;; Copyright (c) 2008 Cyrus Harmon (ch-lisp@bobobeach.com)
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
;;;
;;;   * Redistributions of source code must retain the above copyright
;;;     notice, this list of conditions and the following disclaimer.
;;;
;;;   * Redistributions in binary form must reproduce the above
;;;     copyright notice, this list of conditions and the following
;;;     disclaimer in the documentation and/or other materials
;;;     provided with the distribution.
;;;
;;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS' AND ANY EXPRESSED
;;; OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
;;; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
;;; GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
;;; WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;

(in-package :chemicl)

(defclass element ()
  (;; from XML attributes
   (atomic-number :initarg :atomic-number :accessor atomic-number)
   (id :initarg :id :accessor id)
   (name :initarg :name :accessor name)
   (group :initarg :group :accessor group)
   (period :initarg :period :accessor period)

   ;; from XML elements
   (radii :initarg :radii :accessor :radii :initform nil)
   (max-bond-order :initarg :max-bond-order :accessor :max-bond-order)
   (mass :initarg :mass :accessor mass)
   (electronegativity :initarg :electronegativity :accessor electronegativity))
  (:documentation "A class for representing elements of the periodic table."))

(defgeneric print-element-data (object stream)
  (:method ((object element) stream)
    (format stream
            "~S ~S ~S ~S"
            (atomic-number object)
            (id object)
            (name object)
            (mass object))))

(defmethod print-object ((object element) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (print-element-data object stream)))

(defclass atom ()
  ((name :initarg :name :accessor atom-name :initform nil)
   (element :initarg :element :accessor element :initform nil)
   (charge :initarg :charge :accessor charge :initform 0)
   (isotope-mass :initarg :isotope-mass :accessor isotope-mass :initform nil))
  (:documentation "A class for representing individual atoms. For
  example, a molecule of hydrogen class would contain two atom
  instances, each of whose element slots would contain the (same)
  element instance for the hydrogen. An atom can be associated with at
  most one molecule at a time, as sepcified by its molecule slot,
  which can be NIL, indicating that the atom is not associated with
  any molecule."))

;;; FIXME add a shared-initialize to initialize the exact atomic mass of the
;;; atom etc... !!!

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

(defmethod mass ((atom atom))
  (mass (element atom)))

(defmethod exact-mass ((atom atom))
  (mass (element atom)))


;;; We need a way of explicitly storing configurations around double
;;; bonds and configurations at chiral centers. Hrm...
(defclass molecule (simple-edge-list-graph)
  ((name :initarg :name :accessor name)
   (atom-name-hash :accessor atom-name-hash
                   :initform (make-hash-table :test 'equal)))
  (:documentation "A class for representing molecules."))

(defgeneric get-atom (molecule name)
  (:method ((molecule molecule) (atom atom))
    atom)
  (:method ((molecule molecule) name)
    (gethash name (atom-name-hash molecule))))

(defgeneric add-atom (molecule identifier name)
  (:method ((molecule molecule) identifier name)
    (let ((atom (make-atom identifier :name name)))
      (add-node molecule atom)
      (setf (gethash name (atom-name-hash molecule)) atom)
      atom)))

(defgeneric atom-count (molecule)
  (:method ((molecule molecule))
    (node-count molecule)))

(defun find-atom (molecule atom-identifier)
  (typecase atom-identifier
    (atom atom-identifier)
    (string (get-atom molecule atom-identifier))))

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

(defclass bond (edge)
  ((type :accessor bond-type :initarg :type :initform :single)
   (order :accessor bond-order :initarg :order :initform 1)
   (direction :accessor bond-direction :initarg :direction :initform nil)))

(defmethod print-edge-data :after ((object bond) stream)
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
  (let ((bonds (find-edges-containing molecule atom)))
    (reduce #'+ bonds :key 'bond-order)))

(defmethod atom-bond-order ((molecule molecule) atom-identifier)
  (let ((atom (find-atom molecule atom-identifier)))
    (when atom
      (atom-bond-order molecule atom))))

(defgeneric add-bond (molecule atom-identifier-1 atom-identifier-2
                               &key type order direction))

(defmethod add-bond ((molecule molecule) atom-identifier-1 atom-identifier-2
                     &key (type :single) (order 1) direction)
  (let ((atom-1 (find-atom molecule atom-identifier-1))
        (atom-2 (find-atom molecule atom-identifier-2)))
    (let ((bond (apply #'make-instance 'bond
                       :node1 atom-1
                       :node2 atom-2
                       (append
                        (when type `(:type ,type))
                        (when order `(:order ,order))
                        (when direction `(:direction ,direction))))))
      (add-edge molecule bond)
      bond)))

(defun double-bond-configuration (molecule bond)
  (let ((atom1 (node1 bond))
        (atom2 (node2 bond)))
    (let ((atom1-bonds
           (remove bond (find-edges-containing molecule atom1)))
          (atom2-bonds
           (remove bond (find-edges-containing molecule atom2))))
      (cons atom1-bonds atom2-bonds))))

;;;
;;; B (3), C (4), N (3,5), O (2), P (3,5), S (2,4,6), and 1 for the
;;; halogens

;;; read in the element data from elementdata.xml, parse it and store
;;; in the *elments* array, with the index into the array specified
;;; by the atomic number of the element.
(defmacro parse-integer-if (string-or-nil)
    `(when ,string-or-nil
       (parse-integer ,string-or-nil)))

(defvar *elements*)

(defparameter *element-nodes*
    (cxml:parse-file 
     (asdf:component-pathname
      (let ((path '("chemicl" "elementdata.xml")))
        (reduce #'asdf:find-component (cdr path)
                :initial-value (asdf:find-system (car path)))))
     (stp:make-builder)))

(defvar *element-hash* (make-hash-table :test 'equalp))

(macrolet ((xpath-number (local-name parent-node)
             `(xpath:number-value
               (xpath:evaluate
                ,local-name
                ,parent-node))))
  (defun read-element-data ()
    (let ((element-list
           (xpath:map-node-set->list 
            (lambda (node)
              (stp:with-attributes ((atomic-number "atomicnumber")
                                    id
                                    name
                                    group
                                    period)
                  node
                (let ((max-bond-order (xpath-number "maxbondorder" node))
                      (mass (xpath-number "mass" node))
                      (electronegativity (xpath-number "electronegativity" node)))
                  (make-instance 'element
                                 :atomic-number (parse-integer-if atomic-number)
                                 :id id
                                 :name name
                                 :group group
                                 :period (parse-integer-if period)
                                 :mass mass
                                 :electronegativity electronegativity
                                 :max-bond-order max-bond-order))))
            (xpath:evaluate "/elements/element"
                            *element-nodes*))))
      (let ((max-element (apply #'max (map 'list #'atomic-number element-list))))
        (let ((array (make-array (1+ max-element) :adjustable nil)))
          (loop for l in element-list
             do (setf (aref array (atomic-number l)) l
                      (gethash (id l) *element-hash*) l))
          (setf *elements* array))))))

(read-element-data)

(defun get-element (identifier)
  "Gets the element indicated by identifier. If identifer is a number,
gets the element whose atomic number is identifier. If identifier is a
string, gets the element whose symbol is identifier."
  (etypecase identifier
    (number (aref *elements* identifier))
    (string (gethash identifier *element-hash*))
    (symbol (gethash (symbol-name identifier) *element-hash*))))

(defparameter *element-normal-valences*
  (let ((hash (make-hash-table :test 'eq))
        (valence-list '(("B" 3)
                        ("C" 4)
                        ("N" 3 5)
                        ("O" 2)
                        ("P" 3 5)
                        ("S" 2 4 6)
                        ("F" 1)
                        ("Cl" 1)
                        ("Br" 1)
                        ("I" 1))))
    (map nil (lambda (x)
               (destructuring-bind (symbol &rest valences)
                   x
                 (let ((element (get-element symbol)))
                   (setf (gethash element hash) valences))))
         valence-list)
    hash))

(defgeneric get-normal-valence (element)
  (:method ((element element))
    (gethash element *element-normal-valences*))
  (:method ((atom atom))
    (get-normal-valence (element atom)))
  (:method ((string string))
    (get-normal-valence (get-element string))))

(defun make-molecule (&rest args)
  (apply #'make-instance 'molecule args))

(defun remove-keyword-args (keywords list)
  (if (listp keywords)
      (loop for (x y) on list by #'cddr
         append (unless (member x keywords)
                  (list x y)))
      (loop for (x y) on list by #'cddr
         append (unless (eq x keywords)
                  (list x y)))))

(defun make-atom (identifier &rest args)
  (apply #'make-instance
         'atom
         :element (get-element identifier)
         args))

(defmethod mass ((molecule molecule))
  (let ((mass 0.0d0))
    (map-nodes molecule
             (lambda (atom)
               (incf mass (mass atom))))
    mass))

(defmethod charge ((molecule molecule))
  (let ((charge 0))
    (map-nodes molecule
             (lambda (atom)
               (incf charge (charge atom))))
    charge))

(defun count-element (molecule element-id)
  (let ((element-count 0)
        (element (get-element element-id)))
    (dfs-map molecule (first-node molecule)
             (lambda (atom)
               (when (eq (element atom) element)
                 (incf element-count))))
    element-count))

(defun count-elements (molecule)
   (let ((element-count-hash (make-hash-table)))
    (dfs-map molecule (first-node molecule)
             (lambda (atom)
               (setf (gethash (element atom) element-count-hash)
                     (1+ (or (gethash (element atom)
                                      element-count-hash)
                             0)))))
    (let (l)
      (maphash (lambda (k v)
                 (push (cons k v) l))
               element-count-hash)
      (sort l #'< :key (lambda (x) (atomic-number (car x)))))))

