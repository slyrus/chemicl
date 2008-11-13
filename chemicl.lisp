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

(defmethod print-object ((object element) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (format stream
            "~S ~S ~S ~S"
            (atomic-number object)
            (id object)
            (name object)
            (mass object))))

(defclass atom (node)
  ((element :initarg :element :accessor element :initform nil)
   (molecule :initarg :molecule :accessor molecule :initform nil))
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
    (format stream "~S ~S" 
            (if (slot-boundp object 'element)
                (when (element object)
                  (if (> *atom-print-verbosity* 0)
                      (element object)
                      (id (element object))))
                "Element: unbound")
            (node-name object))))

(defmethod mass ((atom atom))
  (mass (element atom)))

(defclass molecule (simple-edge-list-graph)
  ((name :initarg :name :accessor name))
  (:documentation "A class for representing molecules."))

(defmethod add-atom ((molecule molecule) identifier name)
  (let ((atom (make-atom identifier :name name :molecule molecule)))
    (add-node molecule atom)
    atom))

(defun find-atom (molecule atom-identifier)
  (typecase atom-identifier
    (atom atom-identifier)
    (string (get-node molecule atom-identifier))))

(defclass bond (edge)
  ((bond-type :accessor bond-type :initarg :bond-type :initform :single)))

(defmethod add-bond ((molecule molecule) atom-identifier-1 atom-identifier-2
                     &key (bond-type :single))
  (let ((atom-1 (find-atom molecule atom-identifier-1))
        (atom-2 (find-atom molecule atom-identifier-2)))
    (let ((bond (make-instance 'bond
                               :node1 atom-1
                               :node2 atom-2
                               :bond-type bond-type)))
      (add-edge molecule bond)
      bond)))

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
    (dfs-map molecule
             (first-node molecule)
             (lambda (molecule)
               (incf mass (mass molecule))))
    mass))

(defun smiles->molecule (string))
