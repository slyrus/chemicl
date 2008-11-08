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
  ((element :initarg :element :accessor element :initform nil))
  (:documentation "A class for representing individual atoms. For
  example, a molecule of hydrogen class would contain two atom
  instances, each of whose element slots would contain the (same)
  element instance for the hydrogen."))

(defparameter *atom-print-verbosity* 0)

(defmethod print-object ((object atom) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (format stream "~S ~S" 
            (if (> *atom-print-verbosity* 0)
                (element object)
                (id (element object)))
            (node-name object))))

(defmethod mass ((atom atom))
  (mass (element atom)))

(defclass molecule (edge-list-graph)
  ((name :initarg :name :accessor name)
   (node-class :initform 'atom))
  (:documentation "A class for representing molecules."))


;;; read in the element data from elementdata.xml, parse it and store
;;; in the *elments* array, with the index into the array specified
;;; by the atomic number of the element.
(defmacro parse-integer-if (string-or-nil)
    `(when ,string-or-nil
       (parse-integer ,string-or-nil)))

(defvar *elements*)

(defparameter *element-nodes*
    (cxml:parse-file "elementdata.xml" (stp:make-builder)))

(macrolet ((xpath-number (local-name parent-node)
             `(xpath:number-value
               (xpath:evaluate
                ,local-name
                ,parent-node))))
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
           do (setf (aref array (atomic-number l))
                    l))
        (setf *elements* array)))))

(defun get-element (atomic-number)
  (aref *elements* atomic-number))

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

(defun make-atom (&rest args &key atomic-number &allow-other-keys)
  (apply #'make-instance
         'atom
         :element (get-element atomic-number)
         (remove-keyword-args :atomic-number args)))

(defun make-element-vertex (id)
  (let ((element (find id *elements* :key #'id :test #'string-equal)))
    (make-instance 'element-vertex :vertex-element element)))

(defmethod mass ((molecule molecule))
  (let ((mass 0.0d0))
    (dfs-map molecule
             (car (graph-nodes molecule))
             (lambda (molecule)
               (incf mass (mass molecule))))
    mass))
