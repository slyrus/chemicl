;;; file: package.lisp
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

(in-package #:cl-user)

(defpackage #:chemicl
  (:use #:cl)
  (:shadow #:atom)
  (:nicknames #:chem)
  (:export #:element
           #:id
           #:atomic-number
           #:name
           #:group
           #:period
           #:radii
           #:max-bond-order
           #:mass
           #:electronegativity
           #:charge
           
           #:get-element

           #:isotope
           #:isotope-number
           #:isotope-exact-mass
           #:isotope-relative-abundance

           #:atom
           #:make-atom
           #:atom-name

           #:molecule
           #:make-molecule
           #:copy-molecule
           #:add-atom
           #:find-atom
           #:remove-atom
           #:atoms
           #:map-atoms
           #:map-atoms->list
           #:atom-count

           #:bond
           #:bonds
           #:bond-order
           #:add-bond
           #:remove-bond
           #:find-bonds-containing

           #:exact-mass
           #:get-normal-valence
           #:count-element
           #:count-elements
           #:molecular-formula

           #:parse-smiles-string
           #:write-smiles-string
           
           #:*element-data-xml-pathname*
           #:*isotope-data-xml-pathname*))

(in-package #:chemicl)

(defparameter *element-data-xml-pathname*
  (merge-pathnames #p"data/elementdata.xml" chemicl-config:*base-directory*))

(defparameter *isotope-data-xml-pathname*
  (merge-pathnames #p"data/isotopes.xml" chemicl-config:*base-directory*))

