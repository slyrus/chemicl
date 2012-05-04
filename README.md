
# Introduction

chemicl is a common library for representing chemical structures,
namely molecules, which are composed of atoms and bonds between the
atoms.

# Source

The source code for chemicl can be found on
[github](https://github.com/slyrus/chemicl).

# Examples

## Preliminaries

First, let's make a package for use with chemicl. One thing to watch
out for is that there is an exported symbol in the chemicl package
which conflicts with a symbol in the common-lisp package: atom. One
most take care to import the proper symbol. For my packages that use
chemicl, I usually do a shadowing import from chemicl for the atom
symbol:

    (asdf:load-system 'chemicl)

    (cl:defpackage #:chemicl-user
      (:use #:cl #:chemicl)
      (:shadowing-import-from #:chemicl #:atom))

    (cl:in-package #:chemicl-user)

## A simple example

    (defparameter *cyclohexane*
      (let ((mol (make-molecule :name "cyclohexane")))
        (add-atom mol 6 "C1")
        (add-atom mol 6 "C2")
        (add-atom mol 6 "C3")
        (add-atom mol 6 "C4")
        (add-atom mol 6 "C5")
        (add-atom mol 6 "C6")

        (loop for i from 1 to 12
           do (add-atom mol 1 (format nil "H~A" i)))
        (loop for i from 1 to 6
           do (add-bond mol
                        (format nil "C~A" (1+ (mod (1- i) 6)))
                        (format nil "C~A" (1+ (mod i 6)))))
        (loop for i from 1 to 6
           do (add-bond mol
                        (format nil "C~A" i)
                        (format nil "H~A" (1+ (* (1- i) 2))))
             (add-bond mol
                       (format nil "C~A" i)
                       (format nil "H~A" (+ (* (1- i) 2) 2))))
        mol))

## Parsing a SMILES string:

Creating a molecule using the add-atom/add-bond functions as above is
a bit cumbersone. Another way to do so is to specify the molecule
using a
[SMILES](http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
string:

    (defparameter *acetominophen*
      (chem:parse-smiles-string "CC(=O)NC1=CC=C(C=C1)O"
                                :name "acetominophen"))

Then we can do things like:
    
    CHEMICL-USER> (mass *acetominophen*)
    151.16255999999993d0

# Dictionary

## Element

### [class] element

    (defclass element ()
      (;; modeling elements after the entries in data/elementdata.xml

       ;; from XML attributes
       (atomic-number :initarg :atomic-number :accessor atomic-number)
       (id :initarg :id :accessor id)
       (name :initarg :name :accessor name)
       (group :initarg :group :accessor group)
       (period :initarg :period :accessor period)

       ;; from XML elements
       (radii :initarg :radii :accessor radii :initform nil)
       (max-bond-order :initarg :max-bond-order :accessor max-bond-order)
       (mass :initarg :mass :accessor mass)
       (electronegativity :initarg :electronegativity :accessor electronegativity)

       (isotopes :initarg :isotopes :accessor isotopes :initform nil))
      (:documentation "A class for representing elements of the periodic table."))

### [class] isotope

    (defclass isotope ()
      ((number :initarg :number :accessor isotope-number)
       (exact-mass :initarg :exact-mass :accessor isotope-exact-mass)
       (relative-abundance :initarg :relative-abundance
                           :accessor isotope-relative-abundance)))

### [function] get-element element-identifier => element-object

### [function] get-isotope element-identifier mass-number => isotope-object

### [generic function] get-normal-valence element => valence-list

## Atom

### [class] atom

    (defclass atom ()
      ((name :initarg :name :accessor atom-name :initform nil)
       (element :initarg :element :accessor element :initform nil)
       (charge :initarg :charge :accessor charge :initform 0)
       (isotope :initarg :isotope :accessor isotope :initform nil))
      (:documentation "A class for representing individual atoms. For
      example, a molecule of hydrogen class would contain two atom
      instances, each of whose element slots would contain the (same)
      element instance for the hydrogen. An atom can be associated with at
      most one molecule at a time, as sepcified by its molecule slot,
      which can be NIL, indicating that the atom is not associated with
      any molecule."))

### [method] mass atom => result

### [method] exact-mass atom => result

### [method] get-normal-valence atom => result

### [function] make-atom element-identifier => atom-object

## Molecule

The MOLECULE class is a subclass of the epigraph:graph class and is
used to represent a given species of molecule, such as water, benzene,
alanine, etc...

### [class] molecule

    (defclass molecule (graph:simple-edge-list-graph)
      ((name :initarg :name :accessor name)
       (atom-name-hash :accessor atom-name-hash
                       :initform (make-hash-table :test 'equal)))
      (:documentation "A class for representing molecules."))

### [method] graph:copy-graph molecule &key copy-edges => new-molecule

### [generic function] get-atom molecule atom-identifier => atom-object

### [generic function] add-atom molecule atom => atom

### [generic function] remove-atom molecule atom => atom

### [generic function] atom-count molecule => count

## Bond

The BOND class is a subclass of the graph:edge class and is used to
represent chemical bonds between atoms in a molecule.

### [class] bond

    (defclass bond (graph:edge)
      ((type :accessor bond-type :initarg :type :initform :single)
       (order :accessor bond-order :initarg :order :initform 1)
       (direction :accessor bond-direction :initarg :direction :initform nil)))
