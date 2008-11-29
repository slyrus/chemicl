((:smarkup-metadata
  (:copyright
   "Copyright 2008, Cyrus Harmon. All Rights Reserved. See COPYRIGHT
  file for details")
  (:title "chemicl is a common lisp library for representing and
manipulating chemical entities and concepts such as molecules, atoms,
elements, etc...  nodges and edges.")
  (:author "Cyrus L. Harmon"))
 (:html-metadata (:htmlcss "simple.css") )
 
 (:span
  (:title "chemicl")

  (:h1 "Introduction")

  (:p "chemicl is a common library for representing chemical
  structures, namely molecules, which are composed of atoms and bonds
  between the atoms.")

  (:h1 "Examples")

  (:pre
   (:code #q{(defparameter *cyclohexane*
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
    mol))}))

  (:h1 "Dictionary")

  (:h2 "Element")

  (:list
   (:item "[class] "
     (:code "element")
     (:pre
      (:code #q{(defclass element ()
  (;; modeling elements after the entries in data/elementdata.xml

   ;; from XML attributes
   (atomic-number :initarg :atomic-number :accessor atomic-number)
   (id :initarg :id :accessor id)
   (name :initarg :name :accessor name)
   (group :initarg :group :accessor group)
   (period :initarg :period :accessor period)

   ;; from XML elements
   (radii :initarg :radii :accessor :radii :initform nil)
   (max-bond-order :initarg :max-bond-order :accessor :max-bond-order)
   (mass :initarg :mass :accessor mass)
   (electronegativity :initarg :electronegativity :accessor electronegativity)

   (isotopes :initarg :isotopes :accessor isotopes :initform nil))
  (:documentation "A class for representing elements of the periodic table."))})))

   (:item "[class] "
     (:code "isotope")
     (:pre
      (:code #q{(defclass isotope ()
  ((number :initarg :number :accessor isotope-number)
   (exact-mass :initarg :exact-mass :accessor isotope-exact-mass)
   (relative-abundance :initarg :relative-abundance
                       :accessor isotope-relative-abundance)))})))
   (:item "[function] "
     (:code "get-element element-identifier => element-object")
     (:p #.(documentation 'get-element 'function)))

   (:item "[function] "
     (:code "get-isotope element-identifier mass-number => isotope-object")
     (:p #.(documentation 'get-isotope 'function)))

   (:item "[generic function] "
     (:code "get-normal-valence element => valence-list")
     (:p #.(documentation 'get-normal-valence 'function))))

  (:h2 "Atom")

  (:list
   (:item "[class] "
     (:code "element")
     (:pre
      (:code #q{(defclass atom ()
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
  any molecule."))})))
   
   (:item "[method] "
     (:code "mass atom => result"))

   (:item "[method] "
     (:code "exact-mass atom => result"))

   (:item "[method] "
     (:code "get-normal-valence atom => result"))

   (:item "[function] "
     (:code "make-atom element-identifier => atom-object")
     (:p #.(documentation 'make-atom 'function))))

  (:h2 "Molecule")

  (:p "The " (:code "MOLECULE") " class is a subclass of the
  epigraph:graph class and is used to represent a given species of
  molecule, such as water, benzene, alanine, etc...")

  (:list
   (:item "[class] "
     (:pre
      (:code #q{(defclass molecule (graph:simple-edge-list-graph)
  ((name :initarg :name :accessor name)
   (atom-name-hash :accessor atom-name-hash
                   :initform (make-hash-table :test 'equal)))
  (:documentation "A class for representing molecules."))})))
   
   (:item "[method] "
     (:code "graph:copy-graph molecule &key copy-edges => new-molecule"))
   
   (:item "[generic function] "
     (:code "get-atom molecule atom-identifier => atom-object")
     (:p #.(documentation 'get-atom 'function)))

   (:item "[generic function] "
     (:code "add-atom molecule atom => atom")
     (:p #.(documentation 'add-atom 'function)))

   (:item "[generic function] "
     (:code "remove-atom molecule atom => atom")
     (:p #.(documentation 'remove-atom 'function)))

   (:item "[generic function] "
     (:code "atom-count molecule => count")
     (:p #.(documentation 'atom-count 'function))))

  
  (:h2 "Bond")
  (:p "The " (:code "BOND") " class is a subclass of the graph:edge
  class and is used to represent chemical bonds between atoms in a
  molecule.")

  (:list
   (:item "[class] "
     (:pre
      (:code #q{(defclass bond (graph:edge)
  ((type :accessor bond-type :initarg :type :initform :single)
   (order :accessor bond-order :initarg :order :initform 1)
   (direction :accessor bond-direction :initarg :direction :initform nil)))}))))))
