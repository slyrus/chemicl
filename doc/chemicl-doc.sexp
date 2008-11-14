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
   (:code #q{
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
}))

  (:h1 "Dictionary")

  (:h2 "Element")

  (:h2 "Atom")

  (:h2 "Molecule")

  (:h2 "Bond")))
