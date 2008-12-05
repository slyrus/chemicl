

(asdf:oos 'asdf:load-op :chemicl)

(in-package :chemicl)

(defun join-xpath-result (result)
  (if (xpath:node-set-p result)
      (format nil "~{~a~^~&~}"
              (xpath:map-node-set->list #'xpath:string-value result))
      (xpath:string-value result)))

(join-xpath-result 
 (xpath:evaluate
  (concatenate 'string "/elements"
               "/element[attribute::name=\"hydrogen\"]")
  *element-nodes*))

(xpath:map-node-set->list 
 (lambda (node)
   (stp:with-attributes (id name)
       node
     (list (mapcar #'stp:local-name (stp:list-attributes node)) id name)))
 (xpath:evaluate "/elements/element"
  *element-nodes*))

(join-xpath-result 
 (xpath:evaluate
  (concatenate 'string "/elements"
               "/element[attribute::name=\"hydrogen\"]/maxbondorder")
  *element-nodes*))

(xpath:evaluate
  (concatenate 'string "/elements"
               "/element[attribute::name=\"hydrogen\"]")
  *element-nodes*)

(stp:find-child-if (stp:of-name "maxbondorder")
                   (xpath:first-node
                    (xpath:evaluate
                     (concatenate
                      'string "/elements"
                      "/element[attribute::name=\"hydrogen\"]")
                     *element-nodes*)))

(xpath:map-node-set->list
 (lambda (node)
   (xpath:evaluate "maxbondorder" node))
 (xpath:evaluate
  (concatenate
   'string "/elements"
   "/element[attribute::name=\"hydrogen\"]")
  *element-nodes*))

#+nil
(defun list-vertices (graph)
  (let (l) 
    (iterate-vertexes graph (lambda (x) (push x l)))
    (nreverse l)))

#+nil
(defun list-edges (graph)
  (let (l)
    (iterate-edges graph (lambda (x) (push x l)))
    (nreverse l)))

#+nil
(defun list->graph (l)
  (let ((g (make-graph 'graph-container)))
    (loop for (source . target) in l
       do
	 (add-edge-between-vertexes g source target))
    g))

(defparameter *benzene-manual*
  (let ((mol (make-molecule :name "benzene")))
    (loop for i from 1 to 6
       do
         (add-atom mol 6 (format nil "C~A" i))
         (add-atom mol 1 (format nil "H~A" i)))
    (loop for i from 1 to 6
       do
         (add-bond mol
                   (format nil "C~A" (1+ (mod (1- i) 6)))
                   (format nil "C~A" (1+ (mod i 6)))
                   :type :aromatic
                   :order (if (oddp i) 2 1))
         (add-bond mol
                   (format nil "C~A" i)
                   (format nil "H~A" i)))
    mol))

(find-atom *benzene* "C1")
(graph:find-edges-containing *benzene* (find-atom *benzene* "C1"))
(atom-bond-order *benzene* "C1")
(atom-bond-order *benzene* "H1")

(defparameter *cyclohexane-manual*
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

(defparameter *methane*
  (let ((mol (make-molecule :name "methane")))
    (let ((c1 (add-atom mol :c "C1")))
      (add-bond mol c1 (add-atom mol :h "H1"))
      (add-bond mol c1 (add-atom mol :h "H2"))
      (add-bond mol c1 (add-atom mol :h "H3"))
      (add-bond mol c1 (add-atom mol :h "H4")))
    mol))

(defparameter *ethane*
  (let ((mol (make-molecule :name "ethane")))
    (let ((c1 (add-atom mol :c "C1"))
          (c2 (add-atom mol :c "C2")))
      (add-bond mol "C1" "C2")
      (add-bond mol "C1" (add-atom mol :h "H1"))
      (add-bond mol "C1" (add-atom mol :h "H2"))
      (add-bond mol "C1" (add-atom mol :h "H3"))
      (add-bond mol c2 (add-atom mol :h "H4"))
      (add-bond mol c2 (add-atom mol :h "H5"))
      (add-bond mol c2 (add-atom mol :h "H6")))
    mol))

;;; SMILES examples
(defparameter *methane* (parse-smiles-string "C" :name "methane"))
(defparameter *ethane* (parse-smiles-string "CC" :name "ethane"))
(defparameter *propane* (parse-smiles-string "CCC" :name "propane"))
(defparameter *butane* (parse-smiles-string "CCCC" :name "butane"))
(defparameter *pentane* (parse-smiles-string "CCCCC" :name "pentane"))
(defparameter *hexane* (parse-smiles-string "CCCCCC" :name "hexane"))

(defparameter *heavy-hexane* (parse-smiles-string "[13CH3]CCCCC" :name "hexane"))
(defparameter *cycloehexane* (parse-smiles-string "C1CCCCC1"
                                                  :name "cyclohexane"))
(defparameter *benzene* (parse-smiles-string "c1ccccc1"
                                             :name "benzene"))

(defparameter *pyrrole* (parse-smiles-string "[nH]1cccc1"
                                             :name "pyrrole"))

(defparameter *tamoxifen*
  (parse-smiles-string "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3"
                       :name "tamoxifen"))

(defparameter *z-tamoxifen*
  (parse-smiles-string "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3"
                       :name "z-tamoxifen"))

(defparameter *phentermine*
  (parse-smiles-string "CC(C)(N)Cc1ccccc1"))

(defparameter *Z-1.2-difluoroethene*
  (parse-smiles-string "F/C=C/F"))

(defparameter *E-1.2-difluoroethene*
  (parse-smiles-string "F/C=C\\F"))

(defparameter *furan*
  (parse-smiles-string "C1=COC=C1"))

(multiple-value-bind (mol bonds cycles)
    (parse-smiles-string "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3"
                         :add-implicit-hydrogens nil)
  (mapcar (lambda (x) (atom-bond-order mol x)) (butlast (car cycles))))

(multiple-value-bind (mol bonds cycles)
    (parse-smiles-string "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3"
                         :add-implicit-hydrogens nil)
  (mapcar (lambda (x) (bond-order x)) bonds))


(multiple-value-bind (edges nodes)
    (graph:find-cycle *z-tamoxifen*)
  (member-if (lambda (x) (and (find (graph:node1 x) nodes)
                              (find (graph:node2 x) nodes)))
             edges))

(member "C5" (graph:find-cycle *z-tamoxifen*)
        :key (lambda (x) (atom-name (graph:node1 x)))
        :test 'equal)


(defparameter *isotope-nodes*
  (with-cml-namespace
    (xpath:evaluate "/cml"
                    (cxml:parse-file 
                     (asdf:component-pathname
                      (let ((path '("chemicl" "data" "isotopes.xml")))
                        (reduce #'asdf:find-component (cdr path)
                                :initial-value (asdf:find-system (car path)))))
                     (stp:make-builder)))))

(with-cml-namespace
  (xpath:map-node-set->list
   (lambda (node)
     (xpath:map-node-set->list 
      (lambda (node)
        (stp:with-attributes ((element-type "elementType"))
            node
          (list element-type 
                (xpath:number-value
                 (xpath:evaluate "scalar[attribute::dictRef=\"bo:exactMass\"]" node))
                (xpath:number-value
                 (xpath:evaluate "scalar[attribute::dictRef=\"bo:relativeAbundance\"]" node)))))
      (xpath:evaluate "isotopeList/isotope[attribute::elementType=\"H\"]" node)))
   *isotope-nodes*))

(mapcar (lambda (x) (* (isotope-relative-abundance x)
                       (isotope-exact-mass x)))
        (isotopes (get-element 6)))



;;;
;;; CML stuff

(defparameter *concat* (cxml:parse-file "cml/concatenated.xml" (stp:make-builder)))

(xpath:with-namespaces ((nil "http://www.xml-cml.org/schema/cml2/core")
                           ("bo" "http://www.blueobelisk.org/dict/terminology" ))
        (xpath:evaluate "/list/molecule" *concat*))

(defparameter *benzene-core*
  (parse-smiles-string "C1=CC=CC=C1"
                       :name "benzene-core"
                       :add-implicit-hydrogens nil))

(defparameter *furan-core*
  (parse-smiles-string "O1C=CC=C1"
                       :name "furan-core"
                       :add-implicit-hydrogens nil))

(defparameter *cyclotetraoctene-core*
  (parse-smiles-string "C1=CC=CC=CC=C1"
                       :name "cyclotetraoctene"
                       :add-implicit-hydrogens nil))
