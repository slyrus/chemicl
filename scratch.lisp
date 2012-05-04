

(asdf:load-system 'chemicl)

(cl:defpackage #:chemicl-user
  (:use #:cl #:chemicl)
  (:shadowing-import-from #:chemicl #:atom))

(cl:in-package #:chemicl-user)

(defun join-xpath-result (result)
  (if (xpath:node-set-p result)
      (format nil "~{~a~^~&~}"
              (xpath:map-node-set->list #'xpath:string-value result))
      (xpath:string-value result)))

(defparameter *element-nodes*
  (cxml:parse-file *element-data-xml-pathname*
                   (stp:make-builder)))

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

(defparameter *alanine* (parse-smiles-string "CC(C(=O)O)N"
                                             :name "alanine"))

(defparameter *l-alanine* (parse-smiles-string "C[C@@H](C(=O)O)N"
                                               :name "l-alanine"))

(defparameter |*e-1,2-difluoroethene*| (parse-smiles-string "F/C=C/F"
                                                            :name "e-1,2-difluoroethene"))

(defparameter |*z-1,2-difluoroethene*| (parse-smiles-string "F\\C=C/F"
                                                            :name "z-1,2-difluoroethene"))

#+nil
(defparameter *tamoxifen*
  (parse-smiles-string "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3"
                       :name "tamoxifen"))

(defparameter *tamoxifen*
  (parse-smiles-string "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3"
                       :name "tamoxifen"))

(defparameter *z-tamoxifen*
  (parse-smiles-string "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3"
                       :name "z-tamoxifen"))

(defparameter *quetiapine*
  (parse-smiles-string "C1CN(CCN1CCOCCO)C2=NC3=CC=CC=C3SC4=CC=CC=C42"))

(defparameter *morphine*
  (parse-smiles-string "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"))

#+nil
(defparameter *morphine-isomeric*
  (parse-smiles-string "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"))

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
                     *isotope-data-xml-pathname*
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

(defparameter *smilanol*
  (parse-smiles-string "OCC(CC)CCC(CN)CN"))

(defun molecule-info (molecule)
  (format t "~&Molecular Formula ~{~{~A~A ~}~}~%"
          (mapcar (lambda (x)
                    (list (id (car x))
                          (cdr x)))
                  (count-elements molecule)))
  (format t "~&Molecular Mass: ~A~%" (mass molecule))
  (format t "~&Exact Mass: ~A~%" (exact-mass molecule)))

(defun smiles-info (smiles)
  (let ((molecule (parse-smiles-string smiles)))
    (molecule-info molecule)))

(molecule-info
 (parse-smiles-string
  "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3"
  :name "z-tamoxifen"))


(let* ((molecule *smilanol*)
       (atoms (get-non-h-atoms molecule))
       (invariants (canon-invariants molecule atoms))
       (ranks (rank-order invariants #'list<)))

  (let ((primes (mapcar #'nth-prime ranks)))
    (let ((product-of-primes
           (compute-product-of-primes primes atoms molecule)))
      (list
       (mapcar
        #'1+
        (rank-order
         (mapcar #'list product-of-primes ranks)
         #'list<))
       atoms))))

(defparameter *tamoxifen-inchi*
  "InChI=1/C26H29NO/c1-4-25(21-11-7-5-8-12-21)26(22-13-9-6-10-14-22)23-15-17-24(18-16-23)28-20-19-27(2)3/h5-18H,4,19-20H2,1-3H3/b26-25-")




(defparameter *graph* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (loop for i from 1 below 12
       do (graph:add-node g (concatenate 'string "Node " (princ-to-string i))))

    (graph:add-edge-between-nodes g "Node 5" "Node 1")
    (graph:add-edge-between-nodes g "Node 1" "Node 2")
    (graph:add-edge-between-nodes g "Node 2" "Node 3")
    (graph:add-edge-between-nodes g "Node 3" "Node 4")
    (graph:add-edge-between-nodes g "Node 4" "Node 2")
    g))

(defparameter *graph-2* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(:a :b :c :d :e :f :g :h :i :j :k :l :m :n :o :p))
    (graph:add-edge-between-nodes g :a :b)
    (graph:add-edge-between-nodes g :b :c)
    (graph:add-edge-between-nodes g :b :d)
    (graph:add-edge-between-nodes g :d :e)
    (graph:add-edge-between-nodes g :e :f)
    (graph:add-edge-between-nodes g :f :g)
    (graph:add-edge-between-nodes g :g :e)
    (graph:add-edge-between-nodes g :d :h)
    (graph:add-edge-between-nodes g :h :i)
    (graph:add-edge-between-nodes g :i :j)
    (graph:add-edge-between-nodes g :j :m)
    (graph:add-edge-between-nodes g :m :n)
    (graph:add-edge-between-nodes g :n :o)
    (graph:add-edge-between-nodes g :o :p)
    (graph:add-edge-between-nodes g :p :l)
    (graph:add-edge-between-nodes g :l :m)
    (graph:add-edge-between-nodes g :l :k)
    (graph:add-edge-between-nodes g :k :h)
    g))


(defparameter *graph-3* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(:a :b :c :d :e :f :g :h :i :j :k :l :m :n :o :p :q))
    (graph:add-edge-between-nodes g :a :b)
    (graph:add-edge-between-nodes g :b :c)
    (graph:add-edge-between-nodes g :c :d)
    (graph:add-edge-between-nodes g :d :e)
    (graph:add-edge-between-nodes g :e :f)
    (graph:add-edge-between-nodes g :f :a)
    (graph:add-edge-between-nodes g :c :g)
    (graph:add-edge-between-nodes g :g :h)
    (graph:add-edge-between-nodes g :h :i)
    (graph:add-edge-between-nodes g :i :j)
    (graph:add-edge-between-nodes g :j :d)
    (graph:add-edge-between-nodes g :g :k)
    (graph:add-edge-between-nodes g :k :l)
    (graph:add-edge-between-nodes g :l :m)
    (graph:add-edge-between-nodes g :m :n)
    (graph:add-edge-between-nodes g :n :h)
    (graph:add-edge-between-nodes g :m :o)
    (graph:add-edge-between-nodes g :o :p)
    (graph:add-edge-between-nodes g :p :q)
    (graph:add-edge-between-nodes g :q :n)
    g))

(defparameter *graph-4* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(:a :b :c :d :e :f :g :h))
    (graph:add-edge-between-nodes g :a :b)
    (graph:add-edge-between-nodes g :a :d)
    (graph:add-edge-between-nodes g :a :h)
    (graph:add-edge-between-nodes g :b :c)
    (graph:add-edge-between-nodes g :b :g)
    (graph:add-edge-between-nodes g :c :d)
    (graph:add-edge-between-nodes g :c :f)
    (graph:add-edge-between-nodes g :d :e)
    (graph:add-edge-between-nodes g :e :f)
    (graph:add-edge-between-nodes g :e :h)
    (graph:add-edge-between-nodes g :f :g)
    (graph:add-edge-between-nodes g :g :h)
    g))

(defparameter *graph-5* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(1 2 3 4 5 6 7 8))
    (graph:add-edge-between-nodes g 1 2)
    (graph:add-edge-between-nodes g 2 3)
    (graph:add-edge-between-nodes g 3 4)
    (graph:add-edge-between-nodes g 1 4)
    (graph:add-edge-between-nodes g 4 5)
    (graph:add-edge-between-nodes g 5 6)
    (graph:add-edge-between-nodes g 3 6)
    (graph:add-edge-between-nodes g 6 7)
    (graph:add-edge-between-nodes g 2 7)
    (graph:add-edge-between-nodes g 7 8)
    (graph:add-edge-between-nodes g 1 8)
    (graph:add-edge-between-nodes g 5 8)
    g))

(defparameter *graph-6* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(1 2 3 4))
    (graph:add-edge-between-nodes g 1 2)
    (graph:add-edge-between-nodes g 2 3)
    (graph:add-edge-between-nodes g 3 4)
    (graph:add-edge-between-nodes g 1 4)
    g))

(defparameter *graph-7* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(1 2 3))
    (graph:add-edge-between-nodes g 1 2)
    (graph:add-edge-between-nodes g 2 3)
    (graph:add-edge-between-nodes g 1 3)
    
    g))

(defparameter *graph-8* 
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(1 2 3 4))
    (graph:add-edge-between-nodes g 1 2)
    (graph:add-edge-between-nodes g 2 3)
    (graph:add-edge-between-nodes g 1 3)
    (graph:add-edge-between-nodes g 3 4)
    
    g))

(defparameter *graph-8b*
  (let ((g (graph:make-graph :node-test 'equal)))
    (mapcar (lambda (x) (graph:add-node g x))
            '(1 2 3 4))
    (graph:add-edge-between-nodes g 1 2)
    (graph:add-edge-between-nodes g 2 3)
    (graph:add-edge-between-nodes g 3 1)
    (graph:add-edge-between-nodes g 3 4)
    
    g))

(mapcar #'find-duplicate
        (mapcar #'butlast (hanser-rings *graph-3*)))



(let ((mol *z-tamoxifen*))
  (mapcar (lambda (x) (cons (find-bonds-containing (atom1 x) mol) (atom2 x))) 
          (remove-if-not (lambda (x) (= (bond-order x) 2)) (bonds mol))))
        

(defparameter *chain-test*
  (parse-smiles-string "CCCCCCC(=O)C(CCC(O)CC)CCC(N)CCC"))

(defparameter *chain-test-2* (copy-molecule *chain-test*))

(remove-atoms-of-element *chain-test-2* "H")
(graph:find-longest-path *chain-test-2*)
(map-bonds #'print *chain-test-2*)
(graph:find-connected-components *chain-test*)

(defparameter *chain-test-3* (copy-molecule *chain-test-2*))
(remove-bond *chain-test-3* "C7" "C8")
(graph:find-connected-components *chain-test-3*)

;;; find the distance between two nodes
(let ((mol *chain-test*))
  (let ((atom1 (get-atom mol "C1"))
        (atom2 (get-atom mol "C12")))
    (graph:graph-distance mol atom1 atom2)))

(defparameter q (time (graph::graph-distance-matrix *chain-test*)))

(defparameter *acetominophen*
  (chem:parse-smiles-string "CC(=O)NC1=CC=C(C=C1)O"
                            :name "acetominophen"))
(defparameter *ring-test*
  (chem:parse-smiles-string "NN(N1OO1)(N2OON3OOON3O2)"
                            :name "ring-test"))


;; list the atoms:
(atoms *acetominophen*)

;; list the bonds:
(bonds *acetominophen*)

(let ((mol *ring-test*))
  ;; remove the hydrogens
  (let ((mol2 (remove-atoms-of-element (copy-molecule mol) "H")))
    
    ;; find the terminal atoms
    (let ((atoms-to-place (atoms mol2))
          terminal-atoms)
      (map-atoms (lambda (atom)
                   (when (= (length (graph:neighbors mol2 atom)) 1)
                     (push atom terminal-atoms)))
                 mol2)
      (let ((mol3 (copy-molecule mol2)))
        (map nil (lambda (atom)
                   (remove-atom mol3 atom))
             terminal-atoms)
        (let ((ring-graph (collapse-rings mol3)))
          (list (graph:nodes ring-graph)
                (graph:edges ring-graph)))))))
