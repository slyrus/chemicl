
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

(defparameter *benzene*
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
                   :order 1.5)
         (add-bond mol
                   (format nil "C~A" i)
                   (format nil "H~A" i)))
    mol))

(find-atom *benzene* "C1")
(find-edges-containing *benzene* (find-atom *benzene* "C1"))
(atom-bond-order *benzene* "C1")
(atom-bond-order *benzene* "H1")

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

(defparameter *methane*
  (let ((mol (make-molecule :name "methane")))
    (let ((c1 (make-atom :c :name "C1")))
      (add-bond mol c1 (make-atom :h :name "H1"))
      (add-bond mol c1 (make-atom :h :name "H2"))
      (add-bond mol c1 (make-atom :h :name "H3"))
      (add-bond mol c1 (make-atom :h :name "H4")))
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

