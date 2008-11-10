
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

(defparameter *cyclohexane*
  (make-molecule :name "cyclohexane"))
(defparameter *c1* (make-atom 6 :name "C1"))
(defparameter *c2* (make-atom 6 :name "C2"))
(defparameter *c3* (make-atom 6 :name "C3"))
(defparameter *c4* (make-atom 6 :name "C4"))
(defparameter *c5* (make-atom 6 :name "C5"))
(defparameter *c6* (make-atom 6 :name "C6"))

(add-edge *cyclohexane* *c1* *c2*)
(add-edge *cyclohexane* *c2* *c3*)
(add-edge *cyclohexane* *c3* *c4*)
(add-edge *cyclohexane* *c4* *c5*)
(add-edge *cyclohexane* *c5* *c6*)
(add-edge *cyclohexane* *c6* *c1*)
 
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



