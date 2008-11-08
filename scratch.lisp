
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

(defparameter *cyclohexane*
  (make-molecule :name "cyclohexane"))
(defparameter *c1* (make-atom :atomic-number 6 :name "C1"))
(defparameter *c2* (make-atom :atomic-number 6 :name "C2"))
(defparameter *c3* (make-atom :atomic-number 6 :name "C3"))
(defparameter *c4* (make-atom :atomic-number 6 :name "C4"))
(defparameter *c5* (make-atom :atomic-number 6 :name "C5"))
(defparameter *c6* (make-atom :atomic-number 6 :name "C6"))

(add-edge *cyclohexane* *c1* *c2*)
(add-edge *cyclohexane* *c2* *c3*)
(add-edge *cyclohexane* *c3* *c4*)
(add-edge *cyclohexane* *c4* *c5*)
(add-edge *cyclohexane* *c5* *c6*)
(add-edge *cyclohexane* *c6* *c1*)
 
;;; new stuff

(defparameter *graph* (make-instance 'edge-list-graph))
(defparameter *n1* (make-instance 'node :name "Node 1"))
(defparameter *n2* (make-instance 'node :name "Node 2"))
(defparameter *n3* (make-instance 'node :name "Node 3"))
(defparameter *n4* (make-instance 'node :name "Node 4"))
(defparameter *n5* (make-instance 'node :name "Node 5"))
(defparameter *n6* (make-instance 'node :name "Node 6"))
(defparameter *n7* (make-instance 'node :name "Node 7"))
(defparameter *n8* (make-instance 'node :name "Node 8"))
(defparameter *n9* (make-instance 'node :name "Node 9"))
(defparameter *n10* (make-instance 'node :name "Node 10"))
(defparameter *n11* (make-instance 'node :name "Node 11"))

(add-edge *graph* *n1* *n2*)
(add-edge *graph* *n1* *n3*)
(add-edge *graph* *n3* *n4*)
(add-edge *graph* *n3* *n5*)
(add-edge *graph* *n2* *n6*)
(add-edge *graph* *n6* *n7*)
(add-edge *graph* *n6* *n8*)
(add-edge *graph* *n4* *n9*)
(add-edge *graph* *n5* *n10*)
(add-edge *graph* *n5* *n11*)

(edgep *graph* *n1* *n2*)

(remove-edge *graph* *n1* *n2*)

;;;

(cl-graph:element (make-element-vertex "C"))

(add-vertex *cyclohexane* (make-element-vertex "C") :if-duplicate-do :force)
(add-vertex *cyclohexane* :bogus)

(defun list-vertices (graph)
  (let (l) 
    (iterate-vertexes graph (lambda (x) (push x l)))
    (nreverse l)))

(defun list-edges (graph)
  (let (l)
    (iterate-edges graph (lambda (x) (push x l)))
    (nreverse l)))

(defun list->graph (l)
  (let ((g (make-graph 'graph-container)))
    (loop for (source . target) in l
       do
	 (add-edge-between-vertexes g source target))
    g))

(defparameter *cyclohexane*
  (make-molecule :name "cyclohexane"))
(defparameter *c1* (make-atom :atomic-number 6 :name "C1"))
(add-edge *cyclohexane* *c1* *c2*)

(defparameter *methane*
  (let ((mol (make-molecule :name "methane")))
    (let ((c1 (make-atom :c :name "C1")))
      (add-edge mol c1 (make-atom :h :name "H1"))
      (add-edge mol c1 (make-atom :h :name "H2"))
      (add-edge mol c1 (make-atom :h :name "H3"))
      (add-edge mol c1 (make-atom :h :name "H4")))
    mol))

(defparameter *ethane*
  (let ((mol (make-molecule :name "ethane")))
    (let ((c1 (make-atom :c :name "C1"))
          (c2 (make-atom :c :name "C2")))
      (add-edge mol c1 c2)
      (add-edge mol c1 (make-atom :h :name "H1"))
      (add-edge mol c1 (make-atom :h :name "H2"))
      (add-edge mol c1 (make-atom :h :name "H3"))
      (add-edge mol c2 (make-atom :h :name "H4"))
      (add-edge mol c2 (make-atom :h :name "H5"))
      (add-edge mol c2 (make-atom :h :name "H6")))
    mol))

