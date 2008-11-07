
(in-package :chemicl)

(defclass node ()
  ((name :accessor node-name :initarg :name)
   (data :accessor node-data :initarg :data)))

(defmethod print-object ((object node) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (format stream "~S" (node-name object))))

(defclass graph ()
  ((nodes :accessor graph-nodes :initarg :nodes :initform nil)))

(defclass edge-list-graph (graph)
  ((edges :accessor graph-edges :initarg :edges :initform nil)))

(defgeneric add-edge (graph node1 node2))

(defmethod add-edge ((graph edge-list-graph) (node1 node) (node2 node))
  (pushnew node1 (graph-nodes graph))
  (pushnew node2 (graph-nodes graph))
  (let ((edge (cons node1 node2)))
    (pushnew edge (graph-edges graph) :test 'equalp)))

(defmethod remove-edge ((graph edge-list-graph) (node1 node) (node2 node))
  (let ((edge (cons node1 node2)))
    (setf (graph-edges graph)
          (remove edge (graph-edges graph) :test 'equalp))))

(defmethod edgep ((graph edge-list-graph) (node1 node) (node2 node))
  (find (cons node1 node2) (graph-edges graph) :test 'equalp))

(defmethod find-edges-from ((graph edge-list-graph) node)
  (remove-if-not (lambda (x)
                   (eq node x))
                 (graph-edges graph) :key #'car))

(defmethod find-edges-to ((graph edge-list-graph) node)
  (remove-if-not (lambda (x)
                   (eq node x))
                 (graph-edges graph) :key #'cdr))

(defmethod find-edges-containing ((graph edge-list-graph) node)
  (union (find-edges-from graph node)
         (find-edges-to graph node)))

(defmethod bfs ((graph edge-list-graph) start-node query-node)
  (let ((visited (list start-node)))
    (labels
        ((bfs-visit (node-set-list)
           (let (children)
             (map nil
                  (lambda (node-and-path)
                    (destructuring-bind (node . path)
                        node-and-path
                      (when (eq node query-node)
                        (return-from bfs
                          (nreverse (cons node path))))
                      (let ((edges (find-edges-from graph node)))
                        (let ((neighbors (map (type-of edges) #'cdr edges)))
                          (map nil
                               (lambda (x)
                                 (unless (member x visited)
                                   (push (cons x (cons node path)) children)))
                               neighbors)))))
                  node-set-list)
             (when children
               (bfs-visit children)))))
      (bfs-visit (list (cons start-node nil))))))

(defmethod dfs ((graph edge-list-graph) start-node query-node)
  (let ((visited (list start-node)))
    (labels ((dfs-visit (node path)
               (if (eq node query-node)
                   (return-from dfs (nreverse (cons node path)))
                   (let ((edges (find-edges-from graph node)))
                     (let ((neighbors (map (type-of edges) #'cdr edges)))
                       (map nil
                            (lambda (x)
                              (unless (member x visited)
                                (push x visited)
                                (dfs-visit x (cons node path))))
                            neighbors))))))
      (dfs-visit start-node nil))))

