;;; file: rings.lisp
;;;
;;; Copyright (c) 2008-2009 Cyrus Harmon (ch-lisp@bobobeach.com)
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

(in-package :chemicl)

(defclass ring-node (atom-container)
  ((nodes :accessor graph:nodes :initarg :nodes)
   (edges :accessor graph:edges :initarg :edges)))

;;;; Hanser Ring Perception Algorithm

;;; For details see:
;;; Th. Hanser, Ph. Jauffret, and G. Kaufmann
;;; A New Algorithm for Exhaustive Ring Perception in a Molecular Graph
;;; J. Chem. Inf. Comput. Sci. 1996, 36, 1146-1152

(defun pairs (list)
  (loop for x on list
     append (loop with y = (car x) for z in (cdr x)
               collect (cons y z))))

(defun hanser-rings (graph)
  "Returns a list of condensed cycles contained in the graph."
  (let ((edge-hash (make-hash-table))
        (graph (graph:copy-graph graph))
        rings)
    (graph:map-edges->list (lambda (edge) (setf (gethash edge edge-hash)
                                                (graph:nodes edge)))
                           graph)
    (labels ((append-paths (p1 p2 x)
               (when (> (length p2) (length p1)) (rotatef p1 p2))
               (let ((l (cond ((and (equal x (car p1))
                                    (equal x (car p2)))
                               (append
                                (reverse (graph:node-remove graph x p2))
                                (list x)
                                (graph:node-remove graph x p1)))
                              ((eql x (car p1))
                               (append
                                (graph:node-remove graph x p2)
                                (cons x (graph:node-remove graph x p1))))
                              (t (append
                                  (graph:node-remove graph x p1)
                                  (cons x (graph:node-remove graph x p2)))))))
                 ;; if there is a duplicate (other than the last
                 ;; node), then this isn't a "real" path and we can
                 ;; ignore it.
                 (unless (find-duplicate (butlast l))
                   l)))
             (hanser-remove (x)
               (loop for (first . second)
                  in (pairs (remove-if
                             (lambda (edge) (graph:self-edge-p graph edge))
                             (graph:find-edges-containing graph x)))
                  do 
                    (let ((path (append-paths (gethash first edge-hash)
                                              (gethash second edge-hash)
                                              x)))
                      (when path
                        (let ((new-edge (graph:add-edge-between-nodes
                                         graph
                                         (graph:other-edge-node first x)
                                         (graph:other-edge-node second x))))
                          (setf (gethash new-edge edge-hash) path)))))
               (loop for path in (graph:find-edges-containing graph x)
                  do (when (graph:self-edge-p graph path)
                       (push (gethash path edge-hash) rings))
                    (graph:remove-edge graph path))
               (graph:remove-node graph x)))
      (map nil #'hanser-remove
           (sort (graph:nodes graph)
                 #'< :key (lambda (node)
                            (length (graph:neighbors graph node))))))
    rings))

(defun ring-edges (ring)
  (loop for (a b) on ring
     when b
     append (list (cons a b)
                  (cons b a))))

(defun minimal-rings (rings)
  "Takes a set of rings from (say) Hanser ring perception and returns
only the minimal rings."
  (declare (optimize (debug 2)))
  (loop for (ring . rest) on (sort (copy-list rings) #'< :key #'length)
     do
     (loop for other-ring in rest
        do
        (when (null (set-difference ring other-ring))
          (return-from minimal-rings (minimal-rings
                                      (remove other-ring rings :test 'equal))))))

  rings)

(defun remove-rings-with-common-edges (rings)
  (let (edges)
    (prog1
        (loop for ring in (sort (copy-list rings) #'< :key #'length)
           unless (every #'identity 
                         (mapcar (lambda (edge)
                                   (member edge edges :test 'equal))
                                 (ring-edges ring)))
           collect ring
           do (mapcar (lambda (edge)
                        (pushnew edge edges :test 'equal))
                      (ring-edges ring))))))

(defun get-cons-node-names (ring)
  (mapcar (lambda (edge)
            (cons (atom-name (car edge))
                  (atom-name (cdr edge))))
          (ring-edges ring)))

(defun collapse-rings (graph)
  (let ((ring-graph (graph:copy-graph graph))
        (node-hash (graph::make-node-hash-table graph))
        (rings (remove-rings-with-common-edges
                (minimal-rings (chem::hanser-rings graph)))))
    (graph:map-nodes (lambda (node)
                       (setf (gethash node node-hash) node))
                     ring-graph)
    (loop for (ring . rest) on (sort (copy-list rings) #'< :key #'length)
       do 
       (let ((ring-node (make-instance 'ring-node :nodes ring)))
         (graph:add-node ring-graph ring-node)
         (map nil
              (lambda (node)
                (map nil
                     (lambda (edge)
                       (let ((neighbor (graph:other-edge-node edge node)))
                         (unless (member neighbor ring)
                           (graph:add-edge-between-nodes ring-graph
                                                         ring-node
                                                         neighbor)))
                       (graph:remove-edge ring-graph edge))
                     (graph:find-edges-containing ring-graph node))
                (graph:remove-node ring-graph node))
              ring)))
    ring-graph))
