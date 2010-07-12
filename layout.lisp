
(in-package :chemicl)

(defparameter *benzene* (parse-smiles-string "c1ccccc1"
                                             :name "benzene"))

(defparameter *cyclooctane* (parse-smiles-string "C1CCCCCCC1"
                                                 :name "cyclooctane"))

(defun layout-cycle (nodes &key
                     (ring-center #(0.0d0 0.0d0))
                     (start-angle 0.0d0)
                     (radius 1.0d0))
  (let ((atom-count (length nodes)))
    (let ((phi (/ (* 2 pi) atom-count)))
      (make-array atom-count
                  :initial-contents
                  (loop for i below atom-count
                     for theta from start-angle by phi
                     collect
                       (vector (+ (aref ring-center 0) (* radius (cos theta)))
                               (+ (aref ring-center 1) (* radius (sin theta)))))))))

(let ((mol *cyclooctane*))
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
          ;; let's lay out the first ring manually
          (let ((first-ring (car (graph:nodes ring-graph))))
            first-ring
            (let ((nodes (butlast (graph:nodes first-ring))))
              (layout-cycle nodes))))))))
