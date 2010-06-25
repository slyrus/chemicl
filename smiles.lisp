;;; file: smiles.lisp
;;;
;;; Copyright (c) 2008 Cyrus Harmon (ch-lisp@bobobeach.com)
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms prohibited.
;;;

(in-package :chemicl)

(define-condition smiles-error (error)
  ((description :initarg :description 
                :reader smiles-error-description))
  (:report
   (lambda (condition stream)
     (write-string (smiles-error-description condition) stream))))


(defun available-pi-electrons (molecule atom ring)
  (let* ((bonds (graph:find-edges-containing molecule atom))
         (ring-bonds (intersection bonds ring))
         (non-ring-bonds (set-difference bonds ring-bonds)))
    (let ((ring-bond-order (reduce #'+ (mapcar #'bond-order ring-bonds)))
          (non-ring-bond-order (reduce #'+ (mapcar #'bond-order non-ring-bonds))))
      (cond ((equal (id (element atom)) "C")
             (if (and (<= non-ring-bond-order 1)
                      (<= ring-bond-order 3))
                 1
                 0))
            ((equal (id (element atom)) "N")
             (if (and (<= non-ring-bond-order 1)
                      (<= ring-bond-order 3))
                 1
                 0))))))

(defun add-hydrogens (molecule &optional exclude-list)
  "Adds hydrogens to atoms in a molecule such that each atom has the
lowest normal valence consistent with the number of pre-existing bonds
for that atom. Returns two values, molecule and a list of the newly
added hydrogen ATOMs."
  (let ((hydrogen-count (count-element molecule "H"))
        atoms)
    (graph:map-nodes
     (lambda (atom)
       (unless (member atom exclude-list)
         (let* ((order (atom-bond-order molecule atom))
                (normal-valence 
                 (let ((normal-valence-list (get-normal-valence atom)))
                   (find-if (lambda (x)
                              (>= x order))
                            normal-valence-list))))
           (when normal-valence
             (let ((count (- normal-valence order)))
               (dotimes (i (truncate count))
                 (let ((h-atom (add-atom molecule 1 
                                         (format nil "H~A" (incf hydrogen-count)))))
                   (push h-atom atoms)
                   (add-bond molecule atom h-atom))))))))
     molecule)
    (values molecule (nreverse atoms))))

(defparameter *aromatic-atoms* '("c" "n" "o" "p" "s" "as" "se"))
(defparameter *organic-atoms* '("B" "C" "N" "O" "P" "S" "F" "Cl" "Br" "I"))

(defun aromatic-smiles-atom (string)
  (find string *aromatic-atoms* :test 'equal))

(defun parse-smiles-string (string &key name (add-implicit-hydrogens t))
  "Parses a SMILES description of a molecule and returns an instance
of the MOLECULE class with the appropriate atoms and bonds."
  (let ((mol (apply #'make-instance 'molecule
                    (when name `(:name ,name))))
        (ring-openings (make-hash-table))
        (element-count (make-hash-table))
        first-atom
        explicit-atoms
        aromatic)
    (labels ((read-number (stream)
               (parse-integer
                (coerce (loop for digit = (read-char stream)
                           for next = (peek-char nil stream nil)
                           collect digit
                           while (and next (digit-char-p next)))
                        'string)))
             (read-charge (stream charge-char)
               (let ((next-char (peek-char nil stream)))
                 (or
                  (let ((digit (digit-char-p next-char)))
                    (when digit (read-number stream)))
                  (let ((charge 1))
                    (when (and (not (char-equal next-char #\]))
                               (eql next-char charge-char))
                      (loop for char = (read-char stream)
                         do (incf charge)
                         while (eql charge-char (peek-char nil stream))))
                    charge))))
             (read-bracket-expression (stream)
               (let (isotope-number)
                  (let ((char (peek-char nil stream)))
                   (when (digit-char-p char)
                     (setf isotope-number (read-number stream)))
                   (let* ((element-string 
                           (coerce
                            (cons (read-char stream)
                                  (let ((next (peek-char nil stream)))
                                    (when (and (alpha-char-p next)
                                               (lower-case-p next))
                                      (cons (read-char stream) nil))))
                            'string))
                          
                          (atom (add-molecule-atom
                                (get-element element-string))))
                     (when (aromatic-smiles-atom element-string)
                       (setf aromatic t))
                     (when isotope-number
                       (let ((isotope (get-isotope (element atom) isotope-number)))
                         (setf (isotope atom) isotope)))
                     (loop for char = (peek-char nil stream)
                        do
                          (cond ((eq char #\])
                                 (read-char stream)
                                 (loop-finish))
                                ((char-equal char #\H)
                                 (read-char stream)
                                 (let ((count
                                        (if (digit-char-p
                                             (peek-char nil stream)) 
                                            (read-number stream) 
                                            1)))
                                   (dotimes (i count)
                                     (let ((h (add-molecule-atom (get-element "H"))))
                                       (add-bond mol atom h)))))
                                ((char-equal char #\-)
                                 (read-char stream)
                                 (let ((charge (- (read-charge stream #\-))))
                                   (setf (charge atom) charge)))
                                ((char-equal char #\+)
                                 (read-char stream)
                                 (let ((charge (read-charge stream #\+)))
                                   (setf (charge atom) charge)))
                                
                                ;; FIXME! Add support for @ and @@ here!!!
                                ((char-equal char #\@)
                                 (read-char stream)
                                 (if (char-equal (peek-char nil stream) #\@)
                                     ;; it's clockwise
                                     (progn
                                       (read-char stream)
                                       (push (make-instance 'tetrahedral-center
                                                            :chiral-atom atom
                                                            :orientation :clockwise)
                                             (spatial-arrangements mol)))
                                     ;; it's anticlockwise
                                     (push (make-instance 'tetrahedral-center
                                                          :chiral-atom atom
                                                          :orientation :anticlockwise)
                                           (spatial-arrangements mol))))
                                (t (read-char stream)
                                   (print "unknown SMILES character!"))))
                     (list (cons :explicit-atom atom))))))
             (add-molecule-atom (element)
               (let ((count (setf (gethash element element-count)
                                  (1+ (or (gethash element element-count) 0)))))
                 (let ((atom (add-atom mol (id element)
                                       (format nil "~A~A" (id element) count))))
                   ;; we need to squirrel away the first atom because
                   ;; the chirality w.r.t. implicit H atoms is
                   ;; different for the first atom.
                   (unless first-atom
                     (setf first-atom atom))
                   atom)))
             (read-branch (stream &optional source)
               (list (cons :branch
                           (loop for token = (read-smiles-stream stream source)
                              while token))))
             (read-smiles-tokens (stream &optional source)
               (let ((char (read-char stream nil nil)))
                 (cond
                   ((null char) nil) 
                   ((eql char #\[)
                    (read-bracket-expression stream))
                   ((eql char #\()
                    (read-branch stream source))
                   ((eql char #\)) nil)
                   ((eql char #\-) (list (cons :bond :single)))
                   ((eql char #\=) (list (cons :bond :double)))
                   ((eql char #\#) (list (cons :bond :triple)))
                   ((eql char #\:) (list (cons :bond :aromatic)))
                   ((eql char #\/) (list (cons :bond :up)))
                   ((eql char #\\) (list (cons :bond :down)))
                   ((eql char #\.) (list (cons :disconnected nil)))
                   ((or (digit-char-p char)
                        (eql char #\%))
                    (let ((number (or (digit-char-p char)
                                      (read-number stream))))

                      (if number
                          (let* ((lookup (gethash number ring-openings)))
                            (if lookup
                                (progn
                                  (setf (gethash number ring-openings) nil)
                                  (list (cons :ring lookup)))
                                (list (cons :ring number))))
                          (error 'smiles-error
                                 :description "Couldn't read number!"))))

                   ;; Explicitly handle Cl and Br as they are the only
                   ;; multi-character atoms allowed in the SMILES
                   ;; organic subset. Other atoms have to be specified
                   ;; in brackets.
                   ((and (eql char #\C)
                         (eql (peek-char nil stream nil) #\l))
                    (read-char stream)
                    (list
                     (cons :atom (add-molecule-atom (get-element "Cl")))))
                   ((and (eql char #\B)
                         (eql (peek-char nil stream nil) #\r))
                    (read-char stream)
                    (list
                     (cons :atom (add-molecule-atom (get-element "Br")))))

                   (char
                    (unless (member (string char) *organic-atoms*
                                    :test #'string-equal) 
                      (error 'smiles-error
                             :description "Unexpected element symbol"))
                    (when (aromatic-smiles-atom (coerce (list char) 'string))
                      (setf aromatic t))
                    (list (cons :atom (add-molecule-atom
                                       (get-element
                                        ;; check for *, and use
                                        ;; element 0 for wildcard.
                                        (if (char-equal char #\*)
                                            0
                                            (string char))))))))))
             (read-smiles-stream (stream &optional last)
               (loop for tokens = (read-smiles-tokens stream last)
                  with bond-type = :single
                  while tokens
                  do
                  (loop for token in tokens
                     do (case (car token) 
                          (:disconnected (setf last nil))
                          (:bond (setf bond-type (cdr token)))
                          (:ring
                           (let ((ring (cdr token)))
                             (if (numberp ring)
                                 (setf (gethash ring ring-openings) last)
                                 (progn
                                   (add-bond mol ring last
                                             :type (get-bond-order-keyword bond-type)
                                             :order (get-bond-order-number bond-type))))))
                          ((:atom :explicit-atom)
                           (let ((atom (cdr token))
                                 bond-direction)
                             (let ((last-arrangement
                                    (find last (spatial-arrangements mol) :key 'chiral-atom)))
                               (when last-arrangement
                                 (vector-push atom (neighbors last-arrangement))))
                             (let ((arrangement
                                    (find atom (spatial-arrangements mol) :key 'chiral-atom)))
                               (when arrangement
                                 (vector-push last (neighbors arrangement))))
                             (cond ((eql bond-type :up)
                                    (setf bond-type :single)
                                    (setf bond-direction :up))
                                   ((eql bond-type :down)
                                    (setf bond-type :single)
                                    (setf bond-direction :down)))
                             (when last
                               (apply #'add-bond mol last atom
                                      :type (get-bond-order-keyword bond-type)
                                      :order (get-bond-order-number bond-type)
                                      (when bond-direction
                                        `(:direction ,bond-direction))))
                             (setf bond-type (if aromatic
                                                 :aromatic 
                                                 :single))
                             (setf last atom)
                             (when (eql (car token) :explicit-atom)
                               (push atom explicit-atoms)))))))))
      (with-input-from-string (stream string)
        (read-smiles-stream stream)))
    ;; Now we need to do some post-processing.

    ;; 2. Find aromatic rings and replace Kekule structures with
    ;; explicitly aromatic rings. Unfortunately, we need to do this
    ;; BEFORE we add the hydrogens, I think, 
    
    #+nil 
    (multiple-value-bind (bonds cycles cycles-removed-mol)
        (graph:break-cycles mol)
      (print bonds)
      (print cycles)
      
      )
    ;; 1. Add implicit hydrogens
    (when add-implicit-hydrogens
      (let ((implicit-h-atoms (nth-value 1 (add-hydrogens mol explicit-atoms))))))
      
    ;; 3. Add double-bond-configurations for explicit spatial
    ;; arrangements around double bonds, based on the information
    ;; placed in the bonds when we were parsing them above.
    
    (let ((explicit-double-bonds
           (mapcan (lambda (x)
                     (let ((left-neighbors
                            (remove-if-not
                             #'bond-direction
                             (find-bonds-containing mol (atom1 x))))
                           (right-neighbors
                            (remove-if-not
                             #'bond-direction
                             (find-bonds-containing mol (atom2 x)))))
                       (when (and left-neighbors right-neighbors)
                         (list (list x left-neighbors right-neighbors))))) 
                   (remove-if-not (lambda (x) (= (bond-order x) 2)) (bonds mol)))))
      (mapcar (lambda (x)
                (destructuring-bind (bond left right) x
                  (let ((config (make-instance 'double-bond-configuration
                                               :bond bond
                                               :left-atom (atom1 bond)
                                               :right-atom (atom2 bond))))
                    (flet ((check-and-set-substituent (index atom)
                             (if (elt (substituents config) index)
                                 (error 'smiles-error
                                        :description "Invalid double bond configuration")
                                 (setf (elt (substituents config) index) atom))))
                      (loop for child in left
                         do
                           (cond ((and (eq (atom1 bond) (atom1 child))
                                       (eq (bond-direction child) :up))
                                  (check-and-set-substituent +bottom-left+ child))
                                 ((and (eq (atom1 bond) (atom1 child))
                                       (eq (bond-direction child) :down))
                                  (check-and-set-substituent +top-left+ child))
                               ((and (eq (atom1 bond) (atom2 child))
                                     (eq (bond-direction child) :up))
                                (check-and-set-substituent +bottom-left+ child))
                               ((and (eq (atom1 bond) (atom2 child))
                                     (eq (bond-direction child) :down))
                                (check-and-set-substituent +top-left+ child))))

                      (loop for child in right
                         do
                           (cond ((and (eq (atom2 bond) (atom1 child))
                                       (eq (bond-direction child) :up))
                                  (check-and-set-substituent +top-right+ child))
                                 ((and (eq (atom2 bond) (atom1 child))
                                       (eq (bond-direction child) :down))
                                (check-and-set-substituent +bottom-right+ child))
                                 ((and (eq (atom2 bond) (atom2 child))
                                       (eq (bond-direction child) :up))
                                  (check-and-set-substituent +bottom-right+ child))
                                 ((and (eq (atom2 bond) (atom2 child))
                                       (eq (bond-direction child) :down))
                                  (check-and-set-substituent +top-right+ child)))))

                    ;; FIXME now we should fill in the missing
                    ;; (implicitly specified) bonds

                    (push config (spatial-arrangements mol)))))
              explicit-double-bonds))

    (values mol)))

(defun remove-atoms (element-identifier atom-list)
  (remove
   (get-element element-identifier)
   atom-list
   :key #'element))

(defun remove-edges-containing (element-identifier edge-list)
  (remove-if
   (lambda (x)
     (member (get-element element-identifier)
             (graph:edge-nodes x)
             :key #'element))
   edge-list))

(defun remove-edges-not-containing (element-identifier edge-list)
  (remove-if-not
   (lambda (x)
     (member (get-element element-identifier)
             (graph:edge-nodes x)
             :key #'element))
   edge-list))

;;; The paper "SMILES. 2. Algorithm for Generation of Unique SMILES
;;; Notation" allegedly contains descriptions of how to generate a
;;; canonical SMILES representation of a given atom. See that paper
;;; for details.
(defun canon-invariant (molecule atom)
  "Returns a so-called invariant for an atom in a molecule."
  (let* ((bonds (graph:find-edges-containing molecule atom))
         (non-h-bonds (remove-edges-containing "H" bonds))
         (h-bonds (remove-edges-not-containing "H" bonds)))
    (let ((invariant-list
           (list
            
            ;; (1) number of connections. This is unclear to me. Is this the
            ;; total number of connected atoms or the number of connected
            ;; non-hydrogen atoms? The text would make me think the former
            ;; but the examples suggest the latter.
            (length non-h-bonds)

            ;; (2) number of non-hydrogen bonds. Given the interpretation of
            ;; (1), how could this be different? I think they mean the total
            ;; bond-order of non-h bonds, that is a double bond would count 2
            ;; for this.
            (reduce #'+ non-h-bonds :key #'bond-order)

            ;; (3) atomic number -- this one is fairly non-controversial.
            (atomic-number atom)

            ;; (4) sign of charge.  The paper says to use the sign of the
            ;; formal charge but all of the so-called invariants are positive
            ;; numbers (that eventually get concatenated together). I would
            ;; think that positive numbers would get 1 for this, but the
            ;; examples suggest otherwise. So, give negatively charged atoms
            ;; a 1 and 0 otherwise.
            (if (minusp (charge atom)) 1 0)

            ;; (5) absolute value of the charge. again, this one seems
            ;; straightforward.
            (abs (charge atom))
     
            ;; (6) number of attached hydrogens.
            (length h-bonds))))
      invariant-list)))

;;; to compute the canonical SMILES we're going to need to do a few
;;; things:
;;; 1. compute the invariants for each atom in the molecule
;;; 2. assign a rank order to each atom
;;; 3. convert the rank into the nth prime
;;; 4. compute the product of the neighboring primes
;;; 5. rank the product of the primes using the previous ranks to
;;;    break ties

(defun get-non-h-atoms (molecule &key (start (graph:first-node molecule)))
  (let (l) 
    (graph:dfs-map molecule (get-atom molecule start)
                   (lambda (x)
                     (push x l)))
    (nreverse (remove-atoms "H" l))))

(defun canon-invariants (molecule atoms)
  (mapcar (lambda (x) 
            (canon-invariant molecule x))
          atoms))

(defun rank-order (list predicate &key (key 'identity) (test 'equal))
  "Returns the rank-order of the items in list sorted by predicate
using the key "
  (let (unique)
    (map nil (lambda (x) (pushnew (funcall key x)
                                  unique :test test))
         list)
    (let ((sorted (sort unique predicate)))
      (mapcar (lambda (x)
                (position (funcall key x) sorted
                          :test test))
              list))))

(defun assign-primes (ranked-invariants)
  (mapcar
   (lambda (x)
     (let ((prime (nth-prime (1- (car x)))))
       (cons prime (cdr x))))
   ranked-invariants))

(defun compute-product-of-primes (primes atoms molecule)
  (mapcar
   (lambda (atom)
     (reduce #'*
             (mapcar (lambda (y)
                       (elt primes (position y atoms)))
                     (remove-atoms "H" (graph:neighbors molecule atom)))))
   atoms))

(defun next-ranks (molecule atoms ranks)
  (let ((primes (mapcar #'nth-prime ranks)))
    (let ((product-of-primes
           (compute-product-of-primes primes atoms molecule)))
      (rank-order
       (mapcar #'list ranks product-of-primes)
       #'list<))))

(defun canonicalize-atoms-1 (molecule &optional atoms ranks)
  (let* ((atoms (or atoms (get-non-h-atoms molecule)))
         (invariants (canon-invariants molecule atoms))
         (ranks (or ranks (rank-order invariants #'list<))))
    (list
     (loop with last-ranks = ranks
        for next-ranks = (next-ranks molecule atoms last-ranks)
        for i below 500
        while (not (equal next-ranks last-ranks))
        do (setf last-ranks next-ranks)
        finally (return (mapcar #'1+ next-ranks)))
     atoms)))

(defun canonicalize-atoms (molecule)
  (loop  for (ranks atoms) = (canonicalize-atoms-1 molecule atoms ranks)
     for dup = (find-duplicate ranks)
     while dup
     do (let ((pos (position dup ranks)))
          (setf ranks (mapcar #'(lambda (x) (ash x 1)) ranks))
          (decf (elt ranks pos)))
     finally (return (list ranks atoms))))

(defun remove-implicit-h-atoms (list)
  (remove-atoms "H" list))

(defun write-smiles-string-to-stream (molecule stream &key)
  (declare (optimize (debug 2)))
  (destructuring-bind (ranks atoms)
      (canonicalize-atoms molecule)
    (let ((visited-atoms (make-hash-table))
          (start (elt atoms (position 1 ranks)))
          (rank-hash (make-hash-table))
          (cycle-counter 1)
          (cycle-hash (make-hash-table)))
      (map nil (lambda (rank atom)
                 (setf (gethash atom rank-hash) rank))
           ranks atoms)
      (multiple-value-bind (cycle-edge-lists
                            cycle-node-lists
                            removed-edge-list
                            broken-molecule) 
          (graph:break-cycles molecule
           :pick-function
           (lambda (x y)
             (cond
               ((< (bond-order x) (bond-order y)) x)
               ((> (bond-order x) (bond-order y)) y)
               (t (let ((base (car (intersection
                                    (graph:edge-nodes x)
                                    (graph:edge-nodes y)))))
                    (if (< (gethash (graph:other-edge-node x base)
                                    rank-hash)
                           (gethash (graph:other-edge-node y base)
                                    rank-hash))
                        y
                        x))))))
        (labels
            ((print-atom (atom)
               (princ (id atom) stream))
             (print-bond (bond)
               (case (bond-order bond)
                 (2 (princ "=" stream))
                 (3 (princ "#" stream))))
             (dfs-visit (atom path)
               (setf (gethash atom visited-atoms) atom)
               (let* ((neighbors
                       (remove-if
                        (lambda (x) (gethash x visited-atoms))
                        (remove-implicit-h-atoms
                         (graph:neighbors
                          (or broken-molecule molecule) atom))))
                      (count (length neighbors))
                      (i 0)
                      (cycle-edges
                       (remove-if-not
                        (lambda (edge)
                          (or (equal atom (graph:node1 edge))
                              (equal atom (graph:node2 edge))))
                        removed-edge-list)))
                 (when path
                   (let ((bond (graph:edgep molecule (car path) atom)))
                     (print-bond bond)))
                 (print-atom atom)
                 (let (cycle-openings cycle-closings)
                   (map nil
                        (lambda (cycle-edge)
                          (let ((lookup (gethash cycle-edge cycle-hash)))
                            (if lookup
                                (push (cons lookup cycle-edge) cycle-closings)
                                (progn
                                  (push (cons cycle-counter cycle-edge) cycle-openings)
                                  (setf (gethash cycle-edge cycle-hash)
                                        cycle-counter)
                                  (incf cycle-counter)))))
                        cycle-edges)
                   (map nil
                        (lambda (list)
                          (destructuring-bind (cycle-label . cycle-edge)
                              list
                            (princ cycle-label stream)))
                        (sort cycle-openings #'< :key #'car))
                   (map nil
                        (lambda (list)
                          (destructuring-bind (cycle-label . cycle-edge)
                              list
                            (print-bond cycle-edge)
                            (princ cycle-label stream)))
                        (sort cycle-closings #'< :key #'car)))
                 (map nil
                      (lambda (x)
                        (when (and (> count 1) (< i (1- count)))
                          (princ "(" stream))
                        (dfs-visit x (cons atom path))
                        (when (and (> count 1) (< i (1- count)))
                          (princ ")" stream))
                        (incf i))
                      (sort
                       (remove-if
                        (lambda (x) (gethash x visited-atoms))
                        neighbors)
                       #'list<
                       :key
                       (lambda (x)
                         (list
                          (let ((bond (graph:edgep molecule x atom)))
                            (- (if (member bond
                                           (apply #'append cycle-edge-lists)
                                           :test 'graph:edges-nodes-equal)
                                   (bond-order bond)
                                   1)))
                          (gethash x rank-hash))))))))
          (dfs-visit start nil))))))

(defun write-smiles-string (molecule)
  (with-output-to-string (stream)
    (write-smiles-string-to-stream molecule stream)))
