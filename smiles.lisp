;;; file: chemicl.lisp
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

(defun add-hydrogens (molecule &optional exclude-list)
  "Adds hydrogens to atoms in a molecule such that each atom has the
lowest normal valence consistent with the number of pre-existing bonds
for that atom. Returns two values, molecule and a list of the newly
added hydrogen ATOMs."
  (let ((hydrogen-count (count-element molecule "H"))
        atoms)
    (graph:map-nodes
     molecule
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
                   (add-bond molecule atom h-atom)))))))))
    (values molecule (nreverse atoms))))

(defparameter *aromatic-atoms* '("c" "n" "o" "p" "s" "as" "se" "*"))

(defun aromaticp (string)
  (find string *aromatic-atoms* :test 'equal))

(defun parse-smiles-string (string &optional name)
  "Parses a SMILES description of a molecule and returns an instance
of the MOLECULE class with the appropriate atoms and bonds."
  (declare (optimize (debug 3)))
  (let ((mol (apply #'make-instance 'molecule
                    (when name `(:name ,name))))
        (ring-openings (make-hash-table))
        (element-count (make-hash-table))
        first-atom
        explicit-atoms
        aromatic
        chirality)
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
                    (when (eql next-char charge-char)
                      (loop for char = (read-char stream)
                         do (incf charge)
                         while (eql charge-char (peek-char nil stream))))
                    charge))))
             (read-bracket-expression (stream)
               (let (mass)
                 (let ((char (peek-char nil stream)))
                   (when (digit-char-p char)
                     (setf mass (read-number stream)))
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
                     (when (aromaticp element-string)
                       (setf aromatic t))
                     (when mass (setf (isotope-mass atom) mass))
                     (let ((char (peek-char nil stream)))
                       (cond ((eq char #\]))
                             ((char-equal char #\H)
                              (read-char stream)
                              (let ((count
                                     (if (digit-char-p
                                          (peek-char nil stream)) 
                                         (read-number stream) 
                                         1)))
                                (dotimes (i count)
                                  (let ((h (add-molecule-atom (get-element "H"))))
                                    (add-bond mol atom h)))))))
                     (let ((char (peek-char nil stream)))
                       (cond ((eq char #\]) (read-char stream))
                             ((eq char #\-)
                              (read-char stream)
                              (let ((charge (- (read-charge stream #\-))))
                                (setf (charge atom) charge))
                              (read-char stream))
                             ((eq char #\+)
                              (read-char stream)
                              (let ((charge (read-charge stream #\+)))
                                (setf (charge atom) charge))
                              (read-char stream))))
                     (list (cons :explicit-atom atom))))))
             (add-molecule-atom (element)
               (let ((count (setf (gethash element element-count)
                                  (1+ (or (gethash element element-count) 0)))))
                 (let ((atom (add-atom mol (id element)
                                       (format nil "~A~A" (id element) count))))
                   ;; we need to squirrel away the first atom because
                   ;; the chirality w.r.t. implicit H atoms is
                   ;; different for the first ato.
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
                                (list (cons :ring lookup))
                                (list (cons :ring number))))
                          (error 'smiles-error :description "Couldn't read number!"))))
                   (char
                    (when (aromaticp (coerce (list char) 'string))
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
                                             :order (get-bond-order-number bond-type))
                                   (setf (gethash ring ring-openings) nil)))))
                          ((:atom :explicit-atom)
                           (let ((atom (cdr token))
                                 bond-direction)
                             (cond ((eql bond-type :up)
                                    (setf bond-type :single)
                                    (setf bond-direction :up))
                                   ((eql bond-type :down)
                                    (setf bond-type :single)
                                    (setf bond-direction :down)))
                             (when last (apply #'add-bond mol last atom
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

    ;; 1. Add implicit hydrogens
    (let ((implicit-h-atoms (nth-value 1 (add-hydrogens mol explicit-atoms)))))
    
    mol))

