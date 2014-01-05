;;; file: canonicalize.lisp
;;;
;;; Copyright (c) 2008-2014 Cyrus Harmon (ch-lisp@bobobeach.com)
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

;;; to compute the canonical SMILES we're going to need to do a few
;;; things:
;;; 1. compute the invariants for each atom in the molecule
;;; 2. assign a rank order to each atom
;;; 3. convert the rank into the nth prime
;;; 4. compute the product of the neighboring primes
;;; 5. rank the product of the primes using the previous ranks to
;;;    break ties

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
                     (remove-atoms-from-list "H" (graph:neighbors molecule atom)))))
   atoms))

(defun next-ranks (molecule atoms ranks)
  (let ((primes (mapcar #'nth-prime ranks)))
    (let ((product-of-primes
           (compute-product-of-primes primes atoms molecule)))
      (rank-order
       (mapcar #'list ranks product-of-primes)
       #'list<))))

(defun get-non-h-atoms (molecule &key (start (graph:first-node molecule)))
  (let (l) 
    (graph:dfs-map molecule (get-atom molecule start)
                   (lambda (x)
                     (push x l)))
    (nreverse (remove-atoms-from-list "H" l))))

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
