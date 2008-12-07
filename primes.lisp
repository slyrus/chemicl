;;; file: primes.lisp
;;;
;;; Copyright (c) 2008 Cyrus Harmon (ch-lisp@bobobeach.com)
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms prohibited.
;;;

(in-package :chemicl)

;; we could certainly just store these in a file, but it's more fun to
;; generate it dynamically.
(defun primes-up-to-n (n)
  (let ((not-primes (make-array n
                                :element-type '(unsigned-byte 1)
                                :initial-element 0
                                :adjustable nil))
        (sqrtn (truncate (sqrt n))))
    (loop for i from 2 to sqrtn
       do
         (unless (not (zerop (aref not-primes i)))
           (loop for j from (* 2 i) below n by i
              do
                (progn
                  (setf (aref not-primes j) 1)))))
    (loop for i from 2 below n
       when (zerop (aref not-primes i))
       collect i)))

(defconstant +nth-prime-limit+ 2500)
(defconstant +nth-prime+ 22343)

(let ((primes (primes-up-to-n +nth-prime+)))
  (defparameter *nth-prime-vector*
    (make-array (length primes)
                :element-type '(unsigned-byte 32)
                :initial-contents primes)))

(defun nth-prime (n)
  (when (>= n (length *nth-prime-vector*))
    (error "Not enough primes!"))
  (aref *nth-prime-vector* n))
