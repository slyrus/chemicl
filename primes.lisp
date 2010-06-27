;;; file: primes.lisp
;;;
;;; Copyright (c) 2008 Cyrus Harmon (ch-lisp@bobobeach.com)
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

(eval-when (:compile-toplevel :load-toplevel :execute)
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
                  :initial-contents primes))))

(defun nth-prime (n)
  (when (>= n (length *nth-prime-vector*))
    (error "Not enough primes!"))
  (aref *nth-prime-vector* n))
