
(in-package :chemicl)

(defun list< (&rest lists)
  (let (done ret)
    (apply #'mapcar
           (lambda (&rest vals)
             (unless (or done (apply #'= vals))
               (if (apply #'< vals)
                   (setf done t ret t)
                   (setf done t ret nil))))
           lists)
    ret))

(defun list> (&rest lists)
  (let (done ret)
    (apply #'mapcar
           (lambda (&rest vals)
             (unless (or done (apply #'= vals))
               (if (apply #'> vals)
                   (setf done t ret t)
                   (setf done t ret nil))))
           lists)
    ret))

(defun find-duplicate (sequence &key (test #'eql))
  (let ((hash (make-hash-table :test test)))
    (map nil (lambda (x)
               (if (gethash x hash)
                   (return-from find-duplicate x)
                   (setf (gethash x hash) x)))
         sequence)))

