
(asdf:defsystem :chemicl-test
  :name "chemicl-test"
  :author "Cyrus Harmon <ch-lisp@bobobeach.com>"
  :version #.(with-open-file
                 (vers (merge-pathnames "version.lisp-expr" *load-truename*))
               (read vers))
  :licence "For Internal Use Only"
  :depends-on (chemicl)
  :components
  ((:module
    :test
    :components
    ((:cl-source-file "package")
     (:cl-source-file "chemicl-test" :depends-on (package))
     (:cl-source-file "smiles-test" :depends-on (package))))))
