
(asdf:defsystem #:chemicl
  :name "chemicl"
  :author "Cyrus Harmon <ch-lisp@bobobeach.com>"
  :version #.(with-open-file
                 (vers (merge-pathnames "version.lisp-expr" *load-truename*))
               (read vers))
  :licence "BSD"
  :description "A library for representing chemical structures"
  :depends-on (cxml cxml-stp)
  :components
  ((:static-file "version" :pathname #p"version.lisp-expr")
   (:cl-source-file "package")
   (:cl-source-file "graph" :depends-on (package))
   (:cl-source-file "chemicl" :depends-on (package graph))
   (:static-file "COPYRIGHT")
   (:static-file "README")
   (:static-file "make-dist" :pathname #.(make-pathname :name "make-dist" :type "sh"))))

