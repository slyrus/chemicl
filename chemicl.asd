
(asdf:defsystem #:chemicl
  :name "chemicl"
  :author "Cyrus Harmon <ch-lisp@bobobeach.com>"
  :version #.(with-open-file
                 (vers (merge-pathnames "version.lisp-expr" *load-truename*))
               (read vers))
  :licence "For Internal Use Only"
  :description "A library for representing chemical structures"
  :depends-on (cxml cxml-stp epigraph)
  :components
  ((:static-file "version" :pathname #p"version.lisp-expr")
   (:cl-source-file "package")
   (:cl-source-file "elements" :depends-on (package))
   (:cl-source-file "chemicl" :depends-on (package elements))
   (:cl-source-file "primes" :depends-on (package chemicl))
   (:cl-source-file "smiles" :depends-on (package chemicl primes))
   (:static-file "COPYRIGHT")
   
   (:static-file "README")
   (:static-file "make-dist" :pathname #.(make-pathname :name "make-dist" :type "sh"))
   (:module "data"
            :components ((:static-file "elementdata.xml")
                         (:static-file "isotopes.xml")))))

